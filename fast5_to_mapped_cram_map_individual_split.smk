# This version splits the guppy runs.

import os.path
wildcard_constraints:
    type="pass|fail",
    f5type="pass|fail",
    fast5name="[a-zA-Z0-9_]+",
    suffix="|[.]gz",
    split="[0123]"

report: "run_report/methylmapping.rst"

rule all:
    input:
        config["sampleid"]+".pass.minimap2.sort.cram",
        config["sampleid"]+".fail.minimap2.sort.cram"
    shell:  # This is fairly dangerous step. Might loose several days of work.
        "echo rm [0123]/guppy_output/workspace/*.fast5"




rule copy_fast5:
    input:
        fast5_file=os.path.join(config["fast5_path"],"fast5_{type}/{fast5name}.fast5")
    output:
        local_file=temp("{split}/input_fast5_copy/{type}_{fast5name}.fast5")
    shell:
        "cp {input.fast5_file} {output.local_file};"
        "chmod a+rw {output.local_file}"

rule uncompress_fast5:
    input:
        fast5_file=os.path.join(config["fast5_path"],"fast5_{type}/{fast5name}.fast5.gz")
    output:
        local_file=temp("{split}/input_fast5_copy/{type}_{fast5name}.gz.fast5")
    shell:
        "zcat {input.fast5_file} > {output.local_file};"
        "chmod a+rw {output.local_file}"

rule link_fast5:
    input:
        fast5_file="input_fast5_copy/{new_fast5_name}.fast5"
    output:
        split_file="{split}/input_fast5_copy/{new_fast5_name}.fast5"
    shell:
        "ln {input.fast5_file} {output.split_file}"

ruleorder: link_fast5 > copy_fast5 > uncompress_fast5


def split_fast5s(wildcards):
    import os.path
    fast5_path = os.path.join(config["fast5_path"],"fast5_{type}/{fast5name}.fast5{suffix,|.gz}")
    
    wcards = glob_wildcards(fast5_path)
    
    split_idx = int(wildcards.split)
    s = slice(split_idx,None,4)
    import itertools as it
    
    all_files =  expand("{split}/input_fast5_copy/{type}_{fast5name}{suffix}.fast5",zip,
            split=it.repeat(wildcards.split,len(wcards.type[s])),
            type=wcards.type[s],
            fast5name=wcards.fast5name[s],
            suffix=wcards.suffix[s])
    

    return all_files

checkpoint call_methylation:
    input:
        input_fast5s=split_fast5s
    output:
        #basecalled=directory("guppy_output/"),
        #basecalled_f5=temp(directory("{split}/guppy_output/workspace/")),
        flag=touch("{split,\d}/guppy_output/guppy.done")
    threads: 4
    benchmark: "{split}/guppy_output/guppy.time.txt"
    shell:
        """
        RESUME=""
        if [ -e {wildcards.split}/guppy_output/sequencing_summary.txt ]; then
            RESUME="--resume"
        fi
        
        flock $(which guppy_basecaller).{wildcards.split} /usr/bin/time -v guppy_basecaller $RESUME --gpu_runners_per_device  16 --config {config[guppy_config]} -x 'cuda:{wildcards.split}:100%'  \
        --records_per_fastq 100000 --save_path {wildcards.split}/guppy_output/  --compress_fastq --fast5_out -i {wildcards.split}/input_fast5_copy/
        """






def get_aligned_files(wildcards):
    for s in ["0","1","2","3"]:
        checkpoint_output = checkpoints.call_methylation.get(split=s).output[0]
    
    w = glob_wildcards("{split}/guppy_output/workspace/{new_fast5_name}.fast5")    
    
    fs = expand("mapped/{new_fast5_name}.{{type}}.minimap2.sort.bam",new_fast5_name=w.new_fast5_name)
    return fs

def get_header(wildcards):
    import os.path


    checkpoint_output = checkpoints.call_methylation.get(split="0").output[0]
    from snakemake.exceptions import IncompleteCheckpointException

    out_path = "0/guppy_output/workspace/"
    w = glob_wildcards(os.path.join(out_path+"{new_fast5_name}.fast5"))
    return "unmapped/{new_fast5_name}.sam.header".format(new_fast5_name=w.new_fast5_name[0])



def get_fast5(wildcards):


    #checkpoint_output = checkpoints.call_methylation.get(**wildcards).output[0]
    import glob
    import os.path
    files = list(glob.glob("[0123]/guppy_output/workspace/{new_fast5_name}.fast5".format(new_fast5_name=wildcards.new_fast5_name)))
    splits = list(set(x.split("/")[0] for x in files))
    return {"fast5s":files,
        "guppy_done":[os.path.join(x,"/guppy_output/guppy.done") for x in splits]}

rule extract_methylation_likelihood_filter_fastq:
    input:
        unpack(get_fast5)
        #guppy_done="a?/guppy_output/guppy.done" # This might be good, or not
    output:
        fastq_pass=pipe("unmapped/{new_fast5_name}.pass_fastq"),
        fastq_fail=pipe("unmapped/{new_fast5_name}.fail_fastq"),
        header="unmapped/{new_fast5_name}.sam.header"
    shell:
        "python $(dirname {workflow.snakefile})/scripts/extract_methylation_fast5_to_sam.py --fastq {output.header} "
        "-o {output.fastq_pass} --failed_reads {output.fastq_fail} -V -L -F -- {input.fast5s}"



rule map_fastq:
    input:
        unmapped="unmapped/{new_fast5_name}.{type}_fastq",
        reference=config["reference_fasta"]
    output:
        cram=temp("mapped/{new_fast5_name}.{type,pass|fail}.minimap2.sort.bam")
    threads: 4
    shell:
        "minimap2 -x map-ont -y -a -t 2  {input.reference} {input.unmapped} |"
        "samtools sort -O bam -l 0  -@ 2 -m 15G --reference {input.reference}  -o {output.cram} /dev/stdin"


rule merge_bams:
    input:
        mapped=get_aligned_files,
        guppy_done = expand("{split}/guppy_output/guppy.done",split=list("0123")),
        reference=config["reference_fasta"]
    output:
        cram=temp(config["sampleid"]+".{type}.minimap2.sort.norg.bam")
    threads:
        10
    shell:
        "samtools merge --threads {threads} {output.cram} {input.mapped} "
        

rule add_rg:
    input:
        mapped=rules.merge_bams.output.cram,
        reference=config["reference_fasta"],
        header=get_header
    output:
        cram=protected(config["sampleid"]+".{type}.minimap2.sort.cram")
    threads:
        10
    shell:
        "RGID=$(sed -e '/^@RG/!d' {input.header});"
        "samtools addreplacerg -r \"${{RGID}}\" -o {output.cram} -O cram --reference {input.reference} {input.mapped};"
        "samtools index {output.cram}"
