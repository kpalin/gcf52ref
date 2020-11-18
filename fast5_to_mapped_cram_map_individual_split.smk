# This version splits the guppy runs.
import logging
logging.basicConfig(level=logging.DEBUG)



nsplits = int(config.get("nsplits",4))
SPLITS=[str(i) for i in range(nsplits)]
import os.path
wildcard_constraints:
    type="pass|fail|skip",
    f5type="pass|fail",
    fast5name="[a-zA-Z0-9_-]+",
    suffix="|[.]gz",
    split="[{}]".format("".join(SPLITS))

report: "run_report/methylmapping.rst"


if "sftp_server" in config:
    from snakemake.remote.SFTP import RemoteProvider
    from paramiko import Agent
    STORAGE_SERVER=config["sftp_server"]
    #"10.110.9.8"
    SFTP = RemoteProvider(private_key=Agent().get_keys()[0])
    def input_fast5_fmt(pattern):
        if isinstance(pattern,str):
            pattern=[ pattern ]
        return SFTP.remote([config["sftp_server"]+("" if x.startswith("/") else "/")+x for x in pattern])
else:
    def input_fast5_fmt(pattern):
        return pattern


rule all:
    input:
        config["sampleid"]+".pass.minimap2.sort.stats",
        config["sampleid"]+".fail.minimap2.sort.stats",
        "sequencing_summary.txt.gz"
    run:  # This is fairly dangerous step. Might loose several days of work.
        from pathlib import Path
        for s in SPLITS:
            for f in Path(f"{s}/guppy_output/workspace/").glob("*.fast5"):
                f.unlink()

rule sequencing_summary:
    input:
        guppy_done=expand("{split}/guppy_output/guppy.done",split=SPLITS)
    output:
        sequencing_summary="sequencing_summary.txt.gz"
    run:
        shell(f"cat <(head -1q {SPLITS[0]}/guppy_output/sequencing_summary.txt) <(tail --lines=+2 -q " +
            " ".join( (x+"/guppy_output/sequencing_summary.txt")  for x in SPLITS ) + " )|" +
            "gzip >{output.sequencing_summary}")


#        shell(f"cat {SPLITS[0]}/guppy_output/sequencing_summary.txt <(tail --lines=+1 -q " +
#            " ".join( (x+"/guppy_output/sequencing_summary.txt")  for x in SPLITS[1:] ) + " )|" +
#            "gzip >{output.sequencing_summary}")




rule copy_fast5:
    input:
        fast5_file=input_fast5_fmt(os.path.join(config["fast5_path"],"fast5_{type}/{fast5name}.fast5"))
    output:
        local_file=temp("{split}/input_fast5_copy/{type}_{fast5name}.fast5")
    resources:
        disk_io=10
    shell:
        "cp {input.fast5_file} {output.local_file};"
        "chmod a+rw {output.local_file}"

rule uncompress_fast5:
    input:
        fast5_file=input_fast5_fmt(os.path.join(config["fast5_path"],"fast5_{type}/{fast5name}.fast5.gz"))
    output:
        local_file=temp("{split}/input_fast5_copy/{type}_{fast5name}.gz.fast5")
    resources:
        disk_io=10
    shell:
        "zcat {input.fast5_file} > {output.local_file};"
        "chmod a+rw {output.local_file}"

#rule link_fast5:
#    input:
#        fast5_file="input_fast5_copy/{new_fast5_name}.fast5"
#    output:
#        split_file=temp("{split}/input_fast5_copy/{new_fast5_name}.fast5")
#    shell:
#        "ln {input.fast5_file} {output.split_file}"

ruleorder: copy_fast5 > uncompress_fast5

import functools
@functools.lru_cache(10)
def get_split_wildcards(split):
    import os.path
    info=snakemake.logger.info 
    
    info(f"split {split}")
    fast5_path = os.path.join(config["fast5_path"],"{justFast5,fast5_}{type}/{fast5name}.fast5{suffix,|.gz}")
    info(str(fast5_path))
    
    if "sftp_server" in config:
        _remote_path = config["sftp_server"]+str(fast5_path)
        info(f"Remote path: {_remote_path}")
        wcards = SFTP.glob_wildcards(_remote_path)
    else:
        wcards = glob_wildcards(fast5_path)
    info(str(wcards))

    split_idx = int(split)
    s = slice(split_idx,None,len(SPLITS))
    ret = dict(type=wcards.type[s],
            fast5name=wcards.fast5name[s],
            suffix=wcards.suffix[s])
    
    assert len(ret["fast5name"])>0, fast5_path
    
    return ret

def split_fast5s(wildcards):
    import itertools as it
    
    split_files = get_split_wildcards(wildcards.split)
    all_files =  expand("{split}/input_fast5_copy/{type}_{fast5name}{suffix}.fast5",zip,
            split=it.repeat(wildcards.split,len(split_files["fast5name"])),
            **split_files)
    

    return all_files

checkpoint call_methylation:
    input:
        input_fast5s=split_fast5s
    output:
        #basecalled=directory("guppy_output/"),
        #basecalled_f5=temp(directory("{split}/guppy_output/workspace/")),
#        sequencing_summary="{split}/guppy_output/sequencing_summary.txt",
        flag=touch("{split,\d}/guppy_output/guppy.done")
    threads: 8
    benchmark: "{split}/guppy_output/guppy.time.txt"
    log:"{split}/guppy_output/guppy.snake.log"
    shell:
        """
        RESUME=""
        if [ -e {wildcards.split}/guppy_output/sequencing_summary.txt ]; then
            RESUME="--resume"
        fi
        
        flock $(which guppy_basecaller).{wildcards.split} /usr/bin/time -v guppy_basecaller $RESUME  --config {config[guppy_config]} -x "cuda:$(( {wildcards.split} % 4 )):100%"  \
        --num_callers {threads} --gpu_runners_per_device 16 --chunks_per_runner 512 \
        --records_per_fastq 100000 --save_path $(readlink -f {wildcards.split}/guppy_output/ )  --compress_fastq --fast5_out -i $(readlink -f {wildcards.split}/input_fast5_copy/ ) 2>&1 |tee {log}
        """


import functools
@functools.lru_cache(10000)
def _get_new_fast5_names(split):
    import subprocess as sp
    cmd = f"cut -f1 {split}/guppy_output/sequencing_summary.txt |uniq"
    p = sp.Popen(cmd,shell=True,stdout=sp.PIPE,universal_newlines=True)
    p.stdout.readline()
    rest = set(x.strip()[:-len(".fast5")] for x in p.stdout if x.strip().endswith(".fast5"))
    return tuple(rest)


def get_new_fast5_names(split):

    checkpoint_output = checkpoints.call_methylation.get(split=split)
    rest = _get_new_fast5_names(split)
    return rest


def get_fastqs(wildcards):

    w = get_new_fast5_names(wildcards.split) 
    fs = expand("{{split}}/unmapped/{new_fast5_name}.{{type}}_fastq.gz",new_fast5_name=w)
    return fs

def get_header(wildcards):
    fname = get_new_fast5_names("0")[0]
    w = f"0/unmapped/{fname}.sam.header"
    return w 




rule extract_methylation_likelihood_filter_fastq:
    input:
        fast5s="{split}/guppy_output/workspace/{new_fast5_name}.fast5",
        guppy_done="{split}/guppy_output/guppy.done"
    output:
        fastq_pass=temp("{split}/unmapped/{new_fast5_name}.pass_fastq.gz"),
        fastq_fail=temp("{split}/unmapped/{new_fast5_name}.fail_fastq.gz"),
        header="{split}/unmapped/{new_fast5_name}.sam.header"
    threads: 2  # Not really, at least 3 but we'll never get to use that much CPU
    log: "{split}/unmapped/{new_fast5_name}.log"
    priority: 30
    resources:
        disk_io=3
    shell:
        "python $(dirname {workflow.snakefile})/scripts/extract_methylation_fast5_to_sam.py --fastq {output.header} "
        "-o >(gzip > {output.fastq_pass} ) --failed_reads >(gzip >{output.fastq_fail}) -L -F -- {input.fast5s} 2>{log}"

rule map_index:
    input:
        reference=config["reference_fasta"]
    output:
        ref_idx=temp("tmp/minimap2.map-ont.idx")
    shell:
        "minimap2 -x map-ont -d {output.ref_idx} {input.reference} "


rule write_fastq_fofn:
    input:
        unmapped_fastqs=get_fastqs
    output:
        fofn="{split}/mapped/fastq_{type}.fofn"
    run:
        print("Should write this: "+str(input.unmapped_fastqs))
        open(output["fofn"],"w").write("\n".join(input.unmapped_fastqs))
        print("Wrote:"+str(open(output["fofn"],"rt").read()))


    

rule map_fastq:
    input:
        fastqs=get_fastqs,
        reference=config["reference_fasta"],
        ref_idx=rules.map_index.output.ref_idx
    output:
        cram=temp("{split}/mapped/mapped.{type,pass|fail}.minimap2.sort.split.cram")
    threads: 20
    benchmark:
        "{split}/mapped/mapped.{type,pass|fail}.minimap2.sort.split.time.txt"
    log: 
        sort="{split}/mapped/mapped.{type,pass|fail}.sort.log",
        map="{split}/mapped/mapped.{type,pass|fail}.minimap2.log"
    shell:
        "minimap2 -x map-ont -y -a -t {threads}  {input.ref_idx} {input.fastqs} 2>{log.map} |"
        "samtools sort -O cram -l 9  -@ 6 -m 5G --reference {input.reference}  -o {output.cram} /dev/stdin 2>{log.sort}"


rule merge_bams:
    input:
        mapped=expand("{split}/mapped/mapped.{{type}}.minimap2.sort.split.cram",split=SPLITS),
        guppy_done = expand("{split}/guppy_output/guppy.done",split=SPLITS),
        reference=config["reference_fasta"]
    output:
        cram=temp("{sample_id}.{type}.minimap2.sort.norg.bam")
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
        cram=protected("{sample_id}.{type}.minimap2.sort.cram")
    threads:
        10
    shell:
        "RGID=$(sed -e '/^@RG/!d' {input.header});"
        "samtools addreplacerg --threads {threads} -r \"${{RGID}}\" -o {output.cram} -O cram --reference {input.reference} {input.mapped};"
        "samtools index {output.cram}"




rule cram_stats:
    input:
        cram=rules.add_rg.output.cram,
        reference=config["reference_fasta"]
    output:
        stats="{sample_id}.{type}.minimap2.sort.stats"
    threads:
        15
    shell:
        "samtools stats --threads {threads} --ref-seq {input.reference} {input.cram} > {output.stats}" 