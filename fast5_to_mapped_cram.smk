#REFERENCE="/data/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta"
#import glob 
#FAST5S=list(glob.glob("/data/kpalin/Fam_c745_1_9623TK/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac/workspace/PAD64960_b924920b26f29f6118f165b040093705d5beb46b_1_fast5_*.fast5"))

import os.path
wildcard_constraints:
    type="pass|fail",
    fast5name="[^.]+"
report: "run_report/methylmapping.rst"

rule all:
    input:
        config["sampleid"]+".pass.minimap2.sort.cram",
        config["sampleid"]+".fail.minimap2.sort.cram"



#guppy_basecaller --config $HOME/ont-guppy/data/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg -x 'cuda:0 cuda:1 cuda:2 cuda:3' --records_per_fastq 100000 --save_path $PWD/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac/  --compress_fastq --fast5_out -i $PWD/workspace/

rule copy_fast5:
    input:
        fast5_file=os.path.join(config["fast5_path"],"fast5_{type}/{fast5name}.fast5")
    output:
        local_file=temp("input_fast5_copy/{type}_{fast5name}.fast5")
    shell:
        "cp {input.fast5_file} {output.local_file}"

rule uncompress_fast5:
    input:
        fast5_file=os.path.join(config["fast5_path"],"fast5_{type}/{fast5name}.fast5.gz")
    output:
        local_file=temp("input_fast5_copy/{type}_{fast5name}.gz.fast5")
    shell:
        "zcat {input.fast5_file} > {output.local_file}"




def all_fast5s(wildcards):
    import os.path
    fast5_path = os.path.join(config["fast5_path"],"fast5_{type}/{fast5name}.fast5{suffix,|.gz}")
    
    wcards = glob_wildcards(fast5_path)

    all_files =  expand("input_fast5_copy/{type}_{fast5name}{suffix}.fast5",zip,
            type=wcards.type,
            fast5name=wcards.fast5name,
            suffix=wcards.suffix)
    

    return all_files

rule call_methylation:
    input:
        input_fast5s=all_fast5s
    output:
        #basecalled=directory("guppy_output/"),
        flag=touch("guppy_output/guppy.done")
    threads: 13
    benchmark: "guppy_output/guppy.time.txt"
    shell:
        """
        RESUME=""
        if [ -e guppy_output/sequencing_summary.txt ]; then
            RESUME="--resume"
        fi
        flock $(which guppy_basecaller).{config[cuda]} guppy_basecaller --gpu_runners_per_device 12 --chunks_per_runner 2048 --num_callers {threads} $RESUME --config {config[guppy_config]} -x 'cuda:{config[cuda]}:100%'  \
        --records_per_fastq 100000 --save_path guppy_output/  --compress_fastq --fast5_out -i input_fast5_copy/
        """




rule extract_methylation_header:
    input:
        fast5s_done="guppy_output/guppy.done"
    output:
        header="{basename}.usam_header"
    log:
        "{basename}.header.log"
    threads:
        2
    shell:
        "python $(dirname {workflow.snakefile})/scripts/extract_methylation_fast5_to_sam.py --fastq {output.header} -V -L -F --failed_reads /dev/null -- guppy_output/workspace/*fast5  2>{log} |head -1 >/dev/null ||true;"

rule extract_methylation_likelihood_filter_fastq:
    input:
        fast5s_done="guppy_output/guppy.done"
        #fast5s=rules.call_methylation.output.basecalled
    output:
        fastq_pass=pipe("{basename}.pass_fastq"),
        fastq_fail=pipe("{basename}.fail_fastq"),
    threads:
        2
    shell:
        "python $(dirname {workflow.snakefile})/scripts/extract_methylation_fast5_to_sam.py --fastq -o {output.fastq_pass} --failed_reads {output.fastq_fail} -V -L -F -- guppy_output/workspace/*fast5"



rule map_fastq:
    input:
        unmapped="{basename}.{type}_fastq",
        header=rules.extract_methylation_header.output.header,
        reference=config["reference_fasta"]
    output:
        cram=protected("{basename}.{type,pass|fail}.minimap2.sort.cram")
    threads:
        3
    shell:
        "minimap2 -x map-ont -y -a -R \"$(grep -E ^@RG {input.header}|sed -e 's/\t/\\\\t/g' )\"  -t {threads}  {input.reference} {input.unmapped} |"
        "samtools sort -O cram -@ 5 -m 15G --reference {input.reference}  -o {output.cram} /dev/stdin;"
        "samtools index {output.cram}"

