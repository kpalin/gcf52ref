#REFERENCE="/data/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta"
#import glob 
#FAST5S=list(glob.glob("/data/kpalin/Fam_c745_1_9623TK/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac/workspace/PAD64960_b924920b26f29f6118f165b040093705d5beb46b_1_fast5_*.fast5"))

import os.path
wildcard_constraints:
    type="pass|fail",
    f5type="pass|fail",
    fast5name="[a-zA-Z0-9_]+",
    suffix="|[.]gz"

report: "run_report/methylmapping.rst"

rule all:
    input:
        config["sampleid"]+".pass.minimap2.sort.cram",
        config["sampleid"]+".fail.minimap2.sort.cram"






in_wildcards = glob_wildcards(os.path.join(config["fast5_path"],"fast5_{f5type}/{fast5name}.fast5{suffix,|.gz}"))
#print("\n".join(in_wildcards.fast5name))
rule copy_fast5:
    input:
        fast5_file=os.path.join(config["fast5_path"],"fast5_{f5type}","{fast5name}.fast5")
    output:
        local_file=temp("input_fast5_copy/{f5type}_{fast5name}.fast5")
    shell:
        "cp {input.fast5_file} {output.local_file}"

rule uncompress_fast5:
    input:
        fast5_file=os.path.join(config["fast5_path"],"fast5_{f5type}/{fast5name}.fast5.gz")
    output:
        local_file=temp("input_fast5_copy/{f5type}_{fast5name}.gz.fast5")
    shell:
        "zcat {input.fast5_file} > {output.local_file}"



rule call_methylation:
    input:
        input_fast5s="input_fast5_copy/{f5type}_{fast5name}{suffix}.fast5"
    output:
        basecalled=directory("guppy_output/{f5type}_{fast5name}{suffix}/workspace/")
    threads: 13
    benchmark: "guppy_output/{f5type}_{fast5name}{suffix}/workspace/guppy.time.txt"
    shell:
        """
        CUDAID_=$(( (${{RANDOM}} + $$) % 4 ))
        flock $(which guppy_basecaller).${{CUDAID_}} guppy_basecaller --gpu_runners_per_device 12 --chunks_per_runner 2048 \
            --num_callers {threads} $RESUME --config {config[guppy_config]} -x 'cuda:${{CUDAID_}}:100%'  \
            --records_per_fastq 100000 --save_path guppy_output/{wildcards.f5type}_{wildcards.fast5name}{wildcards.suffix}/ \
                --compress_fastq --fast5_out -i input_fast5_copy/{input.input_fast5s}/
        """




rule extract_methylation_likelihood_filter_fastq:
    input:
        fast5s=rules.call_methylation.output.basecalled
    output:
        fastq_pass=pipe("{f5type}_{fast5name}{suffix}.pass_fastq"),
        fastq_fail=pipe("{f5type}_{fast5name}{suffix}.fail_fastq"),
        header=temp("mapped/{f5type}_{fast5name}{suffix}.sam.header")
    threads:
        2
    shell:
        "python $(dirname {workflow.snakefile})/scripts/extract_methylation_fast5_to_sam.py --fastq {output.header} "
        "-o {output.fastq_pass} --failed_reads {output.fastq_fail} -V -L -F -- {input.fast5s}/*.fast5"



rule map_fastq:
    input:
        unmapped="{f5type}_{fast5name}{suffix}.{type}_fastq",
        reference=config["reference_fasta"]
    output:
        cram=temp("mapped/{f5type}_{fast5name}{suffix}.{type,pass|fail}.minimap2.sort.bam")
    threads:
        3
    shell:
        "minimap2 -x map-ont -y -a -t {threads}  {input.reference} {input.unmapped} |"
        "samtools sort -O bam -l 0  -@ 5 -m 15G --reference {input.reference}  -o {output.cram} /dev/stdin"

rule merge_bams:
    input:
        mapped=expand("mapped/{f5type}_{fast5name}{suffix}.{{type}}.minimap2.sort.bam",zip,f5type=in_wildcards.f5type,
                                                                fast5name=in_wildcards.fast5name,suffix=in_wildcards.suffix),
        header="mapped/{f5type}_{fast5name}{suffix}.sam.header".format(f5type=in_wildcards.f5type[0],
                                                fast5name=in_wildcards.fast5name[0],suffix=in_wildcards.suffix[0]),
        reference=config["reference_fasta"]
    output:
        cram=protected(config["sampleid"]+".{type}.minimap2.sort.cram")
    threads:
        10
    shell:
        """RGID=$(sed -e '/^@RG/!d
s/^@RG.*ID:\([^\t]\+\).*$/\1/' {input.header});"""
        "samtools merge --threads {threads} -h {input.header} -u {input.mapped}|"
        "samtools add addreplacerg -R ${{RGID}} -o {output.cram} -O cram --reference {input.reference}"
