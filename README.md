[# GpC methylation from fast5:s to reference anchored calls (gcf52ref)

Set of scripts to convert guppy produced methylation calls from fast5 files to reference anchored files similar to nanopolish output.

NOTE: This is first quick test for converting guppy methylation call fast5:s to something more useful. Unfortunately the process is not very efficient and the output is not very useful. I'm looking forward for other 
toolchain built around [SAM/BAM/CRAM format](https://github.com/samtools/hts-specs) and [MM and MP tags](https://github.com/samtools/hts-specs/pull/418).

## Requirements

Mostly 'rocksdb' and many other more or less standard bioinformatics tools. Conda environment defined in `environment.yml` will define them all.

## Usage

The tools is very ad-hoc and not in any way recommended for production use and one should really think if this is sensible even for testing.  The process is two step: First convert the methylation calls from guppy 
called fast5 files to rocksdb key-value database. Second reads aligned cram/bam/sam file and outputs a tab separated file similar to [nanopolish](https://github.com/jts/nanopolish) methylcall command.


The first step `extract_methylation_fast5.py` reads the fast5 files and stores the CpG methylation likelihoods (value) for each read-id (key) in a rocksdb key-value store.  It is possible to run this step concurrently 
using multiple processes but that is always slower than running single thread (The hdf5/fast5 format is truly horendious)

Second step `extract_methylation_from_rocks.py` reads pre-aligned cram file of the reads and the rocksdb database from the first step and produces the reference achored output.



# Methylation to sam

Script `scripts/extract_methylation_fast5_to_sam.py` contains a very early code for outputting the methylation calls from guppy to SAM format using MM/MP tags in [pull request](https://github.com/samtools/hts-specs/pull/418) for SAM specification by James Bonfield. 

If you try out this script, do note: IT IS WHOLLY UNTESTED APART FROM OUTPUTTING SAM ACCEPTABLE BY samtools.

Usage for this would be something like: Extract methylation data to SAM/CRAM, map to reference and merge methylation and alignment information (As usual, picard will probably give you most problems)

```bash
python extract_methylation_fast5_to_sam.py -V -- PAD64960_*.fast5 | samtools sort -n -m 5G -@ 10 -O cram -o PAD64960.nsort.cram
samtools fastq -@ 5 PAD64960.nsort.cram|minimap2 -x map-ont -a -t 50 GRCh38_no_alt.fasta - |samtools sort -O cram -@ 5 -m 15G --reference GRCh38_no_alt/GRCh38_no_alt.fasta  -o PAD64960.minimap2.sort.cram
samtools index PAD64960.minimap2.sort.cram
picard -Xmx50G -Xms20G MergeBamAlignment R=GRCh38_no_alt.fasta ALIGNED=PAD64960.minimap2.sort.cram UNMAPPED=PAD64960.nsort.cram O=PAD64960s.minimap2.sort.meth.cram

```


There are various notes and to-do:s involved:

- As output from the script, the MM/MP tags are reported correctly in the '+' strand with respect to the SEQ record in the SAM.  The suggested minimap2 (or essentially any) alignment will reverse complement some of the reads but the MM/MP tags will still be in the original strand coordinates.
- Current (commit `1e10d0dc`) version is quite slow (150 reads/second on PromethION compute unit). I think much of this is because of logging and tracking issues that can be easily and safely removed when ready.
- The so-called likelihoods in guppy fast5 files are not likelihoods but (possibly/likely) some sort of posterior probabilities of modification given a call and nanopore signal:  $P(mod|base call, signal) = P(mod, base call, signal)/ ( P(base call|signal) * P(signal) )$
- Current (commit `1e10d0dc`) version calls methylated C:s in any context but logged statistics are computed on CpG sites.
- The calls are made assuming that ONT gives likelihoods so we get $P(mod|D) = P(D|mod)*P(mod)/P(D) = P(D|mod)/( P(D|mod)P(mod)+P(D|no-mod)P(no-mod))$ and assuming $P(mod)=0.9$ and $P(no-mod) = 0.1$.  This prior results in 66.5% of CpG sites called methylated with phred score>=3 on a single colorectal adenocarcinoma sample.
- It would probably be useful to store the ONT 'likelihoods' in SAM tag verbatim, as there appears to be positive 'likelihood' also for bases that are called something else (e.g. likely 5mC if the base would have been called C, but it was actually called A)]