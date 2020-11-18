# GpC methylation from fast5:s to reference anchored calls (gcf52ref)

Set of scripts to convert guppy produced methylation calls from fast5 files to reference anchored
files similar to nanopolish output.

NOTE: This is first quick test for converting guppy methylation call fast5:s to something more
useful. Unfortunately the process is not very efficient and the output is not very useful. I'm
looking forward for other toolchain built around [SAM/BAM/CRAM
format](https://github.com/samtools/hts-specs) and [MM and MP
tags](https://github.com/samtools/hts-specs/pull/418).

## Requirements

Mostly ONT guppy basecaller and many other more or less standard bioinformatics tools. Conda
environment defined in `environment.yml` will define them all.

# Usage

See `test/test_gcf52ref.sh` for demonstration of downloading reference genome and fast5 file for
4000 reads from http://s3.amazonaws.com/nanopore-human-wgs/rel6/MultiFast5Tars/, rebasecalling,
creating tags and mapping the reads to genome.  The heavy lifting is done in a snakemake pipeline
launched by `scripts/fast5_to_cram.sh`

## Methylation to sam

Script `scripts/extract_methylation_fast5_to_sam.py` contains a very early code for outputting the
methylation calls from guppy to SAM format using MM/MP tags in [pull
request](https://github.com/samtools/hts-specs/pull/418) for SAM specification by James Bonfield.
Specifically the tags follow
[commit](https://github.com/samtools/hts-specs/pull/418/commits/11d7fb900b6d51417f59d7cd2cc8540c1b982590)
Oct 21, 2019.

In addition to MM/MP tags, the script can output the modification likelihoods are provided in hex
string tags such that `ml:Z:C+m,ffff04,A+a,000093`. The syntax is
`ml:Z:([ACGTN][-+][a-z],([0-9a-f][0-9a-f])+);)+` where `ml` is the custom SAM tag, `Z` symbol for character
string value. Next three groupings are as in 'Base modifications' in
https://github.com/samtools/hts-specs/pull/418 Final group is the hexadecimal coded likelihood for
the given type modification for each position in SEQ in the *original* strand. Note that the
likelihoods are 'given the underlying base is called the one defined in the tag'. Each hexadecimal
value ranging 00-ff is 255*P(data|base, modified).



# Commandline

Launch script for processing reads `input/path/for_fast5s/pass_fast5/*.fast5` to aligned cram files.  This requires `guppy_basecaller` in `PATH`

```
usage:
scripts/fast5_to_cram.sh  -i input/path/for_fast5s/ -s SAMPLEID 
-C all     GPU to use (0,1,2,3 or all, default all)
-c guppy_config.cfg  Guppy configuration file. Default /home/kpalin/ont-guppy/data/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg
-w WORKDIR Directory where to make the processes files and output.
-r REF_FASTA   Reference genome to use.
-s SAMPLEID  
-S server.remote.com   Remote SFTP server holding the data.
-K         Keep temporary files (except output )
-n         Dry run
-N 8       Number of splits
-h          Show this message and exit.
```

Python code used by above script to convert fast5:s to fastq:s.

```
usage: extract_methylation_fast5_to_sam.py [-h] [-o OUTPUT] [-f [FASTQ]]
                                           [--failed_reads FAILED_READS] [-L]
                                           [-F [FILTER]] [-V]
                                           input_fast5 [input_fast5 ...]

Extract base modifications from fast5 files called with Guppy 3.3. The
modifications are provided as MM and MP tags conforming to
 If requested,
the modification likelihoods are provided in hex string tags such that
ml:Z:C+m,ffff04,A+a,000093. The syntax is
ml:Z:([ACGTN][-+][a-z],([0-9a-f][0-9a-f])+);)+ ml is the custom SAM tag, Z
symbol for character string value. Next three groupings are as in 'Base
modifications' in https://github.com/samtools/hts-specs/pull/418 Final group
is the hexadecimal coded likelihood for the given type modification for each
position in SEQ in the original strand. Note that the likelihoods are 'given
the underlying base is called the one defined in the tag'. Each hexadecimal
value ranging 00-ff is 255*P(data|base, modified).

positional arguments:
  input_fast5           Input paths of fast5 files [default:None]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output the unsorted SAM or passing reads fastq file
                        file here [default:/dev/stdout]
  -f [FASTQ], --fastq [FASTQ]
                        Produce output in fastq format instead of SAM. Store
                        SAM header to file named here. The SAM tags are stored
                        as read comments that can be copied over to SAM
                        [default:None const:/dev/null]
  --failed_reads FAILED_READS
                        Output failed reads in fastq format here. With SAM
                        output, the filter/fail is marked with flag
                        [default:/dev/stderr]
  -L, --likelihoods     Include also the raw likelihoods as 'ml' tag.
                        [default:False]
  -F [FILTER], --filter [FILTER]
                        Mark reads with average q less than this as vendor
                        failed [default:False]
  -V, --verbose         Be more (and more) verbose with output [default:0]
```



There are various notes and to-do:s involved:

- As *output* from the script, the MM/MP tags are reported correctly in the '+' strand with respect
  to the SEQ record in the SAM.  The suggested minimap2 (or essentially any) alignment will reverse
  complement some of the reads but the MM/MP tags will still be in the original strand coordinates.
- The so-called likelihoods in guppy fast5 files are not likelihoods but (possibly/likely) some sort
  of posterior probabilities of modification given a call and nanopore signal:  $P(mod|base call,
  signal) = P(mod, base call, signal)/ ( P(base call|signal) * P(signal) )$
- Current (commit `1e10d0dc`) version calls methylated C:s in any context but logged statistics are
  computed on CpG sites.
- The calls are made assuming that ONT gives likelihoods so we get $P(mod|D) = P(D|mod)*P(mod)/P(D)
  = P(D|mod)/( P(D|mod)P(mod)+P(D|no-mod)P(no-mod))$ and assuming $P(mod)=0.9$ and $P(no-mod) =
  0.1$.  This prior results in 66.5% of CpG sites called methylated with phred score>=3 on a single
  colorectal adenocarcinoma sample.
