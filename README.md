# GpC methylation from fast5:s to reference anchored calls (gcf52ref)

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




