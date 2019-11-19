#!/bin/bash
# kimmo.palin@helsinki.fi 

# For debugging
#set -o verbose

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail



SCRIPTBASE=$(readlink -f $(dirname $0)/..)
REFERENCE_FASTA="/data/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta"
#FAST5PATH="/mnt/cgnano/projects/promethion/190715/c985_1_9835/20190715_1211_1-A3-D3_PAD64026_1b4dc69c/"


GUPPY_CONFIG=$(readlink -f $HOME/ont-guppy/data/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg )
CUDA="all"
usage()  {
    echo -e "usage:"
echo "$0  -i input/path/for_fast5s/ -s SAMPLEID 
-C all     GPU to use (0,1,2,3 or all, default $CUDA)
-s SAMPLEID  
-h          Show this message and exit." >&2
    exit 1

}

while getopts  "i:s:C:hw:" flag
do
    case "${flag}" in 
        s)
              SAMPLEID="${OPTARG}"
        ;;
        C)
              CUDA="${OPTARG}"
        ;;        
        w)
            mkdir -p "${OPTARG}"
            cd "${OPTARG}"
        ;;
        i)
              FAST5PATH="${OPTARG}"
        ;;
        h|*)
              usage
        ;;
    esac
done
shift $((OPTIND-1)); OPTIND=1


#readlink -f ${@} |jq -nR '{methcalled_fast5:[inputs | select(length>0)]}' |\
#jq ".timestamp=\"$(date -Is)\"|.reference_fasta=\"${REFERENCE_FASTA}\"|.sampleid=\"${SAMPLEID}\"" #>config.json


#readlink -f ${@} |jq -nR  methcalled_fast5:[inputs | select(length>0)],'"


jq -nR "{timestamp:\"$(date -Is)\",
reference_fasta:\"${REFERENCE_FASTA}\",
sampleid:\"${SAMPLEID}\",
guppy_config:\"${GUPPY_CONFIG}\",
fast5_path:\"${FAST5PATH}\",
cuda:\"${CUDA}\"}" >config.json


cat config.json

snakemake --cores=50 --snakefile ${SCRIPTBASE}/fast5_to_mapped_cram.smk -p --configfile config.json 
