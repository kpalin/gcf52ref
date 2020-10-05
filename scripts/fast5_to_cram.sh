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


GUPPY_CONFIG=$(readlink -f $HOME/ont-guppy/data/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg )
CUDA="all"
usage()  {
    echo -e "usage:"
echo "$0  -i input/path/for_fast5s/ -s SAMPLEID 
-C all     GPU to use (0,1,2,3 or all, default $CUDA)
-c guppy_config.cfg  Guppy configuration file. Default ${GUPPY_CONFIG}
-s SAMPLEID  
-K         Keep temporary files (except output )
-h          Show this message and exit." >&2
    exit 1

}
CMDS=""
while getopts  "i:s:C:hc:w:nK" flag
do
    case "${flag}" in 
        s)
              SAMPLEID="${OPTARG}"
        ;;
        C)
              CUDA="${OPTARG}"
        ;;        
        c)
              GUPPY_CONFIG="${OPTARG}"
        ;;        
        w)
            mkdir -p "${OPTARG}"
            cd "${OPTARG}"
        ;;
        i)
              FAST5PATH="${OPTARG}"
        ;;
        K)
            CMDS=${CMDS}" --notemp"
        ;;
        n)
            CMDS=${CMDS}" -n"
        ;;
        h|*)
              usage
        ;;
    esac
done
shift $((OPTIND-1)); OPTIND=1

ls ${FAST5PATH}/sequencing_summary/*sequencing_summary.txt* || ls ${FAST5PATH}/*sequencing_summary.txt* || echo "WARNING: NO OLD sequencing_summary.txt FILES FOUND" >&2

test -e ${GUPPY_CONFIG}


test -e config.json || (jq -nR "{timestamp:\"$(date -Is)\",
reference_fasta:\"${REFERENCE_FASTA}\",
sampleid:\"${SAMPLEID}\",
guppy_config:\"${GUPPY_CONFIG}\",
fast5_path:\"${FAST5PATH}\",
cuda:\"${CUDA}\"}" >config.json )


cat config.json
SLEEPTIME=$(( $RANDOM % 60))s
echo Sleeping $SLEEPTIME
#sleep $SLEEPTIME
snakemake --latency-wait 60 --resources disk_io=100 ${CMDS:-} --rerun-incomplete --config nsplits=8 --cores=60 --snakefile ${SCRIPTBASE}/fast5_to_mapped_cram_map_individual_split.smk -p --configfile config.json
