#!/bin/bash

#  Created on Wednesday, 18 November  2020

# For debugging
#set -o verbose 

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail

CONDA_DIR=$(dirname "$(which conda)" ) 
MAMBA=$(which mamba || echo conda)
TESTDIR=$(dirname "$(readlink -f "$0")")/


echo "Setting up conda" >&2
# Conda init for this script.
source "${CONDA_DIR}/../etc/profile.d/conda.sh"
cd "${TESTDIR}"

test -d "${TESTDIR}/gcf52ref_test_env/conda-meta/"  || \
    "${MAMBA}" env create -p "${TESTDIR}/gcf52ref_test_env" --file "${TESTDIR}/../environment.yml"

set +u
conda activate  "${TESTDIR}/gcf52ref_test_env"
set -u 



echo "Fetching test FAST5 file" >&2
FAST5TAR="http://s3.amazonaws.com/nanopore-human-wgs/rel6/MultiFast5Tars/FAF15630-4244782843_Multi_Fast5.tar"
FAST5NAME="Notts/FAF15630-4244782843_Multi/PLSP61583_20170309_FNFAF15630_MN17073_sequencing_run_Nott_RAD002_wh2_0.fast5"
test -e "${FAST5NAME}" || curl  ${FAST5TAR} |tar -b 1024 -xv ${FAST5NAME}  



echo "Moving fast5 to 'classified' folder ./fast5_pass/.  " >&2
mkdir -p "$(dirname ${FAST5NAME})"/fast5_pass
test -e "$(dirname ${FAST5NAME})"/fast5_pass/"$(basename ${FAST5NAME})" || ln "${FAST5NAME}" \
    "$(dirname ${FAST5NAME})"/fast5_pass/ 



echo "Fetching reference genome and compressing with bgzip" >&2
REF_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
REF_NAME="$(basename ${REF_URL})"
test -e ${REF_NAME} ||  ( curl -o tmp_ref.fasta.gz "${REF_URL}" ;
    gunzip tmp_ref.fasta.gz;
    mv tmp_ref.fasta $(basename ${REF_NAME} .gz);
    bgzip $(basename ${REF_NAME} .gz)
)



echo "Running fast5_to_cram.sh" >&2
bash  "${TESTDIR}"/../scripts/fast5_to_cram.sh -s NA1278 -r "$(readlink -f "${REF_NAME}")" \
    -w "${TESTDIR}/NA1278_gcf52ref" \
    -i "$(dirname "$(readlink -f "${FAST5NAME}")" )" \
    -N 1 # Need to define single split (default 8) due to only single input fast5

