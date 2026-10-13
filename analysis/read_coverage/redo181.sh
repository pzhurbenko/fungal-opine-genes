#!/bin/bash
: "${GMO_ROOT:=$(cd "$(dirname "$0")/../.." && pwd)}"; export GMO_ROOT
cd ${GMO_ROOT}/analysis/read_coverage/reads
wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR285/081/SRR28540181/SRR28540181_1.fastq.gz && gzip -t SRR28540181_1.fastq.gz && gzip -t SRR28540181_2.fastq.gz && echo "REDL_OK $(date)" >> ../timing.txt && cd .. && bsub -q medium -n 32 -R 'span[hosts=1]' -R 'rusage[mem=1500]' -J map_GCA_038093215.1 -e lsf_map_GCA_038093215.1.err -o lsf_map_GCA_038093215.1.out "bash $PWD/map_one.sh GCA_038093215.1 SRR28540181"
