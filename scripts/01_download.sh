#!/usr/bin/env bash
set -euo pipefail
#
# Download raw sequencing reads for Zhuang et al. 2019.
#   GEO GSE139313 / SRA SRP226783 / BioProject PRJNA579227
#
# ~10 GB gzipped total. Requires sra-tools (conda install -c bioconda sra-tools).
#
# fastq-dump --gzip is used rather than fasterq-dump: fasterq-dump is faster
# but writes UNCOMPRESSED fastq (~6 GB/sample, ~45 GB total). Disk was the
# binding constraint. MAGeCK reads .gz directly, so nothing is lost but time.

mkdir -p data && cd data

for SRR in SRR10338906 SRR10338907 SRR10338908 SRR10338909 \
           SRR10338910 SRR10338911 SRR10338912; do
    if [[ -f "${SRR}.fastq.gz" ]]; then
        echo ">>> ${SRR} already present, skipping"
        continue
    fi
    echo ">>> ${SRR}  $(date)"
    fastq-dump --gzip "${SRR}"
done

echo ">>> done"
ls -lh *.fastq.gz
