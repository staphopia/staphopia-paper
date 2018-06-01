#! /bin/bash
REFERENCE=$1
PREFIX=$2
SEED=123456

mkdir -p fastq

echo "Simulating Error Free (100bp) 300x coverage."
art_illumina -i ${REFERENCE} -na -p -rs ${SEED} -l 100 -f 300 -m 200 \
    -s 10 -o fastq/${PREFIX}-EF_R -ir 0 -ir2 0 -dr 0 -dr2 0 -qL 35

echo "Simulating HiSeq 2000 (100bp) 300x coverage."
art_illumina -ss HS20 -i ${REFERENCE} -na -p -rs ${SEED} -l 100 -f 300 -m 200 \
    -s 10 -o fastq/${PREFIX}-HS_R

echo "Simulating MiSeq V3 (250bp) 300x coverage."
art_illumina -ss MSv3 -i ${REFERENCE} -na -p -rs ${SEED} -l 250 -f 300 -m 500 \
    -s 10 -o fastq/${PREFIX}-MS_R

echo "Compressing FASTQ files"
gzip --fast fastq/${PREFIX}*.fq
