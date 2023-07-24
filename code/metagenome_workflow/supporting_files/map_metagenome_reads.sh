#!/bin/bash

module load bwa/0.7.17
module load samtools/1.16.1

## read mapping
inp=`ls <path_to_metagenome_files> | grep _1.fastq.gz""`
for R1 in $inp;
do R2=`echo $R1 | cut -d "_" -f1`_2.fastq.gz;
UB=`echo $R1 | cut -d "_" -f1`_unsorted.bam;
SB=`echo $R1 | cut -d "_" -f1`_sorted.bam;
SB_unique=`echo $R1 | cut -d "_" -f1`_unique_sorted.bam;
U_DEP=`echo $R1 | cut -d "_" -f1`_unique_depth.tab;
U_DEP_P=`echo $R1 | cut -d "_" -f1`_unique_depth-pos.tab;
U_TAB=`echo $R1 | cut -d "_" -f1`_unique.tab;
T_DEP=`echo $R1 | cut -d "_" -f1`_depth.tab;
T_DEP_P=`echo $R1 | cut -d "_" -f1`_depth-pos.tab;
T_TAB=`echo $R1 | cut -d "_" -f1`_total.tab;

# initial mapping
bwa mem -t 36 RAG_DNA_rm.fasta <path_to_metagenome_files>/$R1 <path_to_metagenome_files>/$R2 | samtools view -@ 18 -uS - -o $UB;

# sort bam
samtools sort -@ 18 $UB -o $SB;

# extract unique
samtools view -@ 18 -q 1 -f 2 -u $SB -o $SB_unique;
samtools index -@ 18 $SB_unique;
samtools idxstats $SB_unique > $U_DEP;
samtools depth $SB_unique > $U_DEP_P;
python /data/blausteinra/CLEAN/software/MGS-gut/scripts/parse_bwa-depth.py $U_DEP $U_DEP_P > $U_TAB;

# extract total
samtools index -@ 18 $SB;
samtools idxstats $SB > $T_DEP;
samtools depth $SB > $T_DEP_P;
python /data/blausteinra/CLEAN/software/MGS-gut/scripts/parse_bwa-depth.py $T_DEP $T_DEP_P > $T_TAB;

# remove files
rm -rf $UB $SB_unique $U_DEP $U_DEP_P $T_DEP $T_DEP_P;
done

