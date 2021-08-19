#!/bin/bash
### calculate FRiP scores

# 1. prepare
## convert BAM (BAM used to call peaks) to BED
bedtools bamtobed -i ${sample}.sorted.marked.filtered.shifted.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > ${sample}.sorted.marked.filtered.shifted.PE2SE.tagAlign

# 2. total reads in BAM/BED
samtools view -c ${sample}.sorted.marked.filtered.shifted.bam
wc -l ${sample}.sorted.marked.filtered.shifted.PE2SE.tagAlign

# 3. count reads in peak regions
## 3.1 tagAlign, intersectBed -a tagAlign -b bed
time bedtools sort -i ${sample}_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a ${sample}.sorted.marked.filtered.shifted.PE2SE.tagAlign -b stdin | wc -l

#real    0m27.012s
#user    0m25.726s
#sys 0m1.357s

## 3.2 tagAlign, intersectBed -a bed -b tagAlign
time bedtools sort -i ${sample}_peaks.narrowPeak |bedtools merge -i stdin | bedtools intersect -c -a stdin -b ${sample}.sorted.marked.filtered.shifted.PE2SE.tagAlign | awk '{{ sum+=$4 }} END {{ print sum }}'

#real    0m51.089s
#user    0m39.199s
#sys 0m11.945s

## 3.3 BAM, intersectBed -a bam -b bed
time bedtools sort -i ${sample}_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a ${sample}.sorted.marked.filtered.shifted.bam -b stdin -ubam | samtools view -c

#real    1m12.844s
#user    1m11.979s
#sys 0m0.951s

## 3.4 BAM, intersectBed -a bed -b bam
time bedtools sort -i ${sample}_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -c -a stdin -b ${sample}.sorted.marked.filtered.shifted.bam | awk '{{ sum+=$4 }} END {{ print sum }}'

#real    1m49.981s
#user    1m28.747s
#sys 0m20.837s

## 3.5 featureCoutns

### covert BED (the peaks) to SAF
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ${sample}_peaks.narrowPeak > ${sample}_peaks.saf
### count
featureCounts -p -a ${sample}_peaks.saf -F SAF -o readCountInPeaks.txt ${sample}.sorted.marked.filtered.shifted.bam
