#!/bin/bash

#/usr/users/celldev/luf/test/atac/6.chipseeker/P25.annot.txt
#/usr/users/celldev/luf/test/atac/6.chipseeker/P5.annot.txt
#/usr/users/celldev/luf/test/atac/6.chipseeker/T25.annot.txt
#/usr/users/celldev/luf/test/atac/6.chipseeker/T5.annot.txt




#cat /usr/users/celldev/luf/test/atac/3.bampe/P25-1-A2_TKD180702135_HNFCCCCXY_L3.bowtie.clean.sort.exc.reheader.rmdup.uniqmapping.namesort.bam.bed /usr/users/celldev/luf/test/atac/3.bampe/P25-2-A3_TKD180702136_HNFCCCCXY_L3.bowtie.clean.sort.exc.reheader.rmdup.uniqmapping.namesort.bam.bed /usr/users/celldev/luf/test/atac/3.bampe/P25-3-B4_TKD180702137_merge.bowtie.clean.sort.exc.reheader.rmdup.uniqmapping.namesort.bam.bed | grep ^'chr' | sort -k1,1V -k2,2n -k3,3n > /usr/users/celldev/luf/test/atac/7.fragments/P25.chr3DL.fragments.bed

#cat /usr/users/celldev/luf/test/atac/3.bampe/P5-1-A5_TKD180702138_HNFCCCCXY_L3.bowtie.clean.sort.exc.reheader.rmdup.uniqmapping.namesort.bam.bed /usr/users/celldev/luf/test/atac/3.bampe/P5-2-B9_TKD180702139_HNFCCCCXY_L3.bowtie.clean.sort.exc.reheader.rmdup.uniqmapping.namesort.bam.bed /usr/users/celldev/luf/test/atac/3.bampe/P5-3-B10_TKD180702140_HNFCCCCXY_L3.bowtie.clean.sort.exc.reheader.rmdup.uniqmapping.namesort.bam.bed | grep ^'chr' | sort -k1,1V -k2,2n -k3,3n > /usr/users/celldev/luf/test/atac/7.fragments/P5.chr3DL.fragments.bed

#(grep ^'Chr3' /usr/users/celldev/luf/test/atac/3.bampe/T25-1-C1_TKD180702141_HNFCCCCXY_L3.bowtie.clean.sort.exc.rmdup.nmst.bam.bed; grep ^'Chr3' /usr/users/celldev/luf/test/atac/3.bampe/T25-2-C2_TKD180702142_HNFCCCCXY_L3.bowtie.clean.sort.exc.rmdup.nmst.bam.bed; grep ^'Chr3' /usr/users/celldev/luf/test/atac/3.bampe/T25-3-C3_TKD180702143_HNFCCCCXY_L3.bowtie.clean.sort.exc.rmdup.nmst.bam.bed) | perl -lane 'print if ($F[1]>249000000);' | sort -k1,1V -k2,2n -k3,3n > /usr/users/celldev/luf/test/atac/7.fragments/T25.chr3L.fragments.bed

#cat /usr/users/celldev/luf/test/atac/3.bampe/T5-1-C4_TKD180702144_HNFCCCCXY_L3.bowtie.clean.sort.exc.rmdup.nmst.bam.bed /usr/users/celldev/luf/test/atac/3.bampe/T5-2-C5_TKD180702145_HNFCCCCXY_L3.bowtie.clean.sort.exc.rmdup.nmst.bam.bed /usr/users/celldev/luf/test/atac/3.bampe/T5-3-C6_TKD180702146_merge.bowtie.clean.sort.exc.rmdup.nmst.bam.bed  | grep ^'Chr3' | sort -k1,1V -k2,2n -k3,3n > /usr/users/celldev/luf/test/atac/7.fragments/T5.chr3L.fragments.bed



#/usr/users/celldev/luf/test/atac/6.chipseeker/P25.annot.txt
#/usr/users/celldev/luf/test/atac/6.chipseeker/P5.annot.txt
#/usr/users/celldev/luf/test/atac/6.chipseeker/T25.annot.txt
#/usr/users/celldev/luf/test/atac/6.chipseeker/T5.annot.txt

#/usr/users/celldev/luf/test/atac/7.fragments/P25.chr3DL.fragments.bed
#/usr/users/celldev/luf/test/atac/7.fragments/P5.chr3DL.fragments.bed
#/usr/users/celldev/luf/test/atac/7.fragments/T25.chr3L.fragments.bed
#/usr/users/celldev/luf/test/atac/7.fragments/T5.chr3L.fragments.bed


SizeBin=5
RunDir=/usr/users/celldev/luf/test/atac/7.fragments
#declare -a PfxArr=("P25" "P5" "T25" "T5")
declare -a PfxArr=("P5" "T5")
for OutputPfx in ${PfxArr[@]}; do
	echo $OutputPfx
###file
	AnnotFile=/usr/users/celldev/luf/test/atac/6.chipseeker/$OutputPfx.annot.txt
	FragmentFile=$(ls $RunDir/$OutputPfx.*.fragments.bed)
	FinalSum="$RunDir/$OutputPfx.fragment.sum"
	FinalRange="$RunDir/$OutputPfx.peakrange.sum"
	
	echo "    Annotation file: $AnnotFile"
	echo "    Fragment file: $FragmentFile"
	echo "    Final summary: $FinalSum"
### Extract
	Temp1Annot="$RunDir/$OutputPfx.temp.UTR3downstream"
	Temp2Annot="$RunDir/$OutputPfx.temp.UTR5promoter"
	Temp3Annot="$RunDir/$OutputPfx.temp.Intergenic"
	Temp4Annot="$RunDir/$OutputPfx.temp.CDS"
	perl -lne '@F=split(/\t/); if ($F[12]=~/3.*UTR$/){print;}elsif ($F[12]=~/^Downstream\s+\(<1kb\)$/){print;}elsif ($F[12]=~/^Downstream\s+\(1-2kb\)$/){print;}' $AnnotFile | sort -k1,1V -k2,2n -k3,3n > $Temp1Annot
	if [ $? -ne 0 ] || [ ! -s $Temp1Annot ]; then
		echo "Error1: $Temp1Annot" >&2
		echo "CMD used: perl -lane 'if ($F[12]=~/3.*UTR$/){print;}elsif ($F[12]=~/^Downstream\s+\(<1kb\)$/){print;}elsif ($F[12]=~/^Downstream\s+\(1-2kb\)$/){print;}' $AnnotFile | sort -k1,1V -k2,2n -k3,3n > $Temp1Annot" >&2
		exit 100
	fi
	perl -lne '@F=split(/\t/); if ($F[12]=~/^5.*UTR$/){print;}elsif ($F[12]=~/^Promoter\s+\(<=1kb\)$/){print;}elsif ($F[12]=~/^Promoter\s+\(1-2kb\)$/){print;}' $AnnotFile | sort -k1,1V -k2,2n -k3,3n  > $Temp2Annot
	if [ $? -ne 0 ] || [ ! -s $Temp2Annot ]; then
		echo "Error2: $Temp2Annot" >&2
		exit 100
	fi
	perl -lne '@F=split(/\t/); if ($F[12]=~/^Distal\s+Intergenic$/){print;}elsif ($F[12]=~/^Downstream\s+\(2-3kb\)$/){print;}' $AnnotFile | sort -k1,1V -k2,2n -k3,3n  > $Temp3Annot
	if [ $? -ne 0 ] || [ ! -s $Temp3Annot ]; then
		echo "Error3: $Temp3Annot" >&2
		exit 100
	fi
	
	perl -lne '@F=split(/\t/); if ($F[12]=~/^Exon/){print;}elsif ($F[12]=~/^Intron/){print;}' $AnnotFile  | sort -k1,1V -k2,2n -k3,3n > $Temp4Annot
	if [ $? -ne 0 ] || [ ! -s $Temp4Annot ]; then
		echo "Error4: $Temp4Annot" >&2
		exit 100
	fi
### intersect
	declare -a CategArr=("$Temp1Annot" "$Temp2Annot" "$Temp3Annot" "$Temp4Annot")
	declare -a BinArr=()
	BinStr=''
	declare -a PeakRangeArr=()
	PeakRange=''
	for indcateg in ${CategArr[@]}; do
		echo "cat: $indcateg"
		bedtools intersect -a $FragmentFile -b $indcateg -wa -u -sorted > $indcateg.fil
		if [ $? -ne 0 ] || [ ! -s "$indcateg.fil" ]; then
			echo "Error5: $indcateg.fil" >&2
			echo "CMD used: bedtools intersect -a $FragmentFile -b $indcateg -wa -u > $indcateg.fil" >&2
			exit 100
		fi
		perl -lane 'print $F[2]-$F[1]' $indcateg.fil > $indcateg.fil.frag
		if [ $? -ne 0 ] || [ ! -s "$indcateg.fil.frag" ]; then
			echo "Error6: $indcateg.fil.frag" >&2
			exit 100
		fi
		SizeCollectBin_luf.pl $indcateg.fil.frag $SizeBin > $indcateg.fil.frag.bin$SizeBin
		if [ $? -ne 0 ] || [ ! -s "$indcateg.fil.frag.bin$SizeBin" ]; then
			echo "Error7: $indcateg.fil.frag.bin$SizeBin" >&2
			exit 100
		fi
		BinArr+=("$indcateg.fil.frag.bin$SizeBin")
		BinStr="$BinStr $indcateg.fil.frag.bin$SizeBin"
		perl -lane 'print $F[2]-$F[1]' $indcateg > $indcateg.peakrange
		if [ $? -ne 0 ] || [ ! -s "$indcateg.peakrange" ]; then
			echo "Error8: $indcateg.peakrange" >&2
			exit 100
		fi
		SizeCollectBin_luf.pl $indcateg.peakrange $SizeBin > $indcateg.peakrange.bin$SizeBin
		if [ $? -ne 0 ] || [ ! -s "$indcateg.peakrange.bin$SizeBin" ]; then
			echo "Error9: $indcateg.peakrange.bin$SizeBin" >&2
			exit 100
		fi
		PeakRangeArr+=("$indcateg.peakrange.bin$SizeBin")
		PeakRange="$PeakRange $indcateg.peakrange.bin$SizeBin"
	done
	if [ ${#BinArr[@]} -eq 4 ]; then
		list_merger.pl $FinalSum 0 $BinStr
		if [ $? -ne 0 ] || [ ! -s $FinalSum ]; then
			echo "Error10: $FinalSum" >&2
			exit 100
		fi
	else
		echo "Error11: not 4" >&2
		exit 100
	fi
	if [ ${#PeakRangeArr[@]} -eq 4 ]; then
		list_merger.pl $FinalRange 0 $PeakRange
		if [ $? -ne 0 ] || [ ! -s $FinalRange ]; then
			echo "Error12: $FinalRange" >&2
			exit 100
		fi
	else
		echo "Error13: not 4" >&2
		exit 100
	fi
done
