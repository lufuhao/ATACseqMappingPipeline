#!/bin/bash
#set -o errexit
### Set readonly variable
#readonly passwd_file=�/etc/passwd�
### exit when variable undefined
#set -o nounset
### Script Root
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
### MachType
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi
#export NUM_THREADS=`grep -c \'^processor\' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"
RunPath=$PWD
echo "RunDir: $RunPath"



###echo color
#Black        0;30     Dark Gray     1;30
#Red          0;31     Light Red     1;31
#Green        0;32     Light Green   1;32
#Brown/Orange 0;33     Yellow        1;33
#Blue         0;34     Light Blue    1;34
#Purple       0;35     Light Purple  1;35
#Cyan         0;36     Light Cyan    1;36
#Light Gray   0;37     White         1;37
#RED=\'\033[0;31m\'
#NC=\'\033[0m\' # No Color
#printf "I ${RED}love${NC} Stack Overflow\n"



################# help message ######################################

help() {

cat<<HELP
	
$0 --- Brief Introduction
	
Version: 20210809



Requirements:
    Linux: perl, echo
    BEDtools
    bamaddrg
    samtools
    Rscript
    BEDOPS
    MACS2/3



Descriptions:
  This script is used to call ATAC-seq peaks
    1. merge 3 rmdup replicated BAMs
    2. sort by names for merged BAM and rmdup BAMs
    3. convert BAM to BED by 'bedtools bamtobed'
    4. BAMPE to BEDPE, and adjust 9bp for ATAC-seq
            Forward strand +4 and reverse -5
    5. MACS2/3 to call peaks
    6. BEDtools intersact BEDs to get consensus



Options:
  -h    -------    Print this help message
  -i    <B1,B2>    Replicated BAM file list
  -c    <ctrlB>    BAM control for naked DNA [Optional]
  -g    <SciNum>   Genome Size for MACS
  -p    <p1,p2>    Output prefix for Replicates
  -o    <Pfx>      Output prefix
  -D    <Path>     Running Path
  -m    <Str>      Option --max-mem for sort-bed, default 2G
  -t    <INT>      Number of threads for samtools merge; default: 1
  -S    -------    Suppress intersectBED step as sometime 
                     it fails due to chrom sort
  -macs <2/3>      Use MACS2/MACS3 to call ATAC peaks; default: 3
  


Example:
  $0 -i In1.bam,In2.bam -g 1.5e10  -p In1,In2 -o MyOut \
     -D /path/to/run -t 1



Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: Fu-Hao.Lu@jic.ac.uk

HELP

exit 0

}

[ $# -lt 1 ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH



#################### Initializing ###################################
RunDir=$PWD
declare -a BAMInArr=()
declare -a RepPfxArr=()
opt_D=$PWD
opt_t=1
opt_c=''
opt_S=0;
opt_m="2G";
opt_macs_vers=3;
tmpRunDir=''

#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) BAMInArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -c) opt_c=$2;shift 2;;
    -g) opt_g=$2;shift 2;;
    -p) RepPfxArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -o) opt_o=$2;shift 2;;
    -D) opt_D=$2;shift 2;;
    -m) opt_m=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -S) opt_S=1;shift 1;;
    -macs) opt_macs_vers=$2;shift 2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done



#################### Subfuctions ####################################

###Detect command existence

CmdExists () {
	if command -v $1 >/dev/null 2>&1; then
		return 0
	else
		return 1
	fi
}
### Name-sort BAM and convert BAMPE to BEDPE, apply 9 bp shift for ATAC
### Bampe2Bedpe (in.bam, out.namesort.bam, out.bed)
### Global: numthreads
### Program: macs_bedpe_from_bampe.sh
### Output: 
Bampe2Bedpe () {
	local BBbamin=$1;
	local BBbamout=$2;
	local BBbedout=$3;
	
	local BBsubinfo="FuhaoBash_BashMod(Bampe2Bedpe)"
	
	if [ ! -s "$BBbamout" ]; then
		if [ $samtoolsVers -eq 0 ]; then
			samtools sort -@ $opt_t -n $BBbamin ${BBbamout%.bam}
			if [ $? -ne 0 ] || [ ! -s $BBbamout ]; then
				echo "${BBsubinfo}Error: name sort error" >&2
				exit 100
			fi
		elif [ $samtoolsVers -eq 1 ]; then
			samtools sort -@ $opt_t -n -O BAM -o $BBbamout $BBbamin
			if [ $? -ne 0 ] || [ ! -s $BBbamout ]; then
				echo "${BBsubinfo}Error: name sort error" >&2
				exit 100
			fi
		fi
	fi
	$RootDir/scripts/macs_bedpe_from_bampe.sh -i $BBbamout -o $BBbedout
	if [ $? -ne 0 ] || [ ! -s $BBbedout ]; then
		echo "${BBsubinfo}Error: bam2bed error" >&2
		exit 100
	fi
	
	return 0;
}
randomDir () {
	local RDdir=$1
	
	local BStest=0
	local BSpfx="tmp_"
	
	while [ $BStest -eq 0 ]; do
#		local BSrand=$(cat /dev/urandom | LC_CTYPE=C tr -dc '0-9a-zA-Z' | fold -w 8 | head -n 1);
		local BSrand=$(perl -le 'BEGIN{@chars = ("0".."9", "a".."z", "A".."F");$len = 10;} while($len--){ $string .= $chars[rand @chars] }; print "$string";')
		tmpRunDir="$RDdir/${BSpfx}${BSrand}"
		if [ ! -d "$tmpRunDir" ]; then
			BStest=1;
			break
		fi
	done
	
	return 0
}
bedSort () {
	local BSin=$1
	local BSout=$2
	local BStmpdir=$3
	
	if [ -z "$BStmpdir" ]; then
		BStmpdir=$PWD
	fi
	
	randomDir "$BStmpdir"
	
	sort-bed --max-mem $opt_m --tmpdir $tmpRunDir $BSin > $BSout
	if [ $? -ne 0 ] || [ ! -s "$BSout" ]; then
		echo "Error: sort-bed error" >&2
		echo "CMDs: sort-bed --max-mem $opt_m --tmpdir $tmpRunDir $BSin > $BSout" >&2
		exit 100
	fi
	
	rm -rf $tmpRunDir > /dev/null 2>&1
	
	return 0
}
bedSort2 () {
	local BSin=$1
	local BSout=$2
	
	sort -k1,1 -k2,2n -k3,3n $BSin > $BSout
	if [ $? -ne 0 ] || [ $ -s $BSout ]; then
		echo "Error: BED sort 2 error" >&2
		echo "CMD: sort -k1,1 -k2,2n -k3,3n $BSin > $BSout" >&2
		exit 100
	fi
	
	return 0
}



#################### Command test ###################################

#CmdExists 'macs_bedpe_from_bampe.sh'
#if [ $? -ne 0 ]; then
if [ ! -s $RootDir/scripts/macs_bedpe_from_bampe.sh ]; then
	echo "Error: script '$RootDir/scripts/macs_bedpe_from_bampe.sh' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'bedtools'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'bedtools' in PROGRAM 'BEDtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'bamaddrg'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'bamaddrg' in PROGRAM 'bamaddrg' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'samtools'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
samtoolsVers=$(samtools 2>&1 | grep ^'Version' | sed 's/^Version:\s\+//;s/\..*$//;')
CmdExists 'sort-bed'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'sort-bed' in PROGRAM 'BEDOPS' is required but not found.  Aborting..." >&2 
	exit 127
fi

#################### Defaults #######################################
if [ ! -z "$opt_g" ]; then
	echo "Info: Genome Size: $opt_g"
else
	echo "Error: setup genome size for MACS" >&2
	exit 100
fi
cmd_macs="macs3"
if [ $opt_macs_vers -eq 2 ]; then
	cmd_macs="macs2"
	CmdExists 'macs2'
	if [ $? -ne 0 ]; then
		echo "Error: CMD 'macs2' in PROGRAM 'MACS2' is required but not found.  Aborting..." >&2 
		exit 127
	fi
elif [ $opt_macs_vers -eq 3 ]; then
	cmd_macs="macs3"
	CmdExists 'macs3'
	if [ $? -ne 0 ]; then
		echo "Error: CMD 'macs3' in PROGRAM 'MACS3' is required but not found.  Aborting..." >&2 
		exit 127
	fi
else
	echo "Error: please use -macs 2/3 to specify which MACS version to use:  2=MACS2 3=MACS3" >&2
	exit 100
fi



#################### Input and Output ###############################
if [ ! -d $opt_D ]; then
	mkdir -p $opt_D
fi



#################### Main ###########################################

if [ ! -d $opt_D ]; then
	mkdir -p $opt_D
fi
cd $opt_D



### 1. BEM to BEDPE
declare -a BedArr=()
BedToMerge=''
echo "### 1. BAM to BEDPE"; echo "### 1. BAM to BEDPE" >&2;
for ((BamNum=0;BamNum < ${#BAMInArr[@]};BamNum++)); do
	BedOut="${RepPfxArr[$BamNum]}.bed"
	if [ ! -s $BedOut ]; then
		echo "Info : generate BEDPE: ${BAMInArr[$BamNum]}"
		echo "Info : generate BEDPE: ${BAMInArr[$BamNum]}" >&2
		if Bampe2Bedpe ${BAMInArr[$BamNum]} ${RepPfxArr[$BamNum]}.namesort.bam $BedOut; then
			echo "Info: Bampe2Bedpe successful: ${BAMInArr[$BamNum]}"
		else
			echo "Info: Bampe2Bedpe failed: ${BAMInArr[$BamNum]}"
			exit 100
		fi
	else
		echo "Info : using existing BEDPE: ${BAMInArr[$BamNum]}"
		echo "Info : using existing BEDPE: ${BAMInArr[$BamNum]}" >&2
	fi
	BedSort="${RepPfxArr[$BamNum]}.sort.bed"
	if [ ! -s $BedSort ]; then
		echo "Info : generate sorted BEDPE: ${BAMInArr[$BamNum]}"
		echo "Info : generate sorted BEDPE: ${BAMInArr[$BamNum]}" >&2
		bedSort $BedOut $BedSort $PWD
#		bedSort2 $BedOut $BedSort"
	else
		echo "Info: using existing sorted BEDPE: $BedSort"
	fi
	BedToMerge="$BedToMerge $BedSort"
	BedArr+=("$BedSort")
	if [ ! -s "$BedSort.pdf" ]; then
		Rscript $RootDir/scripts/dist_bin_plot.insert.size.Rscript $BedSort $BedSort.pdf
	fi
done
###Control
BedControl="$opt_o.control.bed"
BedControlSort="$opt_o.control.sort.bed"
if [ -s "$opt_c" ]; then
	if [ ! -s "$BedControl" ]; then
		echo "Info : Generating Control BEDPE for BAM: $opt_c"
		echo "Info : Generating Control BEDPE for BAM: $opt_c" >&2
		if Bampe2Bedpe $opt_c $opt_o.control.namesort.bam $BedControl; then
			echo "Info: Bampe2Bedpe successful: $opt_c"
		else
			echo "Info: Bampe2Bedpe failed: $opt_c"
			exit 100
		fi
	else
		echo "Info : using existing BEDPE: $opt_c"
		echo "Info : using existing BEDPE: $opt_c" >&2
	fi
	
	if [ ! -s "$BedControlSort" ]; then
		echo "Info : generate CONTROL sorted BEDPE: $opt_c"
		echo "Info : generate CONTROL sorted BEDPE: $opt_c" >&2
		bedSort "$BedControl" "$BedControlSort" "$PWD"
#		bedSort2 $BedControl $BedControlSort
	else
		echo "Info: using existing control sorted BEDPE: $BedControlSort"
	fi
	if [ ! -s "$BedControlSort.pdf" ]; then
		Rscript $RootDir/scripts/dist_bin_plot.insert.size.Rscript $BedControlSort $BedControlSort.pdf
	fi
fi



### 2. merge BEDPE
echo "### 2. merging and sorting BEDs";echo "### 2. merging and sorting BEDs" >&2
OutBedpe1Merge="$opt_o.1.merge.bed"
OutBedpe2Sort="$opt_o.2.sort.bed"
if [ ! -s $OutBedpe1Merge ]; then
	echo "Info: merging BEDPE"
	echo "    CMD: cat $BedToMerge > $OutBedpe1Merge"
	cat $BedToMerge > $OutBedpe1Merge
	if [ $? -ne 0 ] || [ ! -s $OutBedpe1Merge ]; then
		echo "Error: BEDPE merge failed: $OutBedpe1Merge" >&2
		exit 100
	fi
else
	echo "Info: using existing merged BEDPE: $OutBedpe1Merge"
fi
if [ ! -s $OutBedpe2Sort ]; then
	echo "Info: sorting merged BEDPE"
	echo "Info: sorting merged BEDPE" >&2
	bedSort "$OutBedpe1Merge" "$OutBedpe2Sort" "$PWD"
#	bedSort2 $OutBedpe1Merge $OutBedpe2Sort"
else
	echo "Info: using existing sorted merged BEDPE: $OutBedpe2Sort"
fi
if [ ! -s "$OutBedpe2Sort.pdf" ]; then
	Rscript $RootDir/scripts/dist_bin_plot.insert.size.Rscript $OutBedpe2Sort $OutBedpe2Sort.pdf
fi


### 3. MACS
echo "### 3. Call peaks by $cmd_macs"; echo "### 3. Call peaks by $cmd_macs" >&2
if [ ! -s $OutBedpe2Sort ]; then
	echo "Error: invalid sorted merged BEDPE: $OutBedpe2Sort" >&2
	exit 100;
fi

MacsOptions="--verbose 3 --format BEDPE -g $opt_g --bdg --keep-dup all --nomodel --shift -37 --extsize 73"
#MacsOptions="--verbose 3 --format BEDPE -g $opt_g --bdg --keep-dup all --nomodel --shift -100 --extsize 200 -SPMR -n peak"
### -B, --bdg    save extended fragment pileup, and local lambda tracks (two files) at every bp into a bedGraph file.
### --SPMR       If True, MACS will SAVE signal per million reads for fragment pileup profiles. It won't interfere with computing pvalue/qvalue during peak calling, since internally MACS2 keeps using the raw pileup and scaling factors between larger and smaller dataset to calculate statistics measurements. If you plan to use the signal output in bedGraph to call peaks using bdgcmp and bdgpeakcall, you shouldn't use this option because you will end up with different results. However, this option is recommended for displaying normalized pileup tracks across many datasets. Require -B to be set. Default: False
### -n NAME, --name NAME  Experiment name, which will be used to generate output file names. DEFAULT: "NA"
if [ -s "$opt_c" ] && [ -s "$BedControlSort" ]; then
	MacsOptions="$MacsOptions --control $BedControlSort"
fi

MacsMergedPeaks="$opt_D/$opt_o/${opt_o}_peaks.narrowPeak"
if [ ! -d $opt_D/$opt_o ] || [ ! -s $MacsMergedPeaks ]; then
	randomDir "$opt_D"
	mkdir -p $tmpRunDir
	echo "Info: call peaks for merged BEDPE: $OutBedpe2Sort"
	echo "      CMD: $cmd_macs callpeak $MacsOptions --treatment $OutBedpe2Sort --outdir $opt_D/$opt_o --name $opt_o --tempdir $tmpRunDir > $opt_o.$cmd_macs.log 2>&1"
	$cmd_macs callpeak $MacsOptions --treatment $OutBedpe2Sort --outdir $opt_D/$opt_o --name $opt_o --tempdir $tmpRunDir > $opt_o.$cmd_macs.log 2>&1
	if [ $? -ne 0 ] && [ ! -s $MacsMergedPeaks ]; then
		echo "Info: $cmd_macs running error" >&2
		echo "      CMD: $cmd_macs callpeak $MacsOptions --treatment $OutBedpe2Sort --outdir $opt_D/$opt_o --name $opt_o --tempdir $tmpRunDir > $opt_o.$cmd_macs.log 2>&1" >&2
		exit 100
	fi
	rm -rf $tmpRunDir > /dev/null 2>&1
else
	echo "Info: using existing PEAKs for merged BEDPE: $OutBedpe2Sort"
fi
MacsMergedPeakSort="$opt_D/$opt_o/${opt_o}_peaks.narrowPeak.sorted.bed"
if [ ! -s $MacsMergedPeakSort ]; then
	echo "Info: sorting peaks for merged BEDPE: $MacsMergedPeaks"
	echo "Info: sorting peaks for merged BEDPE: $MacsMergedPeaks" >&2
	bedSort "$MacsMergedPeaks" "$MacsMergedPeakSort" "$PWD"
#	bedSort2 $MacsMergedPeaks $MacsMergedPeakSort
else
	echo "Info: using existing sorted merged Peaks: $MacsMergedPeakSort"
fi

if [ ${#BedArr[@]} -lt 1 ]; then
	echo "Error: invalid replicated BEDPEs" >&2
	exit 100
fi
FinalBed=$MacsMergedPeakSort
for ((BedpeNum=0; BedpeNum<${#BedArr[@]};BedpeNum++)); do
	MacsPeaks="$opt_D/${RepPfxArr[$BedpeNum]}/${RepPfxArr[$BedpeNum]}_peaks.narrowPeak"
	MacsPeakSort="$opt_D/${RepPfxArr[$BedpeNum]}/${RepPfxArr[$BedpeNum]}_peaks.narrowPeak.sorted.bed"
	if [ ! -d $opt_D/${RepPfxArr[$BedpeNum]} ]; then
		mkdir -p $opt_D/${RepPfxArr[$BedpeNum]}
	fi
	if [ ! -s $MacsPeaks ]; then
		randomDir "$opt_D"
		mkdir -p $tmpRunDir
		echo "Info: call peaks for BEDPE: ${BedArr[$BedpeNum]}"
		echo "      CMD: $cmd_macs callpeak $MacsOptions --treatment ${BedArr[$BedpeNum]} --outdir $opt_D/${RepPfxArr[$BedpeNum]} --name ${RepPfxArr[$BedpeNum]} --tempdir $tmpRunDir > $opt_o.$cmd_macs.log 2>&1"
		$cmd_macs callpeak $MacsOptions --treatment ${BedArr[$BedpeNum]} --outdir $opt_D/${RepPfxArr[$BedpeNum]} --name ${RepPfxArr[$BedpeNum]} --tempdir $tmpRunDir > $opt_o.$cmd_macs.log 2>&1
		if [ $? -ne 0 ] && [ ! -s $MacsPeaks ]; then
			echo "Error: $cmd_macs running error for BEDPE: ${BedArr[$BedpeNum]}" >&2
			echo "    CMD:$cmd_macs callpeak $MacsOptions --treatment ${BedArr[$BedpeNum]} --outdir $opt_D/${RepPfxArr[$BedpeNum]} --name ${RepPfxArr[$BedpeNum]} --tempdir $tmpRunDir > $opt_o.$cmd_macs.log 2>&1" >&2
			exit 100
		fi
		rm -rf $tmpRunDir > /dev/null 2>&1
	else
		echo "Info: using existing PEAKs for BEDPE: ${BedArr[$BedpeNum]}"
	fi
	if [ ! -s $MacsPeakSort ]; then
		echo "Info: sorting peaks for BEDPE: $MacsPeaks"
		bedSort "$MacsPeaks" "$MacsPeakSort" "$PWD"
#		bedSort2 $MacsPeaks $MacsPeakSort
	else
		echo "Info: using existing sorted merged Peaks: $MacsMergedPeakSort"
	fi
	if [ ! -s $MacsPeakSort.FRiP ]; then
		###
		echo "Analyzing FRiP: this function not finish yet"
	fi
	if [ $opt_S -eq 0 ]; then
		IntersectBed="$opt_o.bedtools.intersect.$BedpeNum.bed"
		if [ ! -s $IntersectBed ]; then
			echo "Info: intersect peaks for BEDPE: $MacsPeakSort"
			echo "      CMD: bedtools intersect -a $FinalBed -b $MacsPeakSort -wa -u -sorted > $IntersectBed"
			bedtools intersect -a $FinalBed -b $MacsPeakSort -wa -u -sorted > $IntersectBed
			if [ $? -ne 0 ] || [ ! -s $IntersectBed ]; then
				echo "Error: bedtools intersect running failed" >&2
				echo "    CMD: bedtools intersect -a $FinalBed -b $MacsPeakSort -wa -u -sorted > $IntersectBed" >&2
				exit 100
			fi
		else
			echo "Info: using existing intersect peaks: $IntersectBed"
		fi
		FinalBed=$IntersectBed
	else
		echo "Wanings: intersectBED ignored due to '-S' setting; Need to manually intersect the sorted peaks"
	fi
done


echo "#################### END #######################"
echo "#################### END #######################" >&2

exit 0
