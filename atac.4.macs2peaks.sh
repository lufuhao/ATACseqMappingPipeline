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
	
Version: 20190723



Requirements:
    Linux: perl, echo
    BEDtools
    bamaddrg
    macs2_bedpe_from_bampe.sh



Descriptions:
  This script is used to call ATAC-seq peaks
    1. merge 3 rmdup replicated BAMs
    2. sort by names for merged BAM and rmdup BAMs
    3. convert BAM to BED by 'bedtools bamtobed'
    4. BAMPE to BEDPE, and adjust 9bp for ATAC-seq
            Forward strand +4 and reverse -5
    5. MACS2 to call peaks
    6. BEDtools intersact BEDs to get consensus



Options:
  -h    -------    Print this help message
  -i    <B1,B2>    Replicated BAM file list
  -c    <ctrlB>    BAM control for naked DNA
  -g    <SciNum>   Genome Size for MACS2
  -p    <p1,p2>    Output prefix for Replicates
  -o    <Pfx>      Output prefix
  -D    <Path>     Running Path
  -t    <INT>      Number of threads for samtools merge
  -S    -------    Supress intersectBED step as sometime 
                     it fails due to chrom sort


Example:
  $0 -i In1.bam,In2.bam -g 1.5e10  -p In1,In2 -o MyOut \
     -D /path.to/run -t 1



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
    -t) opt_t=$2;shift 2;;
	-S) opt_S=1;shift 1;;
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
### Program: macs2_bedpe_from_bampe.sh
### Output: 
Bampe2Bedpe () {
	local BBbamin=$1;
	local BBbamout=$2;
	local BBbedout=$3;
	
	local BBsubinfo="FuhaoBash_BashMod(Bampe2Bedpe)"
	
	if [ ! -s "$BBbamout" ]; then
		samtools sort -@ $opt_t -n $BBbamin ${BBbamout%.bam}
		if [ $? -ne 0 ] || [ ! -s $BBbamout ]; then
			echo "${BBsubinfo}Error: name sort error" >&2
			exit 100
		fi
	fi
	macs2_bedpe_from_bampe.sh -i $BBbamout -o $BBbedout
	if [ $? -ne 0 ] || [ ! -s $BBbedout ]; then
		echo "${BBsubinfo}Error: bam2bed error" >&2
		exit 100
	fi
	
	return 0;
}


#################### Command test ###################################

CmdExists 'macs2_bedpe_from_bampe.sh'
if [ $? -ne 0 ]; then
	echo "Error: script 'macs2_bedpe_from_bampe.sh' is required but not found.  Aborting..." >&2 
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


#################### Defaults #######################################
if [ ! -z "$opt_g" ]; then
	echo "Info: Genome Size: $opt_g"
else
	echo "Error: setup genome size for MACS2" >&2
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
		echo "      CMD: sort -k1,1 -k2,2n -k3,3n $BedOut > $BedSort"
		sort -k1,1 -k2,2n -k3,3n $BedOut > $BedSort
		if [ ! -s $BedSort ]; then
			echo "Error: BEDPE sort failed: $BedSort" >&2
			exit 100
		fi
	else
		echo "Info: using existing sorted BEDPE: $BedSort"
	fi
	BedToMerge="$BedToMerge $BedSort"
	BedArr+=("$BedSort")
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
		echo "Info : generate sorted BEDPE: $opt_c" >&2
		echo "     CMD: sort -k1,1 -k2,2n -k3,3n $BedControl > $BedControlSort"
		sort -k1,1 -k2,2n -k3,3n $BedControl > $BedControlSort;
		if [ ! -s $BedControlSort ]; then
			echo "Error: BEDPE sort failed: $BedControlSort" >&2
			exit 100;
		fi
	else
		echo "Info: using existing control sorted BEDPE: $BedControlSort"
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
	echo "    CMD: sort -k1,1 -k2,2n -k3,3n $OutBedpe1Merge > $OutBedpe2Sort"
	sort -k1,1 -k2,2n -k3,3n $OutBedpe1Merge > $OutBedpe2Sort
	if [ $? -ne 0 ] || [ ! -s $OutBedpe2Sort ]; then
		echo "Error: BEDPE sort failed: $OutBedpe2Sort" >&2
		exit 100
	fi
else
	echo "Info: using existing sorted merged BEDPE: $OutBedpe2Sort"
fi



### 3. MACS2
echo "### 3. Call peaks by MACS2"; echo "### 3. Call peaks by MACS2" >&2
if [ ! -s $OutBedpe2Sort ]; then
	echo "Error: invalid sorted merged BEDPE: $OutBedpe2Sort" >&2
	exit 100;
fi

Macs2Options="--verbose 3 --format BEDPE -g $opt_g --bdg --keep-dup all --nomodel --shift -37 --extsize 73"
if [ -s "$opt_c" ] && [ -s "$BedControlSort" ]; then
	Macs2Options="$Macs2Options --control $BedControlSort"
fi

MacsMergedPeaks="$opt_D/$opt_o/${opt_o}_peaks.narrowPeak"
if [ ! -d $opt_D/$opt_o ] || [ ! -s $MacsMergedPeaks ]; then
	echo "Info: call peaks for merged BEDPE: $OutBedpe2Sort"
	echo "      CMD: macs2 callpeak $Macs2Options --treatment $OutBedpe2Sort --outdir $opt_D/$opt_o --name $opt_o > $opt_o.macs2.log 2>&1"
	macs2 callpeak $Macs2Options --treatment $OutBedpe2Sort --outdir $opt_D/$opt_o --name $opt_o > $opt_o.macs2.log 2>&1
	if [ $? -ne 0 ] && [ ! -s $MacsMergedPeaks ]; then
		echo "Info: MACS2 running error" >&2
		exit 100
	fi
else
	echo "Info: using existing PEAKs for merged BEDPE: $OutBedpe2Sort"
fi
MacsMergedPeakSort="$opt_D/$opt_o/${opt_o}_peaks.narrowPeak.sorted.bed"
if [ ! -s $MacsMergedPeakSort ]; then
	echo "Info: sorting peaks for merged BEDPE: $MacsMergedPeaks"
	echo "      CMD: sort -k1,1 -k2,2n -k3,3n $MacsMergedPeaks > $MacsMergedPeakSort"
	sort -k1,1 -k2,2n -k3,3n $MacsMergedPeaks > $MacsMergedPeakSort
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
		echo "Info: call peaks for BEDPE: ${BedArr[$BedpeNum]}"
		echo "      CMD: macs2 callpeak $Macs2Options --treatment ${BedArr[$BedpeNum]} --outdir $opt_D/${RepPfxArr[$BedpeNum]} --name ${RepPfxArr[$BedpeNum]} > $opt_o.macs2.log 2>&1"
		macs2 callpeak $Macs2Options --treatment ${BedArr[$BedpeNum]} --outdir $opt_D/${RepPfxArr[$BedpeNum]} --name ${RepPfxArr[$BedpeNum]} > $opt_o.macs2.log 2>&1
		if [ $? -ne 0 ] && [ ! -s $MacsPeaks ]; then
			echo "Error: MACS2 running error for BEDPE: ${BedArr[$BedpeNum]}" >&2
			echo "    CMD:macs2 callpeak $Macs2Options --treatment ${BedArr[$BedpeNum]} --outdir $opt_D/${RepPfxArr[$BedpeNum]} --name ${RepPfxArr[$BedpeNum]} > $opt_o.macs2.log 2>&1" >&2
			exit 100
		fi
	else
		echo "Info: using existing PEAKs for BEDPE: ${BedArr[$BedpeNum]}"
	fi
	if [ ! -s $MacsPeakSort ]; then
		echo "Info: sorting peaks for BEDPE: $MacsPeaks"
		echo "      CMD: sort -k1,1 -k2,2n -k3,3n $MacsPeaks > $MacsPeakSort"
		sort -k1,1 -k2,2n -k3,3n $MacsPeaks > $MacsPeakSort
		if [ $? -ne 0 ] || [ ! -s $MacsPeakSort ]; then
			echo "Error: Peak sorting error: ${RepPfxArr[$BedpeNum]}" >&2
			echo "     CMD: sort -k1,1 -k2,2n -k3,3n $MacsPeaks > $MacsPeakSort" >&2
			exit 100
		fi
	else
		echo "Info: using existing sorted merged Peaks: $MacsMergedPeakSort"
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

