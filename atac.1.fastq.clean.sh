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
	
Version: 20200909



Requirements:
    Linux: perl, echo
    JAVA/JDK
    FuhaoBash_FastqMod
    Trimmomatic
        TRIMMOMATIC_ADAPTERS: set to the adapters folder in trimmomatic
        CLASSPATH=path/to/trimmomatic.jar:\$CLASSPATH
            Add trimmomatic.jar to CLASSPATH
    Trim_Galore



Descriptions:
  This script is used to call ATAC-seq peaks
    1. FastQC
    2. trim_galore
    3. trimmomatic



Options:
  -h    -------    Print this help message
  -1    <Fq1,Fq2>  Fastq list for R1
  -2    <Fq1,Fq2>  Fastq list for R2
  -p    <Px1,Px2>  Output prefix list
  -D    <Path>     Running Path
                      Default: $PWD
  -D1   <Path>     Running Path for FastQC
                      Default: $PWD/0.fastqc
  -D2   <Path>     Running Path for trimgalore
                      Default: $PWD/1.trim_galore
  -D3   <Path>     Running Path for trimmomatic
                      Default: $PWD/2.trimmomatic
  -t    <INT>      Number of threads for samtools merge
                      Default: 1
  -s    <1,2,3>    to run which step
                      Default: 1,2,3
  -l    <INT>      Minimum length to keep a trimmed reads
                      Default: 70
  -q    <INT>      Minimum quality to keep a trimmed reads
                      Default: 15
  -opt  <STR>      Trim_Galore option, default:
                      '--paired --gzip --output_dir ./ --quality 3 --phred33 --nextera --trim1'
  -adp  <STR>      Trimmomatic Adaptors file in full path, default:
                      '\$TRIMMOMATIC_ADAPTERS/NexteraPE-PE.fa'

Example:
  $0 -1 Fq1.R1.gz,Fq2.R1.gz -2 Fq1.R2.gz,Fq2.R2.gz \
     -p Fq1,Fq2 \
     -D /path/to/run \
     -D1 /path/to/run/0.fastqc \
     -D2 /path/to/run/1.trim_galore \
     -D3 /path/to/run/2.trimmomatic \
     -t 1 -s 1,2,3 -l 70 -q 15


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
source FuhaoBash_FastqMod
RunDir=$PWD
declare -a Fastq1Arr=()
declare -a Fastq2Arr=()
declare -a RepPfxArr=()
opt_Drun=$PWD
opt_Dfqc=''
opt_Dgal=''
opt_Dtrm=''
opt_t=1
opt_c=''
StepArr=(1 2 3)
opt_l=70
opt_q=15
opt_TGopt=""
opt_TGadp=""

#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -1) FastQR1Arr=($(echo $2 | tr ',' "\n"));shift 2;;
    -2) FastQR2Arr=($(echo $2 | tr ',' "\n"));shift 2;;
    -p) PfxArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -D) opt_Drun=$2;shift 2;;
    -D1) opt_Dfqc=$2;shift 2;;
    -D2) opt_Dgal=$2;shift 2;;
    -D3) opt_Dtrm=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -s) StepArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -l) opt_l=$2;shift 2;;
    -q) opt_q=15;shift 2;;
    -opt) opt_TGopt=$2;shift 2;;
    -adp) opt_TGadp=$2;shift 2;;
    --) shift;break;;
    -*) echo "Error: no such option $1. -h for help" > /dev/stderr;exit 1;;
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



#################### Command test ###################################

CmdExists 'FuhaoBash_FastqMod'
if [ $? -ne 0 ]; then
	echo "Error: script 'FuhaoBash_FastqMod' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'fastqc'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'fastqc' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'trim_galore'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'trim_galore' is required but not found.  Aborting..." >&2 
	exit 127
fi
    


#################### Defaults #######################################
echo "Info: parameters"

if [ -z "$opt_Dfqc" ]; then
	opt_Dfqc="$opt_Drun/0.fastqc"
fi

if [ -z "$opt_Dgal" ]; then
	opt_Dgal="$opt_Drun/1.trim_galore"
fi
if [ -z "$opt_Dtrm" ]; then
	opt_Dtrm="$opt_Drun/2.trimmomatic"
fi
step1=0
step2=0
step3=0
for IndStep in ${StepArr[@]}; do
	if [ $IndStep -gt 0 ] && [ $IndStep -lt 4 ]; then
		if [ $IndStep -eq 1 ]; then
			step1=1
		fi
		if [ $IndStep -eq 2 ]; then
			step2=1
		fi
		if [ $IndStep -eq 3 ]; then
			step3=1
		fi
	else
		echo "Error: invalid step number: $IndStep" >&2
		exit 100
	fi
done

if [ -z "$opt_TGopt" ]; then
	TrimgloreOptions="--paired --gzip --output_dir ./ --quality 3 --phred33 --nextera --length $opt_l --trim1"
else
	TrimgloreOptions="${opt_TGopt} --length ${opt_l}"
fi
if [ -z "$opt_TGadp" ]; then
	if [ -z "$TRIMMOMATIC_ADAPTERS" ] || [ ! -d "$TRIMMOMATIC_ADAPTERS" ]; then
		echo "Error: Please set TRIMMOMATIC_ADAPTERS to the folder where trimmomatics adaptors locate" >&2
		exit 100
	fi
	TrimmomaticAdapters="$TRIMMOMATIC_ADAPTERS/NexteraPE-PE.fa"
else
	TrimmomaticAdapters=$opt_TGadp
fi
if [ ! -s $TrimmomaticAdapters ]; then
	echo "Error: Trimmomatic adaptors not found: $TrimmomaticAdapters" >&2
	exit 100
fi

echo "    FastQC DIR:          $opt_Dfqc"
echo "    trim_galore DIR:     $opt_Dgal"
echo "    trimmomatic DIR:     $opt_Dtrm"
echo "    Trim_galore options: $trimgloreOptions"
echo "    Trimmomatic adaptor: $TrimmomaticAdapters"



#################### Input and Output ###############################
if [ ! -d $opt_D ]; then
	mkdir -p $opt_D
fi
if [ ! -d $opt_Dfqc ]; then
	mkdir -p $opt_Dfqc
fi
if [ ! -d $opt_Dgal ]; then
	mkdir -p $opt_Dgal
fi
if [ ! -d $opt_Dtrm ]; then
	mkdir -p $opt_Dtrm
fi



#################### Main ###########################################



cd $RunDir
for (( indnum=0; indnum < ${#PfxArr[@]}; indnum++ )); do

	cd $RunDir
	fastqR1=${FastQR1Arr[$indnum]}
	fastqR2=${FastQR2Arr[$indnum]}
	OutPrefix=${PfxArr[$indnum]}
	echo "Prefix: $OutPrefix"; echo "Prefix: $OutPrefix" >&2

	if [ ! -s $fastqR1 ] || [ ! -s $fastqR2 ]; then
		echo "Error: fq1 or fq2 error" >&2
		exit 100
	else
		echo "Fastq1: $fastqR1 fastq2: $fastqR2"
	fi
	echo $idvlib; echo $idvlib >&2

# FastQC
	if [ $step1 -eq 1 ]; then
		echo "### Step1: FastQC"
		cd $opt_Dfqc
		fastqc $fastqR1 -o ./ -f fastq -t $opt_t --noextract
		
		fastqc $fastqR2 -o ./ -f fastq -t $opt_t --noextract
		
	fi
	



###Trim_glore
	OutFq1="$opt_Dgal/$OutPrefix.R1.fq.gz"
	OutFq2="$opt_Dgal/$OutPrefix.R2.fq.gz"
	InFq1=$fastqR1
	InFq2=$fastqR2
	if [ $step2 -eq 1 ]; then
		echo "### Step2: Trim_galore"
		if [ ! -d $opt_Dgal ]; then
			mkdir -p $opt_Dgal
		fi
		cd $opt_Dgal
		if RunTrimGalore $InFq1 $InFq2 "$TrimgloreOptions" $opt_t $OutFq1 $OutFq2; then
#trim_galore -q 25 --phred33 --length 25 -e 0.1 --stringency 4 -o $analysis_dir/clean
			echo "Info: TrimGalore successful: $idvlib"
		else
			echo "Info: TrimGalore error: $idvlib" >&2;
			exit 100
		fi
	fi
	fastqR1=$OutFq1
	fastqR2=$OutFq2

### Trimmomatic
	InFq1=$OutFq1
	InFq2=$OutFq2
	OutFq1="$opt_Dtrm/$OutPrefix.R1.trim.fq.gz"
	OutFq2="$opt_Dtrm/$OutPrefix.R2.trim.fq.gz"
	if [ $step3 -eq 1 ]; then
		echo "### Step3: trimmomatic"
		if [ ! -d $opt_Dtrm ]; then
			mkdir -p $opt_Dtrm
		fi
		cd $opt_Dtrm
		if RunTrimmomatic2 $InFq1 $InFq2 $OutFq1 $OutFq2 $OutPrefix $TrimmomaticAdapters $opt_q $opt_l $opt_t; then
			echo "Info: trimmomatic successful: $idvlib"
		else
			echo "Info: trimmomatic error: $idvlib" >&2;
			exit 100
		fi
	fi
done

echo "### $ProgName finished"

exit 0
