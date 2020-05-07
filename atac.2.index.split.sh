#!/bin/bash
#set -o errexit
### Set readonly variable
#readonly passwd_file=”/etc/passwd”
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
	
$0
	
Version: 20190719



Requirements:
    Linux: grep, cut, sort, echo, cat, perl
    fasta_splitter.pl
    bowtie

Descriptions:
  This script is used to split fasta and create bowtie index
    1. fasta_splitter.pl to split fasta
    2. bowtie-build to create index
    3. Bowtie-inspect to check the index



Options:
  -h    -------    Print this help message
  -i    <FASTA>    Fasta file to be splited
  -m    <Fasta>    Mitochondria Fasta
  -c    <Fasta>    Chloroplast Fasta
  -s    <INT>      Maximum size each split
                     Longest sequence in fasta <=this value
  -p    <STR>      Output BOWTIE index prefix
  -D    <Path>     Path to store index

Example:
  $0 -i ./chr1.fa -s 4000000000 -p myIndex -D ./path



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
opt_i=''
opt_m=''
opt_c=''
opt_s=4000000000
opt_p='myIndex'
opt_D=$PWD

#################### Parameters #####################################

while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -m) opt_m=$2;shift 2;;
    -c) opt_c=$2;shift 2;;
    -s) opt_s=$2;shift 2;;
    -p) opt_p=$2;shift 2;;
    -D) opt_D=$2;shift 2;;
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



#################### Command test ###################################

CmdExists 'fasta_splitter.pl'
if [ $? -ne 0 ]; then
	echo "Error: script 'fasta_splitter.pl' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'bowtie-build'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'bowtie-build' in PROGRAM 'bowtie' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'bowtie-inspect'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'bowtie-inspect' in program 'BOWtie is required but not found.  Aborting..." >&2 
	exit 127
fi
	



#################### Defaults #######################################



#################### Input and Output ###############################
if [ ! -s "$opt_i" ]; then
	echo "Error: invalid fasta file input '-i'" >&2
	exit 100
else
	echo "Info: Fasta file accepted: $opt_i"
fi
if [ $opt_s -gt 0 ] && [ $opt_s -lt 4294967296 ]; then
	echo "Info: split size accepted: $opt_s"
else
	echo "Error: invalid split size: $opt_s" >&2
	exit 100
fi
if [ ! -d $opt_D ]; then
	mkdir -p $opt_D
fi
echo "Info: Prefix $opt_p"
echo "Info: Index PATH: $opt_D"



#################### Main ###########################################


cd $opt_D
### Merge fasta first
FileToBeMerged=''
if [ ! -z "$opt_m" ]; then
	if [ -s "$opt_m" ]; then
		FileToBeMerged="$FileToBeMerged $opt_m"
	else
		echo "Error: Mitochondrial DNA was specified but not exists" >&2
		exit 100
	fi
fi
if [ ! -z "$opt_c" ]; then
	if [ -s "$opt_c" ]; then
		FileToBeMerged="$FileToBeMerged $opt_c"
	else
		echo "Error: Cloroplast DNA was specified but not exists" >&2
		exit 100
	fi
fi
if [ ! -z "$FileToBeMerged" ]; then
	cat $opt_i $FileToBeMerged > $opt_p.merged.fa
	if [ $? -ne 0 ] || [ ! -s "$opt_p.merged.fa" ]; then
		echo "Error: cat not merge genome + mtDNA + ctDNA" >&2
		exit 100
	fi
	opt_i="$opt_p.merged.fa"
fi



### Split fasta
fasta_splitter.pl  -i $opt_i -p $opt_p -l $opt_s
if [ $? -ne 0 ]; then
	echo "Error: Spliting Fasta failed" >&2
	exit 100;
fi

### Bowtie
declare -a FastaArr=($(ls $opt_p.*.fa))
echo "Info: Running Bowtie-build and bowtie-inspect"
NumFa=1;
for IndFa in ${FastaArr[@]}; do
	echo "      creating Index for fasta: $IndFa"
	bowtie-build -f $IndFa "${opt_p}${NumFa}" > ${opt_p}${NumFa}.bowtie-build.log 2>&1
	if [ $? -ne 0 ]; then
		echo "Error: bowtie-build running failed" >&2
		echo "CMD used: bowtie-build -f $IndFa ${opt_p}${NumFa}" >&2
		exit 100
	fi
	echo "      Inspecting Index for fasta: $IndFa"
	bowtie-inspect -s "${opt_p}${NumFa}" > ${opt_p}${NumFa}.bowtie-inspect.log 2>&1
	if [ $? -ne 0 ]; then
		echo "Error: bowtie-inspect running failed" >&2
		echo "CMD used: bowtie-inspect -s ${opt_p}${NumFa}" >&2
		exit 100
	fi
	((NumFa++))
done

echo "####################ALL DONE###################"
