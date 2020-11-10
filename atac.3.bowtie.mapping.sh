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
	
Version: 20200912



Requirements:
    Linux: grep, cut, sort, echo, cat perl
    samtools
    bowtie/bowtie2
    picard
        PICARD_JAR=/path/to/picard.jar
    bam_filter_by_readname_file.pl
    bam_restore_splited_coords.pl

Descriptions:
  This script is used to map ATAC-seq fastqR1 and R2 to large genome



Options:
  -h    -------    Print this help message
  -1    <Fq1,Fq2>  FastqR1 array, separated by comma
  -2    <Fq1,Fq2>  FastqR2 array, separated by comma
  -p    <Px1,Px2>  File ouput Prefix array, separated by comma
  -x    <Ix1,Ix2>  Index array, separated by comma
  -b1              Using Bowtie for Mapping
  -b2              Using Bowtie2 for mapping [default]
  -s    <BED>      BED file to restore splited chom coordinates
  -D1   <PATH>     Running Folder/Path
  -D2   <PATH>     Folder/Path for bowtie[2] mapping BAMs
  -D3   <PATH>     Folder/Path for Cleaning BAMs
  -exu  <FILE>     Chromosome Name to be excluded
                     Mitochondria and Chloroplast genomes
                     Must exclude reads mapped to these two naked DNAs
  -t    <INT>      Num of threads;
  -d    -----      Delete temporary files

Example:
  $0 -1 In1.R1.Fq.gz,In2.R1.Fq.gz \
     -2 In1.R2.Fq.gz,In2.R2.Fq.gz \
     -p In1,In2 -x /path/Index1,/path/Index2,/path/Index3 \
     -D1 /Path/to/run \
     -D2 /Path/to/run/5.mapping \
     -D3 /Path/to/run/8.merge \
     -exu /path/to/mit.cp.names \
     -t 1 -b2
 


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

declare -a PrefixArr=()
declare -a FastqR1Arr=()
declare -a FastqR2Arr=()
declare -a IndexArr=()
RunDir=$PWD
MappingDir=''
CleanDir=''
ExcludeChromList=''
opt_t=1
opt_d=0
opt_s=''
opt_b1=0
opt_b2=1
#################### Parameters #####################################

while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -1) FastqR1Arr=($(echo $2 | tr ',' "\n"));shift 2;;
    -2) FastqR2Arr=($(echo $2 | tr ',' "\n"));shift 2;;
    -p) PrefixArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -x) IndexArr=($(echo $2 | tr ',' "\n"));shift 2;;
    -b1) opt_b1=1;opt_b2=0;shift 1;;
    -b2) opt_b1=0;opt_b2=1;shift 1;;
    -s) opt_s=$2;shift 2;;
    -D1) RunDir=$2;shift 2;;
    -D2) MappingDir=$2;shift 2;;
    -D3) CleanDir=$2; shift 2;;
    -exu) ExcludeChromList=$2; shift 2;;
    -t) opt_t=$2;shift 2;;
    -d) opt_d=1;shift 1;;
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

CmdExists 'samtools'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
samtoolsVers=$(samtools 2>&1 | grep ^'Version' | sed 's/^Version:\s\+//;s/\..*$//;')
CmdExists 'bowtie'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'bowtie' in PROGRAM 'bowtie' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'bam_filter_by_readname_file.pl'
if [ $? -ne 0 ]; then
	echo "Error: script 'bam_filter_by_readname_file.pl' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'bam_restore_splited_coords.pl'
if [ $? -ne 0 ]; then
	echo "Error: script 'bam_restore_splited_coords.pl' is required but not found.  Aborting..." >&2
	exit 127
fi
if [ -z "$PICARD_JAR" ]; then
	echo "Error: PICARD_JAR Not set, please export PICARD_JAR=/path/to/picard.jar"
	exit 127
fi



#################### Defaults #######################################
declare -a TempFiles=()
if [ -z "$RunDir" ]; then
	RunDir=$PWD;
fi
if [ -z "$MappingDir" ]; then
	MappingDir="$RunDir/mapping"
fi
if [ -z "$CleanDir" ]; then
	CleanDir="$RunDir/merge"
fi
opt_b=0
mappingProg="unMap"
if [ $opt_b1 -eq 0 ] && [ $opt_b2 -eq 1 ]; then
	echo "### mapping using bowtie2"
	opt_b=2
	mappingProg="Bowtie2"
elif [ $opt_b1 -eq 1 ] && [ $opt_b2 -eq 0 ]; then
	echo "### mapping using bowtie"
	opt_b=1
	mappingProg="Bowtie"
else
	echo "Error: do not know using bowtie or bowtie2 for mapping" >&2
	exit 100
fi


#################### Input and Output ###############################
if [ ${#PrefixArr[@]} -gt 0 ] && [ ${#PrefixArr[@]} -eq ${#FastqR1Arr[@]} ] && [ ${#PrefixArr[@]} -eq ${#FastqR2Arr[@]} ]; then
	echo "Info: Number of Read Pairs detected: ${#PrefixArr[@]}"
else
	echo "Error: Please check the number of Prefix, R1 and R2 due to the unequal number" >&2
	exit 100;
fi
if [ ${#IndexArr[@]} -gt 0 ]; then
	echo "Info: Number of Index: ${#IndexArr[@]}"
	echo "Info: Number of Index: ${#IndexArr[@]}" >&2
	for IndIdx in ${IndexArr[@]}; do
		echo "Info: Index: $IndIdx"
		echo "Info: Index: $IndIdx" >&2
	done
else
	echo "Error: Empty Index" >&2
	exit 100;
fi



#################### Main ###########################################
if [ ! -d $RunDir ]; then
	mkdir -p $RunDir;
fi
if [ ! -d $MappingDir ]; then
	mkdir -p $MappingDir;
fi
if [ ! -d $CleanDir ]; then
	mkdir -p $CleanDir;
fi
echo "Info: Running folder: $RunDir";
echo "Info: Mapping folder: $MappingDir";
echo "Info: Cleaning folder: $CleanDir";

cd $RunDir

for ((IndNum=0;IndNum<${#PrefixArr[@]};IndNum++));do

	echo "#########################START##########################"; 
	echo "#########################START##########################" >&2

	cd $RunDir
	InFq1=${FastqR1Arr[$IndNum]};
	InFq2=${FastqR2Arr[$IndNum]};
	OutPrefix=${PrefixArr[$IndNum]};
	echo "Info: Prefix: $OutPrefix"; echo "Info: Prefix: $OutPrefix" >&2
	echo "Info: $InFq1"; echo "Info: $InFq1" >&2
	echo "Info: $InFq2"; echo "Info: $InFq2" >&2


###### Mapping
	cd $MappingDir;
	mergeBAMs="";
	NumBams=0;
	BamNum=1;
	declare -a BAMarr=()
	for IndIndex in ${IndexArr[@]}; do
		OutBam="$MappingDir/$OutPrefix.$mappingProg.$BamNum.bam"
		if [ ! -s $OutBam ]; then
			if [ $opt_b -eq 1 ]; then
				bowtie -q -p $opt_t -X 2000 --fr -m 1 --sam --chunkmbs 500 $IndIndex -1 $InFq1 -2 $InFq2 2> $OutPrefix.err | samtools view -@ $opt_t -b -S -h -F 12 - > $OutBam
				if [ $? -ne 0 ]; then
					echo "Error: BOWTIE running error1: $OutPrefix" >&2
					echo "CMD used: bowtie -q -p $opt_t -X 2000 --fr -m 1 --sam --chunkmbs 500 $IndIndex -1 $InFq1 -2 $InFq2 2> $OutPrefix.err | samtools view -@ $opt_t -b -S -h -F 12 - > $OutBam" >&2
					exit 100;
				fi
			elif [ $opt_b -eq 2 ]; then
				bowtie2 -q -p $opt_t -X 2000 --fr $IndIndex -1 $InFq1 -2 $InFq2 2> $OutPrefix.err | samtools view -@ $opt_t -b -S -h -F 12 - > $OutBam
				if [ $? -ne 0 ]; then
					echo "Error: BOWTIE2 running error1: $OutPrefix" >&2
					echo "CMD used: bowtie -q -p $opt_t -X 2000 --fr $IndIndex -1 $InFq1 -2 $InFq2 2> $OutPrefix.err | samtools view -@ $opt_t -b -S -h -F 12 - > $OutBam" >&2
					exit 100;
				fi
			fi
		else
			echo "Info: using existing BAM: $OutBam"
		fi
		if [ -s $OutBam ]; then
			mergeBAMs="$mergeBAMs $OutBam"
			BAMarr+=("$OutBam");
			((NumBams++));
			TempFiles+=("$OutBam");
		else
			echo "Error: BAM not existed: $OutBam" >&2
			exit 100;
		fi
		((BamNum++));
	done

###### Cleaning
	if [ $NumBams -eq 0 ] || [ $NumBams -ne ${#IndexArr[@]} ]; then
		echo "Error: Invalid number of BAMs: Prefix $OutPrefix; Number Index: ${#IndexArr[@]}; Number BAMs: $NumBams" >&2
		exit 100
	fi

###### merge BAMs if more than 1 index
	echo "Info: $mappingProg running succeswsful: $OutPrefix" >&2
	echo "Test: BAMs to be merged:  ${#BAMarr[@]}";
	for InTestBam in ${BAMarr[@]}; do
			echo "        $InTestBam"
	done
	cd $CleanDir
	if [ $NumBams -eq 1 ]; then
		OutMerge=${BAMarr[0]}
	elif [ $NumBams -gt 1 ]; then
		cd $CleanDir
		OutMerge="$CleanDir/$OutPrefix.$mappingProg.bam"
		
		if [ ! -s $OutMerge ]; then
			perl -e 'print "\@HD\tVN:1.0\tSO:unsorted\n"' > $CleanDir/$OutPrefix.reheader
			if [ -e "$CleanDir/$OutPrefix.reheader2" ]; then
				rm -rf $CleanDir/$OutPrefix.reheader2
			fi
			for indoutbam in ${BAMarr[@]}; do
				echo "Retrieving header: BAM $indoutbam"
				samtools view -@ $opt_t -H $indoutbam | grep ^'@SQ' >> $CleanDir/$OutPrefix.reheader2
			done
			cat $CleanDir/$OutPrefix.reheader2 | sort -u | sort -k2,2 >> $CleanDir/$OutPrefix.reheader
			if [ -e "$CleanDir/$OutPrefix.reheader2" ]; then
				rm -rf $CleanDir/$OutPrefix.reheader2
			fi
			(cat $CleanDir/$OutPrefix.reheader; for indoutbam in ${BAMarr[@]}; do samtools view -@ $opt_t $indoutbam; done) | samtools view -@ $opt_t -S -b -h -F 12 - > $OutMerge
			if [ $? -ne 0 ]; then
				echo "Error: samtools merge failed: (cat $CleanDir/$OutPrefix.reheader; for indoutbam in ${BAMarr[@]}; do samtools view -@ $opt_t $indoutbam; done) | samtools view -@ $opt_t -S -b -h - > $OutMerge" >&2
				exit 100
			fi
		else
			echo "Info: using existsing merged BAM: $OutMerge"
		fi
		TempFiles+=("$OutMerge");
	fi
	if [ ! -s $OutMerge ]; then
		echo "Error: Merged BAMs not exists: Prefix $OutPrefix" >&2
		exit 100;
	fi

	cd $CleanDir
	echo "Info: merged BAMs: $OutMerge"
	OutMerge1Clean="$CleanDir/$OutPrefix.$mappingProg.clean.bam"
	OutMerge1R2E="$CleanDir/$OutPrefix.reads2exclude"
	OutMerge2Exclude="$CleanDir/$OutPrefix.$mappingProg.clean.exc.bam"
	OutMerge2reCoord="$CleanDir/$OutPrefix.$mappingProg.clean.recoord.bam"
	OutMerge3Sort="$CleanDir/$OutPrefix.$mappingProg.clean.exc.sort.bam"
	OutMerge4Reheader="$CleanDir/$OutPrefix.$mappingProg.clean.exc.sort.reheader.bam"
	OutMerge5Rmdup="$CleanDir/$OutPrefix.$mappingProg.clean.sort.exc.sort.reheader.rmdup.bam"
###### Clean reads mapped to mtDNA and ctDNA
	if [ -s $ExcludeChromList ]; then
		if [ ! -s $OutMerge1Clean ]; then
			samtools view -@ $opt_t -h -F 12 -f 2 $OutMerge | grep -v -f $ExcludeChromList | samtools view -@ $opt_t -b -h -S - > $OutMerge1Clean
			if [ $? -ne 0 ]; then
				echo "Error: paired error: $OutPrefix" >&2
				exit 100
			fi
			TempFiles+=("$OutMerge1Clean");
		else
			echo "Info: using existsing chromosome-excluded BAM: $OutMerge1Clean"
		fi
		if [ ! -s "$OutMerge1R2E" ]; then
			samtools view -@ $opt_t -F 12 -f 2 $OutMerge | grep -f $ExcludeChromList | cut -f 1 | sort -u > $OutMerge1R2E
			if [ $? -ne 0 ] || [ ! -s "$OutMerge1R2E" ]; then
				echo "Error: extract error: $OutPrefix" >&2
				exit 100
			fi
			TempFiles+=("$OutMerge1R2E");
		else
			echo "Info: using existsing reads to be excluded: $OutMerge1R2E"
		fi
		if [ ! -s $OutMerge2Exclude ]; then
			bam_filter_by_readname_file.pl $OutMerge1Clean $OutMerge1R2E 0 $OutMerge2Exclude
			if [ $? -ne 0 ] || [ ! -s $OutMerge2Exclude ]; then
				echo "Error: clean error: $OutPrefix" >&2
				exit 100
			fi
			TempFiles+=("$OutMerge2Exclude");
		else
			echo "Info: using existsing reads-excluded BAM: $OutMerge2Exclude"
		fi
	else
		echo "Info: ignored cleaning because the chromosome excluding list not provided"
		OutMerge2Exclude=$OutMerge
	fi
	if [ ! -s "$OutMerge2Exclude" ]; then
		echo "Error: BAMs after excluding reads not exists: $OutMerge2Exclude" >&2
		exit 100
	fi
### recoordinate
	if [ ! -z "$opt_s" ]; then
		echo "Info: re-calculate BAM coordinates: $OutMerge2Exclude using BED: $opt_s"
		if [ -s "$opt_s" ]; then
			if [ ! -s $OutMerge2reCoord ]; then
				bam_restore_splited_coords.pl $OutMerge2Exclude $opt_s $OutMerge2reCoord
				if [ $? -ne 0 ] || [ ! -s $OutMerge2reCoord ]; then
					echo "Error: recoordinate failed" >&2
					echo "    CMD: bam_restore_splited_coords.pl $OutMerge2Exclude $opt_s $OutMerge2reCoord"
					exit 100
				fi
			else
				echo "Info: using existsing re-coordinated BAMs: $OutMerge2reCoord"
			fi
			OutMerge2Exclude="$OutMerge2reCoord"
		else
			echo "Error: invalid BED file to recoordinate BAM: -s $opt_s" >&2
			exit 100
		fi
	fi
### Sort
	if [ ! -s $OutMerge3Sort ]; then
		if [ $samtoolsVers -eq 0 ]; then
			samtools sort -@ $opt_t $OutMerge2Exclude ${OutMerge3Sort%.bam}
			if [ $? -ne 0 ] || [ ! -s $OutMerge3Sort ]; then
				echo "Error: sort error: $OutPrefix" >&2
				exit 100
			fi
		elif [ $samtoolsVers -eq 1 ]; then
			samtools sort -@ $opt_t -O BAM -o $OutMerge3Sort $OutMerge2Exclude 
			if [ $? -ne 0 ] || [ ! -s $OutMerge3Sort ]; then
				echo "Error: sort error: $OutPrefix" >&2
				exit 100
			fi
		fi
		TempFiles+=("$OutMerge3Sort");
	else
		echo "Info: using existsing sorted BAM: $OutMerge3Sort"
	fi
###### remove @PG
	if [ ! -s $OutMerge4Reheader ]; then
		(samtools view -@ $opt_t -H $OutMerge3Sort | grep -v ^'@PG'; samtools view -@ $opt_t $OutMerge3Sort) | samtools view -@ $opt_t -h -S -b - > $OutMerge4Reheader
		if [ $? -ne 0 ] || [ ! -s $OutMerge4Reheader ]; then
			echo "Error: reheader error: $OutPrefix" >&2
			exit 100
		fi
		TempFiles+=("$OutMerge4Reheader");
	else
		echo "Info: using existsing reheader BAM: $OutMerge3Sort"
	fi
	if [ ! -s "$OutMerge4Reheader.bai" ]; then
		samtools index -@ $opt_t $OutMerge4Reheader
		TempFiles+=("$OutMerge4Reheader");
	fi
###### MarkDuplicates
	if [ ! -s $OutMerge5Rmdup ]; then
		java  -jar $PICARD_JAR MarkDuplicates INPUT=$OutMerge4Reheader OUTPUT=$OutMerge5Rmdup REMOVE_DUPLICATES=TRUE METRICS_FILE=$OutPrefix.$mappingProg.clean.sort.exc.rmdup.metrix ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
		if [ $? -ne 0 ] || [ ! -s $OutMerge5Rmdup ]; then
			echo "Error: rmdup error: $OutPrefix" >&2
			exit 100
		fi
	else
		echo "Info: using existsing deduplicated BAM: $OutMerge3Sort"
	fi

	echo "##########################END###########################"; 
	echo "##########################END###########################" >&2
done

if [ $opt_d -eq 1 ]; then
	echo "Info: Cleaning temporary files"
	for IndTemp in ${TempFiles[@]}; do
		rm -rf "$IndTemp" > /dev/null 2>/dev/null
	done
fi

echo "Info: $ProgramName successfully finished"
