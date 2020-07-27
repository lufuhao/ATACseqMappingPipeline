# ATACseqMappingPipeline

RRID: [SCR_017558](https://scicrunch.org/resolver/RRID:SCR_017558)

*    This Pipeline is designed to map ATAC-seq data to large genome, for example, for wheat.

*    It splits large genome files into parts and do the mapping and then finally merge them

---

## Descriptions

    1. Fastqc, trimglore, trimmomatic to clean reads

    2. Merge Genome+CtDNA+mtDNA, Split Genome fasta, and bowtie index

    3. Bowtie mapping, merge, re-coordinate, sort and rmdup the BAMs

    4. BAM to BAMPE, then to BEDPE, MACS2 callpeaks

---

## Requirements

* Linux: grep, perl, cut, sort, echo, cat, echo

* trimmomatic

    export TRIMMOMATIC_ROOT=${PROGPATH}/trimmomatic/v0.39/x86_64

    export TRIMMOMATIC_ADAPTORS=${PROGPATH}/trimmomatic/v0.39/x86_64/adaptors

    export TRIMMOMATIC_JAR=${PROGPATH}/trimmomatic/v0.39/x86_64/trimmomatic-0.39.jar

    export CLASSPATH=${PROGPATH}/trimmomatic/v0.39/x86_64/trimmomatic-0.39.jar:$CLASSPATH

* trim_galore

* picard_tools

* Bowtie

* BEDtools

* SAMtools

* bamaddrg

* FuhaoBin: 

    FuhaoBashModulesbam_filter_by_readname_file.pl

    bam_restore_splited_coords.pl

    fasta_splitter.pl

    macs2_bedpe_from_bampe.sh


## Options

    See '-h' for help

1. Trim and clean reads

2. atac.2.index.split.sh

3. atac.3.bowtie.mapping.sh

4. atac.4.macs2peaks.sh


## Examples

---

## Author:
>
>  Fu-Hao Lu
>
>  Post-Doctoral Scientist in Micheal Bevan laboratory
>
>  Cell and Developmental Department, John Innes Centre
>
>  Norwich NR4 7UH, United Kingdom
>
>  E-mail: Fu-Hao.Lu@jic.ac.uk
