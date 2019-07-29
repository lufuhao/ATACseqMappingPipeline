# ATACseqMappingPipeline

*    This Pipeline is designed to map ATAC-seq data to large genome, for example, for wheat.

*    It splits large genome files into parts and do the mapping and then finally merge them

---

## Descriptions

    1. Fastqc, trimglore, trimmomatic to clean reads

	2. Merge Genome+CtDNA+mtDNA, Split Genome fasta, and bowtie index

	3. Bowtie mapping, merge, re-coordinate, sort and rmdup the BAMs

	4. BAM to BAMPE, then to BEDPE, MACS2 callpeaks

## ## Author:
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
