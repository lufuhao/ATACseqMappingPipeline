#!/bin/bash

### source: https://www.jianshu.com/p/5bce43a537fd
### source: https://mp.weixin.qq.com/s/7wNRrpkqcuQmJ7ASlpytqw


###文章是 :The landscape of accessible chromatin in mammalian preimplantation embryos. Nature 2016 Jun 30;534(7609):652-7. PMID: 27309802
###数据GEO:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66581
###在SRA数据库可以下载原始测序数据 , 从文章找到数据的ID： https://www.ncbi.nlm.nih.gov/sra?term=SRP055881 把下面的内容保存到文件，命名为 srr.list 就可以使用prefetch这个函数来下载。



### 1. linux环境及软件安装
###这里首推conda
# https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/
# https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/ 
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc 
## 安装好conda后需要设置镜像。
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
conda config --set show_channel_urls yes

conda  create -n atac -y python=2 bwa
conda info --envs
source activate atac
# 可以用search先进行检索
conda search trim_galore
## 保证所有的软件都是安装在 wes 这个环境下面
conda install -y sra-tools  
conda install -y trim-galore  samtools bedtools
conda install -y deeptools homer  meme
conda install -y macs2 bowtie bowtie2 
conda install -y  multiqc 
conda install -y  sambamba


###ATAC-seq软件环境移植到另外一台电脑!
###首先通过activate target_env要分享的环境target_env，然后输入下面的命令会在当前工作目录下生成一个environment.yml文件，
#conda env export > environment.yml
###小伙伴拿到environment.yml文件后，将该文件放在工作目录下，可以通过以下命令从该文件创建环境
#conda env create -f environment.yml



### 2.下载数据

###前面提到的SRA数据库，该文章配套数据太多，我们节选部分作为练习，文件config.sra 如下：
#2-cell-1 SRR2927015
#2-cell-2 SRR2927016
#2-cell-5 SRR3545580
#2-cell-4 SRR2927018

### 2.1 download with prefetch in SRA-toolkit
cat srr.list |while read id;do
	nohup $prefetch $id -X 100G  &;
done

### 2.2 download with ascp in Aspera-connect
### EBI
ascp -QT -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l200m -P 33001 era-fasp@fasp.sra.ebi.ac.uk:/vol1/srr/SRR292/005/SRR2927015
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/005/SRR2927015/SRR2927015_1.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/005/SRR2927015/SRR2927015_2.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/006/SRR2927016/SRR2927016_1.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/006/SRR2927016/SRR2927016_2.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/008/SRR2927018/SRR2927018_1.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/008/SRR2927018/SRR2927018_2.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR354/000/SRR3545580/SRR3545580_1.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR354/000/SRR3545580/SRR3545580_2.fastq.gz
cat fq.txt | while read id; do
	ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  era-fasp@$id ./
done
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/005/SRR2927015/SRR2927015_1.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/005/SRR2927015/SRR2927015_2.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/006/SRR2927016/SRR2927016_1.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/006/SRR2927016/SRR2927016_2.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/008/SRR2927018/SRR2927018_1.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR292/008/SRR2927018/SRR2927018_2.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR354/000/SRR3545580/SRR3545580_1.fastq.gz
###fasp.sra.ebi.ac.uk:/vol1/fastq/SRR354/000/SRR3545580/SRR3545580_2.fastq.gz
cat fq.txt | while read id;do
	ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  era-fasp@$id ./
done
ascp -k 1 -QT -l 300m -P33001 -i ~/miniconda3/envs/rna/etc/asperaweb_id_dsa.openssh era-fasp@${id} ${outputdir}

###NCBI
###/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2927015/SRR2927015.sra
###/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2927016/SRR2927016.sra
###/sra/sra-instant/reads/ByRun/sra/SRR/SRR354/SRR3545580/SRR3545580.sra
###/sra/sra-instant/reads/ByRun/sra/SRR/SRR292/SRR2927018/SRR2927018.sra
#nohup /home/xlfang/.aspera/connect/bin/ascp -T -l 200M -P 33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -T --mode recv --host  ftp-private.ncbi.nlm.nih.gov --user anonftp --file-list ./srrlist.txt  ./ & ##3803
nohup /home/xlfang/.aspera/connect/bin/ascp -T -P 33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -l 200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR354/SRR3545580/SRR3545580.sra ./ &




### 3.组织项目
mkdir -p  ~/project/atac/
cd ~/project/atac/
mkdir {sra,raw,clean,align,peaks,motif,qc}
cd sra 
## vim 或者cat命令创建 srr.list 文件, 里面保存着作为练习使用的4个数据ID 
source activate atac 
cat srr.list |while read id;do ( nohup  prefetch $id & );done
## 默认下载目录：~/ncbi/public/sra/ 
ls -lh ~/ncbi/public/sra/
## 下载耗时，自行解决，学员使用现成数据：/public/project/epi/atac/sra 

## 假如提前下载好了数据。
cd ~/project/atac/ 
ln -s /public/project/epi/atac/sra  sra

###总之数据如下：
###-rw-r--r-- 1 stu stu 4.2G Aug 25 11:10 SRR2927015.sra
###-rw-r--r-- 1 stu stu 5.5G Aug 25 11:13 SRR2927016.sra
###-rw-r--r-- 1 stu stu 2.0G Aug 25 11:12 SRR2927018.sra
###-rw-r--r-- 1 stu stu 7.0G Aug 25 11:13 SRR3545580.sra


### 4. 下载的sra数据转换为fq格式
cd ~/project/atac/
source activate atac
dump=fastq-dump
analysis_dir=raw
mkdir -p $analysis_dir
## 下面用到的 config.sra 文件，就是上面自行制作的。

# $fastq-dump sra/SRR2927015.sra  --gzip --split-3  -A 2-cell-1 -O clean/
cat config.sra |while read id; do
	echo $id
	arr=($id)
	srr=${arr[1]}
	sample=${arr[0]}
#  测序数据的sra转fasq
	nohup $dump -A $sample -O $analysis_dir  --gzip --split-3  sra/$srr.sra & 
done 

### 如果不只是4个文件，需要使用shell脚本批处理。
cut -f 10,13 SRP055881/SraRunTable.txt|\
sed 's/Embryonic stem cell/ESC/'|sed 's/early 2-cell/e2-cell/' |\
perl -alne '{$h{$F[1]}++;print "$_-$h{$F[1]}"}' |tail -n+2|awk '{print $2"\t"$1}'> config.sra

###得到的原始fq数据如下:
###-rw-rw-r-- 1 jmzeng jmzeng 2.6G Aug 24 23:10 2-cell-1_1.fastq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 2.6G Aug 24 23:10 2-cell-1_2.fastq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 3.4G Aug 24 23:31 2-cell-2_1.fastq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 3.7G Aug 24 23:31 2-cell-2_2.fastq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 1.2G Aug 24 22:46 2-cell-4_1.fastq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 1.2G Aug 24 22:46 2-cell-4_2.fastq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 4.4G Aug 24 23:52 2-cell-5_1.fastq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 4.9G Aug 24 23:52 2-cell-5_2.fastq.gz

###conda activate rna
###fastqc
nohup fastqc -t 2 -o ./ ./*.fastq.gz &
###multiqc
multiqc ./*zip -o ./


### 5. 测序数据的质量控制
###选择trim_galore软件进行过滤，双端测序数据的代码如下:
###需要自行制作 config.raw 文件， 是3列，第一列占位用，没有意义，第二列是fq1的地址，第3列是fq2的地址。
cd ~/project/atac/
mkdir -p clean 
source activate atac  
# trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o clean/ raw/2-cell-1_1.fastq.gz raw/2-cell-1_2.fastq.gz
cat config.raw  |while read id; do
	echo $id
	arr=($id)
	fq2=${arr[2]}
	fq1=${arr[1]}
	sample=${arr[0]}
	nohup  trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o  clean  $fq1   $fq2  & 
done 
ps -ef |grep trim

###trim
for i in `ls *_1.fastq.gz`
do
i=${i/_1.fastq.gz/}
echo "trim_galore --phred33 -q 25 --length 35 --stringency 4 -e 0.1 --fastqc --paired -o ../2.clean ${i}_1.fastq.gz ${i}_2.fastq.gz" 
done > trim_galore.command
cat trim_galore.command

###得到过滤后的fq文件如下：
###-rw-rw-r-- 1 jmzeng jmzeng 2.4G Aug 25 09:35 2-cell-1_1_val_1.fq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 2.3G Aug 25 09:35 2-cell-1_2_val_2.fq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 3.1G Aug 25 10:10 2-cell-2_1_val_1.fq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 3.3G Aug 25 10:10 2-cell-2_2_val_2.fq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 1.1G Aug 25 08:52 2-cell-4_1_val_1.fq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 1.1G Aug 25 08:52 2-cell-4_2_val_2.fq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 3.7G Aug 25 10:27 2-cell-5_1_val_1.fq.gz
###-rw-rw-r-- 1 jmzeng jmzeng 3.9G Aug 25 10:27 2-cell-5_2_val_2.fq.gz

###质量控制前后都需要可视化，肯定是fastqc+multiqc，代码如下；
### FastQC
cd ~/project/atac/qc 
mkdir -p clean 
fastqc -t 5  ../clean/2-cell-*gz -o clean 
mkdir -p raw 
fastqc -t 5  ../raw/2-cell-*gz -o clean 
# https://zh.wikipedia.org/wiki/ASCII
## 还有很多其它工具，比如：
qualimap='/home/jianmingzeng/biosoft/Qualimap/qualimap_v2.2.1/qualimap'
$qualimap bamqc --java-mem-size=20G  -bam $id  -outdir ./ 



### 6. alignment
###比对需要的index，看清楚物种，根据对应的软件来构建，这里直接用bowtie2进行比对和统计比对率, 需要提前下载参考基因组然后使用命令构建索引，或者直接就下载索引文件：下载小鼠参考基因组的索引和注释文件, 这里用常用的mm10
# 索引大小为3.2GB， 不建议自己下载基因组构建，可以直接下载索引文件，代码如下：
mkdir referece && cd reference
wget -4 -q ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip
unzip mm10.zip

###解压后的索引如下：
###848M Jul  5 05:03 /public/reference/index/bowtie/mm10.1.bt2
###633M Jul  5 05:00 /public/reference/index/bowtie/mm10.2.bt2
###6.0K Jul  5 05:05 /public/reference/index/bowtie/mm10.3.bt2
###633M Jul  5 05:05 /public/reference/index/bowtie/mm10.4.bt2 
###848M Jul  5 04:52 /public/reference/index/bowtie/mm10.rev.1.bt2
###633M Jul  5 04:49 /public/reference/index/bowtie/mm10.rev.2.bt2

###单端测序数据的比对代码如下：
###首先可以对测试样本走流程，完善代码:
zcat ../clean/2-cell-1_1_val_1.fq.gz |head -10000 > test1.fq 
zcat ../clean/2-cell-1_2_val_2.fq.gz  |head -10000 > test2.fq
bowtie2 -x /public/reference/index/bowtie/mm10 -1 test1.fq  -2 test2.fq
bowtie2 -x /public/reference/index/bowtie/mm10 -1 test1.fq  -2 test2.fq  -S test.sam
bowtie2 -x /public/reference/index/bowtie/mm10 -1 test1.fq  -2 test2.fq  | samtools sort -@ 5 -O bam -o test.bam -

for i in `ls *_1_val_1.fq.gz`
do
i=${i/_1_val_1.fq.gz/}
echo "bowtie2 -p 2  --very-sensitive -X 2000 -x  ${index} -1 ${i}_1_val_1.fq.gz -2 ${i}_2_val_2.fq.gz |samtools sort  -O bam  -T ${i} -@ 2 -o - > ${i}.raw.bam"
done > mapping.command
## 建议抛弃 samtools markdup功能，避免麻烦。
## https://www.jianshu.com/p/1e6189f641db
samtools markdup -r test.bam test.samtools.rmdup.bam

conda activate rna
for i in `ls *.raw.bam`
do
i=${i/.raw.bam/}
echo " samtools index ${i}.raw.bam | bedtools bamtobed -i ${i}.raw.bam  > ${i}.raw.bed|samtools flagstat ${i}.raw.bam  > ${i}.raw.stat|sambamba markdup --overflow-list-size 600000  --tmpdir='./'  -r ${i}.raw.bam ${i}.rmdup.bam |samtools index ${i}.rmdup.bam"
done > rmdup.command




### 7. Remove chrM chrC
## ref:https://www.biostars.org/p/170294/ 
## Calculate %mtDNA:
mtReads=$(samtools idxstats  ${i}.rmdup.bam | grep 'chrM' | cut -f 3)
totalReads=$(samtools idxstats  ${i}.rmdup.bam | awk '{SUM += $3} END {print SUM}')
echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'
echo "samtools flagstat  ${i}.rmdup.bam > ${i}.rmdup.stat|samtools view  -h  -f 2 -q 30    ${i}.rmdup.bam |grep -v chrM |samtools sort  -O bam  -@ 2 -o - > ${i}.last.bam|samtools index  ${i}.last.bam|samtools flagstat  ${i}.last.bam > ${i}.last.stat|bedtools bamtobed -i ${i}.last.bam  > ${i}.bed"
## 把报错信息在谷歌搜索后，在两个网页上找到了答案。
https://github.com/samtools/samtools/issues/765
https://www.biostars.org/p/288496/


for i in `ls *.raw.bam`
do
i=${i/.raw.bam/}
mtReads=$(samtools idxstats  ${i}.rmdup.bam | grep 'chrM' | cut -f 3)
totalReads=$(samtools idxstats  ${i}.rmdup.bam | awk '{SUM += $3} END {print SUM}')
echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'
echo " samtools flagstat  ${i}.rmdup.bam > ${i}.rmdup.stat|samtools view  -h  -f 2 -q 30 ${i}.rmdup.bam   |grep -v chrM |samtools sort -O bam  -T ${i} -@ 2 -o - > ${i}.last.bam|samtools index ${i}.last.bam|samtools flagstat ${i}.last.bam > ${i}.last.stat"
done > quality.command
for i in `ls *.raw.bam`
do
i=${i/.raw.bam/}
echo "bedtools bamtobed -i ${i}.last.bam> ${i}.bed "
done > bed.command




## gatk 可以在GitHub下载
/public/biosoft/GATK/gatk-4.0.3.0/gatk  MarkDuplicates \
-I test.bam -O test.picard.rmdup.bam  --REMOVE_SEQUENCING_DUPLICATES true -M test.log 

### picards 被包装在GATK里面：
### sambamba 文档： http://lomereiter.github.io/sambamba/docs/sambamba-markdup.html
conda install -y  sambamba
sambamba --help
sambamba markdup --help
sambamba markdup -r test.bam  test.sambamba.rmdup.bam
samtools flagstat test.sambamba.rmdup.bam
samtools flagstat test.bam
## 接下来只保留两条reads要比对到同一条染色体(Proper paired) ，还有高质量的比对结果(Mapping quality>=30)
## 顺便过滤 线粒体reads
samtools view -f 2 -q 30  test.sambamba.rmdup.bam | grep -v chrM |wc
samtools view -f 2 -q 30  test.sambamba.rmdup.bam | wc
samtools view -h -f 2 -q 30  test.sambamba.rmdup.bam |grep -v chrM| samtools sort  -O bam  -@ 5 -o - > test.last.bam
bedtools bamtobed -i test.last.bam  > test.bed 
ls *.bam  |xargs -i samtools index {}

探索好了整个流程，就可以直接写批处理，代码如下：
ls /home/jmzeng/project/atac/clean/*_1.fq.gz > 1
ls /home/jmzeng/project/atac/clean/*_2.fq.gz > 2
ls /home/jmzeng/project/atac/clean/*_2.fq.gz |cut -d"/" -f 7|cut -d"_" -f 1  > 0
paste 0 1 2  > config.clean ## 供mapping使用的配置文件

cd ~/project/epi/align
## 相对目录需要理解
bowtie2_index=/public/reference/index/bowtie/mm10
## 一定要搞清楚自己的bowtie2软件安装在哪里，以及自己的索引文件在什么地方！！！
#source activate atac 
cat config.clean |while read id;
do echo $id
arr=($id)
fq2=${arr[2]}
fq1=${arr[1]}
sample=${arr[0]}
## 比对过程15分钟一个样本
bowtie2  -p 5  --very-sensitive -X 2000 -x  $bowtie2_index -1 $fq1 -2 $fq2 |samtools sort  -O bam  -@ 5 -o - > ${sample}.raw.bam 
samtools index ${sample}.raw.bam 
bedtools bamtobed -i ${sample}.raw.bam  > ${sample}.raw.bed
samtools flagstat ${sample}.raw.bam  > ${sample}.raw.stat
# https://github.com/biod/sambamba/issues/177
sambamba markdup --overflow-list-size 600000  --tmpdir='./'  -r ${sample}.raw.bam  ${sample}.rmdup.bam
samtools index   ${sample}.rmdup.bam 

## ref:https://www.biostars.org/p/170294/ 
## Calculate %mtDNA:
mtReads=$(samtools idxstats  ${sample}.rmdup.bam | grep 'chrM' | cut -f 3)
totalReads=$(samtools idxstats  ${sample}.rmdup.bam | awk '{SUM += $3} END {print SUM}')
echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'

samtools flagstat  ${sample}.rmdup.bam > ${sample}.rmdup.stat
samtools view  -h  -f 2 -q 30    ${sample}.rmdup.bam   |grep -v chrM |samtools sort  -O bam  -@ 5 -o - > ${sample}.last.bam
samtools index   ${sample}.last.bam 
samtools flagstat  ${sample}.last.bam > ${sample}.last.stat 
bedtools bamtobed -i ${sample}.last.bam  > ${sample}.bed
done 

其中bowtie2比对加入了-X 2000 参数，是最大插入片段，宽泛的插入片段范围(10-1000bp)

第一步得到的bam文件如下：
-rw-rw-r-- 1 stu stu 3.7G Aug 25 14:17 2-cell-1.bam
-rw-rw-r-- 1 stu stu 4.6G Aug 25 15:32 2-cell-2.bam
-rw-rw-r-- 1 stu stu 1.8G Aug 25 15:47 2-cell-4.bam
-rw-rw-r-- 1 stu stu 5.5G Aug 25 16:49 2-cell-5.bam

过滤后的bam文件是：
3.7G Aug 25 21:08 2-cell-1.bam
490M Aug 25 21:14 2-cell-1.last.bam
776M Aug 25 21:13 2-cell-1.rmdup.bam
4.6G Aug 25 23:51 2-cell-2.bam
678M Aug 25 23:58 2-cell-2.last.bam
1.1G Aug 25 23:57 2-cell-2.rmdup.bam
1.8G Aug 26 00:41 2-cell-4.bam
427M Aug 26 00:43 2-cell-4.last.bam
586M Aug 26 00:43 2-cell-4.rmdup.bam
5.5G Aug 26 03:26 2-cell-5.bam
523M Aug 26 03:33 2-cell-5.last.bam
899M Aug 26 03:32 2-cell-5.rmdup.bam

上述脚本的步骤都可以拆分运行，比如bam文件构建index或者转为bed的:
ls *.last.bam|xargs -i samtools index {} 
ls *.last.bam|while read id;do (bedtools bamtobed -i $id >${id%%.*}.bed) ;done
ls *.raw.bam|while read id;do (nohup bedtools bamtobed -i $id >${id%%.*}.raw.bed & ) ;done

最后得到的bed文件是：
237M Aug 26 08:00 2-cell-1.bed
338M Aug 26 08:01 2-cell-2.bed
203M Aug 26 08:01 2-cell-4.bed
254M Aug 26 08:01 2-cell-5.bed





第4步，使用macs2找peaks
# macs2 callpeak -t 2-cell-1.bed  -g mm --nomodel --shift -100 --extsize 200  -n 2-cell-1 --outdir ../peaks/
ls *.bed | while read id ;do (macs2 callpeak -t $id  -g mm --nomodel --sHit  -100 --extsize 200  -n ${id%%.*} --outdir ../peaks/) ;done 
## shell 13问
conda activate atac
cat >macs2.command
ls *.bed | while read id ;do (macs2 callpeak -t $id  -g mm --nomodel --shift  -100 --extsize 200  -n ${id%%.*} --outdir ../3.macs2/) ;done 
###参数含义
macs2软件说明书详见：https://www.jianshu.com/p/21e8c51fca23



第5步，计算插入片段长度，FRiP值，IDR计算重复情况
非冗余非线粒体能够比对的fragment、比对率、NRF、PBC1、PBC2、peak数、无核小体区NFR、TSS富集、FRiP 、IDR重复的一致性！
名词解释：https://www.encodeproject.org/data-standards/terms/
参考：https://www.encodeproject.org/atac-seq/
看 sam文件第9列，在R里面统计绘图
cmd=commandArgs(trailingOnly=TRUE); 
input=cmd[1]; output=cmd[2]; 
a=abs(as.numeric(read.table(input)[,1])); 
png(file=output);
hist(a,
main="Insertion Size distribution",
ylab="Read Count",xlab="Insert Size",
xaxt="n",
breaks=seq(0,max(a),by=10)
); 

axis(side=1,
at=seq(0,max(a),by=100),
labels=seq(0,max(a),by=100)
);

dev.off()  

有了上面的绘图R脚本就可以在批量检验bam文件进行出图。

还有NFR:https://github.com/GreenleafLab/NucleoATAC/issues/18

FRiP值的计算：fraction of reads in called peak regions
bedtools intersect -a ../new/2-cell-1.bed -b 2-cell-1_peaks.narrowPeak |wc -l
148928
wc ../new/2-cell-1.bed
5105844
wc ../new/2-cell-1.raw.bed
5105844
### 搞清楚 FRiP值具体定义：


ls *narrowPeak|while  read id;
do 
echo $id
bed=../new/$(basename $id "_peaks.narrowPeak").raw.bed
#ls  -lh $bed 
Reads=$(bedtools intersect -a $bed -b $id |wc -l|awk '{print $1}')
totalReads=$(wc -l $bed|awk '{print $1}')
echo $Reads  $totalReads 
echo '==> FRiP value:' $(bc <<< "scale=2;100*$Reads/$totalReads")'%'
done 
 

Fraction of reads in peaks (FRiP) - Fraction of all mapped reads that fall into the called peak regions, i.e. usable reads in significantly enriched peaks divided by all usable reads. In general, FRiP scores correlate positively with the number of regions. (Landt et al, Genome Research Sept. 2012, 22(9): 1813–1831)

文章其它指标：https://www.nature.com/articles/sdata2016109/tables/4

可以使用R包看不同peaks文件的overlap情况。
if(F){
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
  options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  source("http://bioconductor.org/biocLite.R") 

  library('BiocInstaller')
  biocLite("ChIPpeakAnno")
  biocLite("ChIPseeker")
  
}

library(ChIPseeker)
library(ChIPpeakAnno)
list.files('project/atac/peaks/',"*.narrowPeak")
tmp=lapply(list.files('project/atac/peaks/',"*.narrowPeak"),function(x){
  return(readPeakFile(file.path('project/atac/peaks/', x))) 
})

ol <- findOverlapsOfPeaks(tmp[[1]],tmp[[2]])
png('overlapVenn.png')
makeVennDiagram(ol)
dev.off()

也可以使用专业软件，IDR 来进行计算出来，同时考虑peaks间的overlap，和富集倍数的一致性 。
source activate atac
# 可以用search先进行检索
conda search idr
source  deactivate
## 保证所有的软件都是安装在 wes 这个环境下面
conda  create -n py3 -y   python=3 idr
conda activate py3
idr -h 
idr --samples  2-cell-1_peaks.narrowPeak 2-cell-2_peaks.narrowPeak  --plot

结果如下：
Initial parameter values: [0.10 1.00 0.20 0.50]
Final parameter values: [0.00 1.06 0.64 0.87]
Number of reported peaks - 5893/5893 (100.0%)

Number of peaks passing IDR cutoff of 0.05 - 674/5893 (11.4%)

参考：https://www.biostat.wisc.edu/~kendzior/STAT877/SLIDES/keles3.pdf
第6步，deeptools的可视化

具体仍然是见：https://mp.weixin.qq.com/s/a4qAcKE1DoukpLVV_ybobA 在ChiP-seq 讲解。

首先把bam文件转为bw文件，详情：http://www.bio-info-trainee.com/1815.html
cd  ~/project/atac/new
source activate atac
#ls  *.bam  |xargs -i samtools index {} 
ls *last.bam |while read id;do
nohup bamCoverage -p 5 --normalizeUsing CPM -b $id -o ${id%%.*}.last.bw & 
done 

cd dup 
ls  *.bam  |xargs -i samtools index {} 
ls *.bam |while read id;do
nohup bamCoverage --normalizeUsing CPM -b $id -o ${id%%.*}.rm.bw & 
done 

查看TSS附件信号强度：
## both -R and -S can accept multiple files 
mkdir -p  ~/project/atac/tss
cd   ~/project/atac/tss 
source activate atac
computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 10000 -a 10000    \
-R /public/annotation/CHIPseq/mm10/ucsc.refseq.bed  \
-S ~/project/atac/new/*.bw  \
--skipZeros  -o matrix1_test_TSS.gz  \
--outFileSortedRegions regions1_test_genes.bed

##     both plotHeatmap and plotProfile will use the output from   computeMatrix
plotHeatmap -m matrix1_test_TSS.gz  -out test_Heatmap.png
plotHeatmap -m matrix1_test_TSS.gz  -out test_Heatmap.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m matrix1_test_TSS.gz  -out test_Profile.png
plotProfile -m matrix1_test_TSS.gz  -out test_Profile.pdf --plotFileFormat pdf --perGroup --dpi 720 

### 如果要批处理 ，需要学习好linux命令。

下载 bed文件：https://genome.ucsc.edu/cgi-bin/hgTables 只需要3列坐标格式文件。

查看基因body的信号强度
source activate atac
computeMatrix scale-regions  -p 15  \
-R /public/annotation/CHIPseq/mm10/ucsc.refseq.bed  \
-S ~/project/atac/new/*.bw  \
-b 10000 -a 10000  \
--skipZeros -o matrix1_test_body.gz
plotHeatmap -m matrix1_test_body.gz  -out ExampleHeatmap1.png 

plotHeatmap -m matrix1_test_body.gz  -out test_body_Heatmap.png
plotProfile -m matrix1_test_body.gz  -out test_body_Profile.png

ngsplot也是可以的。
第7步，peaks注释

统计peak在promoter，exon，intron和intergenic区域的分布
if(F){
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
source("http://bioconductor.org/biocLite.R") 

library('BiocInstaller')
biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
biocLite("org.Mm.eg.db")
}

bedPeaksFile         = '8WG16_summits.bed'; 
bedPeaksFile
## loading packages
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
require(clusterProfiler) 
peak <- readPeakFile( bedPeaksFile )  
keepChr= !grepl('_',seqlevels(peak))
seqlevels(peak, pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Mm.eg.db") 
peakAnno_df <- as.data.frame(peakAnno)


promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter) 
# 然后查看这些peaks在所有基因的启动子附近的分布情况，热图模式
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
# 然后查看这些peaks在所有基因的启动子附近的分布情况，信号强度曲线图
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), 
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAnnoPie(peakAnno)

可以载入IGV看看效果，检测软件找到的peaks是否真的合理，还可以配合rmarkdown来出自动化报告。 https://ke.qq.com/course/274681

也可以使用其它代码进行下游分析； https://github.com/jmzeng1314/NGS-pipeline/tree/master/CHIPseq

Homer 可以做
# perl ~/miniconda3/envs/atac/share/homer-4.9.1-5/configureHomer.pl  -install mm10 
# ln -s /home/jmzeng/miniconda3/envs/chipseq/share/homer-4.9.1-5/data/genomes/ genomes
# cp  /home/jmzeng/miniconda3/envs/chipseq/share/homer-4.9.1-5/config.txt  /home/stu/miniconda3/envs/atac/share/homer-4.9.1-5/config.txt
## 保证数据库下载是OK
ls -lh  ~/miniconda3/envs/atac/share/homer-4.9.1-5/data/genomes  
mkdir -p  ~/project/atac/peaks
source activate atac
cd   ~/project/atac/peaks  
ls *.narrowPeak |while read id;
do 
echo $id
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $id >{id%%.*}.homer_peaks.tmp
annotatePeaks.pl  {id%%.*}.homer_peaks.tmp mm10  1>${id%%.*}.peakAnn.xls
  2>${id%%.*}.annLog.txt
done 

Bedtools 也可以做 ：https://bedtools.readthedocs.io/en/latest/content/tools/annotate.html
第8步，motif寻找及注释

Homer 可以做
ls -lh  ~/miniconda3/envs/atac/share/homer-4.9.1-5/data/genomes  
mkdir -p  ~/project/atac/motif
cd   ~/project/atac/motif  
source activate atac
ls ../peaks/*.narrowPeak |while read id;
do 
file=$(basename $id )
sample=${file%%.*} 
echo $sample 
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $id > ${sample}.homer_peaks.tmp
nohup findMotifsGenome.pl ${sample}.homer_peaks.tmp  mm10 ${sample}_motifDir -len 8,10,12  & 
done 

meme 也可以做 ， 获取 序列：https://github.com/jmzeng1314/NGS-pipeline/blob/master/CHIPseq/step7-peaks2sequence.R

R包，比如 motifmatchr包 也可以做。 https://bioconductor.org/packages/release/bioc/html/motifmatchr.html
第9步，差异peaks分析

diffbind 使用R包DiffBind进行chip-seq差异峰的分析-表观组-生信技能树

自行摸索R包用法
第10步，为什么我不用 esATAC

新鲜出炉的一篇文章，esATAC: an easy-to-use systematic pipeline for ATAC-seq data analysis 发表于 Bioinformatics.
# https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
source("http://bioconductor.org/biocLite.R") 

library('BiocInstaller')
biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
biocLite("org.Mm.eg.db")
biocLite("esATAC")
install.packages('idr')
library(esATAC)
printMap()

因为这个R包就是包装了前面我们讲解的多个分析步骤。
第11步，多组学整合分析

待续，高级课程制作中。

然后，如果你付费参加了我的这个课程，就应该是需要教学视频配套PPT，综述以及文献，及各大公司的结题报告，思维导图和教学示例数据的，现在就扫码加入我们的答疑微信聊天群吧！（附上你的购买截图，ID，小助手会核实，谢谢合作！）

课程购买链接：https://study.163.com/course/introduction.htm?courseId=1006013004
"小礼物走一走，来简书关注我"
生信技能树
生信技能树，生信菜鸟团，jimmy
总资产134 (约11.10元)共写了55.1W字获得2,809个赞共5,033个粉丝
全部评论8
小鸿_ed96
5楼 2018.10.11 09:52
执行过程遇到ImportError: numpy.core.multiarray failed to import，问题的原因应该是numpy版本太低，执行conda update numpy
小鸿_ed96
4楼 2018.10.03 21:42
在执行批处理代码的时候，出现syntax error: unexpected end of file，这个问题该怎么修正呢？
代码如下：
bowtie2_index=/home/jihong/reference/mm10

#source activate atac

cat config.clean |while read id;

do echo $id

arr=($id)

fq2=${arr[2]}

fq1=${arr[1]}

sample=${arr[0]}

bowtie2 -p 5 --very-sensitive -X 2000 -x $bowtie2_index -1 $fq1 -2 $fq2 |samtools sort -O bam -@ 5 -o - > ${sample}.raw.bam

samtools index ${sample}.raw.bam

bedtools bamtobed -i ${sample}.raw.bam > ${sample}.raw.bed

samtools flagstat ${sample}.raw.bam > ${sample}.raw.stat

sambamba markdup --overflow-list-size 600000 --tmpdir='./' -r ${sample}.raw.bam ${sample}.rmdup.bam

samtools index ${sample}.rmdup.bam

mtReads=$(samtools idxstats ${sample}.rmdup.bam | grep 'chrM' | cut -f 3)

totalReads=$(samtools idxstats ${sample}.rmdup.bam | awk '{SUM += $3} END {print SUM}')

echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'

samtools flagstat ${sample}.rmdup.bam > ${sample}.rmdup.stat

samtools view -h -f 2 -q 30 ${sample}.rmdup.bam |grep -v chrM |samtools sort -O bam -@ 5 -o - > ${sample}.last.bam

samtools index ${sample}.last.bam

samtools flagstat ${sample}.last.bam > ${sample}.last.stat

bedtools bamtobed -i ${sample}.last.bam > ${sample}.bed

done

