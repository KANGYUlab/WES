# T2T-YAO Improves the Analysis of Whole Exome Sequencing Data from Chinese Samples

This repository contains up-to-date workflows of Analysis of Whole Exome Sequencing Data from Chinese Samples

# Alignment

## FastQC&MultiQC
```
nohup fastqc --outdir $FASTQC_OUTPUT --threads $CORES $FASTQ_PATH/$FASTQ_PATTERN > $FASTQC_OUTPUT/fastqc.log 2>&1 &
nohup multiqc $FASTQC_OUTPUT -o $FASTQC_OUTPUT > $FASTQC_OUTPUT/multiqc.log 2>&1 &
```
## TrimGalore
```
nohup bash -c '
for FILE in $FASTQ_PATH/*_1.fastq.gz
do
  BASE=$(basename $FILE | sed 's/_1\.fastq\.gz//')
  echo "Processing $BASE"
  mkdir -p $TRIMGALORE_OUTPUT
  F1=$FASTQ_PATH/$BASE"_1.fastq.gz"
  F2=$FASTQ_PATH/$BASE"_2.fastq.gz"
  $TRIMGALORE_COMMAND \
    --quality 30 \
    --length 50 \
    --output_dir $TRIMGALORE_OUTPUT/$BASE \
    --path_to_cutadapt $CUTADAPT_COMMAND \
    --cores 8 \
    --paired \
    --fastqc \
    --trim-n $F1 $F2
done' > $TRIMGALORE_OUTPUT/trim_galore.log 2>&1 &
```
## BWA-MEM
```
nohup sh -c '
for FASTQ_PATH in *
do
  for FILE in `ls $FASTQ_PATH/*_1.fastq.gz`
  do
     BASE=`basename $FILE | sed 's/_1\.fastq\.gz//'`
     F1=$FASTQ_PATH/$BASE"_1.fastq.gz"
     F2=$FASTQ_PATH/$BASE"_2.fastq.gz"
     RG="@RG\tID:"$BASE"\tSM:"$BASE"\tLB:WES\tPL:ILLUMINA"
     bwa mem -Y -K 100000000 -t $THREADS -R"$RG" "$BWA_INDEX" "$F1" "$F2" | \
           samtools view -bS -o"$BAM_PATH/$BASE.bam" -
  done
done
' > bwa_YAO.log 2>&1 &
```
## SortSam&MarkDuplicates&Samtools stats
```
nohup java -jar picard.jar SortSam I=$BAM_PATH O=$Sort_PATH SORT_ORDER=coordinate > $SORTSAM_OUTPUT/sortsam.log 2>&1 &
nohup java -jar picard.jar MarkDuplicates I=$BAM_PATH O=$MarkDuplicates_PATH > $MARKDUPLICATES_OUTPUT/markduplicates.log 2>&1 &
nohup samtools stats $BAM_PATH > $STATS_OUTPUT > $STATS_OUTPUT/stats.log 2>&1 &
```
## EMBOSS needle
```
for ((i=1; i<=$region_num; i++))
do
  a_line=$(sed -n "${i}p" yao_seq.txt)  
  b_line=$(sed -n "${i}p" hg38_seq.txt)
  
  needle -asequence <(echo "$a_line") -bsequence <(echo "$b_line") -gapopen 10 -gapextend 0.5 -outfile $output
  grep "^# Identity:" $output >> identity.txt
done
```
# Variant Calling
## DNAScope
```
sentieon driver -t NUMBER_THREADS -r REFERENCE -i DEDUPED_BAM \
    [--interval INTERVAL_FILE] --algo DNAscope [-d dbSNP] \
    --pcr_indel_model none --model DNASCOPE_MODEL/dnascope.model \
    TMP_VARIANT_VCF
```
##TNScope
```
sentieon driver -t NUMBER_THREADS -r REFERENCE \
  -i TUMOR_DEDUPED_BAM [-q TUMOR_RECAL_DATA.TABLE] \
  --algo TNscope --tumor_sample TUMOR_SAMPLE_NAME \
  [--dbsnp DBSNP] OUT_TN_VCF
```
## FPfilter
```
FPfilter -v vcf-file -p FPfilter_path
```
# Variant Annotation
## clinvar datebase
```
vt decompose clinvar_20231126.vcf.gz -o temp.split.vcf
perl prepare_annovar_user.pl -dbtype clinvar_preprocess2 temp.split.vcf -out temp.split2.vcf
vt normalize temp.split2.vcf \
	-r $reference_genomes \
	-o temp.norm.vcf \
	-w 2000000
perl prepare_annovar_user.pl -dbtype clinvar2 temp.norm.vcf -out yao_clinvar_20231126.txt
```
## annovar annotation
```
perl convert2annovar.pl -format vcf4 $VCF_FILE -out output
perl annotate_variation.pl -filter -dbtype clinvar_20231126 -out file_name yao $avinput_file/$file_name.avinput humandb/
```
# Liftover
## Create Chain File
```
minimap2 -cx asm5 --cs QUERY_FASTA.fa REFERENCE_FASTA.fa > PAF_FILE.paf
transanno minimap2-to-chain PAF_FILE.paf --output CHAINFILE.chain
```
## Transanno
```
transanno liftvcf -m --noswap --chain CHAINFILE.chain  -o SUCCEEDED.vcf.gz --query QUERY_FASTA.fa --reference REFERENCE_FASTA.fa --vcf INPUT_VCF.vcf.gz --fail FAILED.vcf.gz
transanno liftbed --chain CHAINFILE.chain $INPUT_bed -o SUCCEEDED.bed.gz --faile FAILED.vcf.gz 
```
