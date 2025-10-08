#!/usr/bin/env bash
# sh transcriptome_pipeline.sh [/path/to/input] [sample_list.txt] [/path/to/output] [path/to/genome.fasta] [memory in Gb]

cd $3

# quality control
if [ ! -d "01_filtered_read" ]; then
  mkdir 01_filtered_read
fi

raw_dir=$1

for SAMPLE in $(cat $2)
do
	r1=${raw_dir}/${SAMPLE}/${SAMPLE}_raw_1.fq.gz
	r2=${raw_dir}/${SAMPLE}/${SAMPLE}_raw_2.fq.gz

	fastp -i $r1 -I $r2 -o 01_filtered_read/${SAMPLE}.fastp.R1.fq.gz -O 01_filtered_read/${SAMPLE}.fastp.R2.fq.gz --json 01_filtered_read/${SAMPLE}.fastp.json -w 16

	rm fastp.html
done


# merged paired-end reads
for SAMPLE in $(cat $2)
do
	pear -f 01_filtered_read/${SAMPLE}.fastp.R1.fq.gz -r 01_filtered_read/${SAMPLE}.fastp.R2.fq.gz -o 01_filtered_read/${SAMPLE} -j 16 2> 01_filtered_read/${SAMPLE}.pear.merge.log
	rm 01_filtered_read/${SAMPLE}.pear.merge.log    
	rm 01_filtered_read/${SAMPLE}.discarded.fastq    

	ls 01_filtered_read/${SAMPLE}.assembled.fastq | awk '{print "pigz -p 16 "$1}' | sh
	ls 01_filtered_read/${SAMPLE}.unassembled.forward.fastq | awk '{print "pigz -p 16 "$1}' | sh
	ls 01_filtered_read/${SAMPLE}.unassembled.reverse.fastq | awk '{print "pigz -p 16 "$1}' | sh
done


# cat reads from all samples
zcat 01_filtered_read/*fastp.R1.fq.gz | pigz > SampleMerge.forward.fq.gz
zcat 01_filtered_read/*fastp.R2.fq.gz | pigz > SampleMerge.reverse.fq.gz
zcat 01_filtered_read/*assembled.fastq.gz | pigz > SampleMerge.assembled.fastq.gz


# co-assembly
if [ ! -d "02_rnaspades_out" ]; then
  mkdir 02_rnaspades_out
fi

rnaspades.py -1 SampleMerge.forward.fq.gz -2 SampleMerge.reverse.fq.gz --merged SampleMerge.assembled.fastq.gz -o 02_rnaspades_out/rnaspades_outs -t 16 --trusted-contigs $4 -m $5


# visualize basic assembly quality
if [ ! -d "02_rnaspades_out/rnaspades_outs/checkm" ]; then
  mkdir -p 02_rnaspades_out/rnaspades_outs/checkm 
fi

cp 02_rnaspades_out/rnaspades_outs/transcripts.fasta 02_rnaspades_out/rnaspades_outs/checkm/transcripts.fasta

checkm nx_plot 02_rnaspades_out/rnaspades_outs/checkm 02_rnaspades_out/rnaspades_outs/checkm -x fasta 
checkm len_hist 02_rnaspades_out/rnaspades_outs/checkm 02_rnaspades_out/rnaspades_outs/checkm -x fasta


# gene prediction and comprehensive annotation
if [ ! -d "03_function_annotation" ]; then
  mkdir 03_function_annotation
fi

bakta --db /path/to/bakta/db 02_rnaspades_out/rnaspades_outs/transcripts.fasta --output 03_function_annotation/bakta_out --threads 16 --skip-plot


# eggnog-mapper, further functional annotation for proteins
if [ ! -d "03_function_annotation/emapper_out" ]; then
  mkdir 03_function_annotation/emapper_out
fi

emapper.py -i 03_function_annotation/bakta_out/transcripts.faa --output 03_function_annotation/emapper_out/eggnog -m diamond --cpu 16 -d /path/to/eggnog-mapper/data


# kegg KO, further functional annotation for proteins 
if [ ! -d "03_function_annotation/kegg_out" ]; then
  mkdir 03_function_annotation/kegg_out
fi

exec_annotation -f mapper -o 03_function_annotation/kegg_out/KEGG.res 03_function_annotation/bakta_out/transcripts.faa -p /path/to/database/kofam/profiles -k /path/to/database/kofam/ko_list --cpu 16 --tmp-dir 03_function_annotation/kegg_out/tmp


# DeepFRI, deep learning-based annotation for proteins
if [ ! -d "03_function_annotation/deepfri_out" ]; then
  mkdir 03_function_annotation/deepfri_out
fi

python /path/to/DeepFRI/predict.py --fasta_fn 03_function_annotation/bakta_out/transcripts.faa -ont mf -o 03_function_annotation/deepfri_out/deepfri -v
python /path/to/DeepFRI/predict.py --fasta_fn 03_function_annotation/bakta_out/transcripts.faa -ont bp -o 03_function_annotation/deepfri_out/deepfri -v
python /path/to/DeepFRI/predict.py --fasta_fn 03_function_annotation/bakta_out/transcripts.faa -ont cc -o 03_function_annotation/deepfri_out/deepfri -v
python /path/to/DeepFRI/predict.py --fasta_fn 03_function_annotation/bakta_out/transcripts.faa -ont ec -o 03_function_annotation/deepfri_out/deepfri -v



# align reads to genes
if [ ! -d "04_FPKM" ]; then
  mkdir 04_FPKM
fi

bowtie2-build 03_function_annotation/bakta_out/transcripts.ffn 03_function_annotation/bakta_out/transcripts.ffn
grep "ribosomal" 03_function_annotation/bakta_out/transcripts.ffn | awk '{print $1"\t1\t99999999"}' | sed 's/>//g' > 03_function_annotation/bakta_out/ribosomal.bed
bioawk -cfastx '{print $name"\t"length($seq)}' 03_function_annotation/bakta_out/transcripts.ffn > 04_FPKM/rnaspades_bakta_transcripts.ffn.length

for SAMPLE in $(cat $2)
do
	bowtie2 -x 03_function_annotation/bakta_out/transcripts.ffn -1 01_filtered_read/${SAMPLE}.fastp.R1.fq.gz -2 01_filtered_read/${SAMPLE}.fastp.R2.fq.gz --sensitive --threads 16 --no-unal 2> 04_FPKM/${SAMPLE}.bowtie.log | samtools sort -O BAM -o 04_FPKM/${SAMPLE}.sort.bam
	samtools view -R 03_function_annotation/bakta_out/ribosomal.bed 04_FPKM/${SAMPLE}.sort.bam | wc -l | awk '{print "'${SAMPLE}', "$1}' > 04_FPKM/${SAMPLE}_ribosomal.count.csv   
	samtools view 04_FPKM/${SAMPLE}.sort.bam | awk '{print $3}' | uniq -c | awk '{print $2"\t"$1}' > 04_FPKM/${SAMPLE}.sort.bam.count
	csvtk join 04_FPKM/${SAMPLE}.sort.bam.count 04_FPKM/rnaspades_bakta_transcripts.ffn.length -f 1 -t | awk '{frkm=$2*150/($3*1000);print $0"\t"frkm}' | sort -k4 -n -r >> 04_FPKM/${SAMPLE}.bakta.fpkm.tsv 
	awk -F "\t" 'BEGIN{print "tag_id\tread_count\tgene_length\tFRKM\tcontig_id\tType\tStart\tStop\tStrand\ttag_id\tGene\tProduct_DbXrefs" }
	{printf("%s",$0"\t");c="grep "$1" 03_function_annotation/bakta_out/transcripts.tsv";system(c)}' 04_FPKM/${SAMPLE}.bakta.fpkm.tsv > 04_FPKM/${SAMPLE}.fpkm.annotation.tsv
	awk 'BEGIN{print "tag_id\tread_count\tgene_length\tFRKM\tquery\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_koKEGG_PathwayKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITEKEGG_TCCAZy\tBiGG_Reaction\tPFAMs" }
	{printf("%s",$0"\t");c="grep "$1" 03_function_annotation/emapper_out/eggnog.emapper.annotations";system(c)}' 04_FPKM/${SAMPLE}.bakta.fpkm.tsv > 04_FPKM/${SAMPLE}.GO.fpkm.tsv
done



# align reads to proteins
diamond makedb --in 03_function_annotation/bakta_out/transcripts.faa --db 03_function_annotation/bakta_out/transcripts.faa

for SAMPLE in $(cat $2)
do

	diamond blastx -d 03_function_annotation/bakta_out/transcripts.faa --query 01_filtered_read/${SAMPLE}.assembled.fastq.gz --sensitive --threads 16 --out 04_FPKM/${SAMPLE}.assembled.fastq.gz.diamond.out
	diamond blastx -d 03_function_annotation/bakta_out/transcripts.faa --query 01_filtered_read/${SAMPLE}.unassembled.forward.fastq.gz --sensitive --threads 16 --out 04_FPKM/${SAMPLE}.unassembled.forward.fastq.gz.diamond.out
	diamond blastx -d 03_function_annotation/bakta_out/transcripts.faa --query 01_filtered_read/${SAMPLE}.unassembled.reverse.fastq.gz --sensitive --threads 16 --out 04_FPKM/${SAMPLE}.unassembled.reverse.fastq.gz.diamond.out
	awk '{if($3> 80 && h[$1]==""){h[$1]=$2;print $2}}' 04_FPKM/${SAMPLE}.assembled.fastq.gz.diamond.out 04_FPKM/${SAMPLE}.unassembled.forward.fastq.gz.diamond.out 04_FPKM/${SAMPLE}.unassembled.reverse.fastq.gz.diamond.out | sort | uniq -c | sort -k1 -n -r | awk '{print $2"," $1}' > 04_FPKM/${SAMPLE}.all.diamond.count.csv

done


# end info
echo 'Pipeline finished successfully!'
