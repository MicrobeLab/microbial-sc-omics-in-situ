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


# get polyA at the contig, gene, and read levels
if [ ! -d "05_polyA" ]; then
  mkdir 05_polyA
fi

bioawk -c fastx '{tmp=0;printf("%s",$name);while(1){ t=index(substr($seq,tmp) ,"AAA");printf("%s",","tmp); if(t=="0"){break} tmp=tmp+t+4} ; t=index(substr($seq,tmp) ,"AAA"); printf("%s","\n")   }' 03_function_annotation/bakta_out/transcripts.fna > 05_polyA/contig.poly_position.csv
bioawk -c fastx '{tmp=0;printf("%s",$name);while(1){ t=index(substr($seq,tmp) ,"AAA");printf("%s",","tmp); if(t=="0"){break} tmp=tmp+t+4} ; t=index(substr($seq,tmp) ,"AAA"); printf("%s","\n")   }' 03_function_annotation/bakta_out/transcripts.ffn > 05_polyA/gene.poly_position.csv

awk -F "," '{for(t=NF;t>3;t--){print $t-$(t-1)}}' 05_polyA/contig.poly_position.csv > 05_polyA/contig.poly_position.csv.diff.txt
awk -F "," '{for(t=NF;t>3;t--){print $t-$(t-1)}}' 05_polyA/gene.poly_position.csv > 05_polyA/gene.poly_position.csv.diff.txt

csvtk plot hist 05_polyA/contig.poly_position.csv.diff.txt --bins 300 --xlab polyAAA_diff_bp --ylab count --title "contig polyA position diff bp" --format svg > 05_polyA/contig.poly_position.csv.diff.svg
csvtk plot hist 05_polyA/gene.poly_position.csv.diff.txt --bins 300 --xlab polyAAA_diff_bp --ylab count --title "annotated region polyA position diff bp" --format svg > 05_polyA/gene.poly_position.csv.diff.svg


for SAMPLE in $(cat $2)
do

	bioawk -c fastx '{tmp=0;printf("%s",$name);while(1){ t=index(substr($seq,tmp) ,"AAA");printf("%s",","tmp); if(t=="0"){break} tmp=tmp+t+4} ; t=index(substr($seq,tmp) ,"AAA"); printf("%s","\n")   }' 01_filtered_read/${SAMPLE}.fastp.R1.fq.gz | awk -F "," '{if(NF>2){print $0}}' > 05_polyA/${SAMPLE}.R1.poly_position.csv
	bioawk -c fastx '{tmp=0;printf("%s",$name);while(1){ t=index(substr($seq,tmp) ,"AAA");printf("%s",","tmp); if(t=="0"){break} tmp=tmp+t+4} ; t=index(substr($seq,tmp) ,"AAA"); printf("%s","\n")   }' 01_filtered_read/${SAMPLE}.fastp.R2.fq.gz | awk -F "," '{if(NF>2){print $0}}' > 05_polyA/${SAMPLE}.R2.poly_position.csv

	awk -F "," '{for(t=NF;t>3;t--){print $t-$(t-1)}}' 05_polyA/${SAMPLE}.R1.poly_position.csv > 05_polyA/${SAMPLE}.R1.poly_position.csv.diff.txt
	awk -F "," '{for(t=NF;t>3;t--){print $t-$(t-1)}}' 05_polyA/${SAMPLE}.R2.poly_position.csv > 05_polyA/${SAMPLE}.R2.poly_position.csv.diff.txt
	
	csvtk plot hist 05_polyA/${SAMPLE}.R1.poly_position.csv.diff.txt --bins 150 --xlab polyAAA_diff_bp --ylab count --title "R1 polyA position diff bp" --format svg > 05_polyA/${SAMPLE}_R1.poly_position.csv.diff.svg
	csvtk plot hist 05_polyA/${SAMPLE}.R2.poly_position.csv.diff.txt --bins 150 --xlab polyAAA_diff_bp --ylab count --title "R2 polyA position diff bp" --format svg > 05_polyA/${SAMPLE}_R2.poly_position.csv.diff.svg

done


# further play with polyA at contig and gene levels
bioawk -cfastx '{print $name"\t"length($seq)}' 03_function_annotation/bakta_out/transcripts.fna > 05_polyA/contig.length
bioawk -cfastx '{print $name"\t"length($seq)}' 03_function_annotation/bakta_out/transcripts.ffn > 05_polyA/gene.length

awk -F "," '{print $1"\t"NF-2}' 05_polyA/contig.poly_position.csv > 05_polyA/contig.poly.count
awk -F "," '{print $1"\t"NF-2}' 05_polyA/gene.poly_position.csv > 05_polyA/gene.poly.count

csvtk join 05_polyA/contig.length 05_polyA/contig.poly.count -f 1 -t | awk -F "\t" 'BEGIN{print "length,count" }{print $2 "," $3}' > 05_polyA/contig.poly_position.csv.length_count.csv
csvtk join 05_polyA/gene.length 05_polyA/gene.poly.count -f 1 -t | awk -F "\t" 'BEGIN{print "length,count" }{print $2 "," $3}' > 05_polyA/gene.poly_position.csv.length_count.csv

csvtk plot line 05_polyA/contig.poly_position.csv.length_count.csv --xlab contig_length --ylab polyA_occurance_count --title "contig_length vs polyA_occurance_count" --scatter -x length -y count --format svg > 05_polyA/contig.poly_position.csv.length_count.svg
csvtk plot line 05_polyA/gene.poly_position.csv.length_count.csv --xlab gene_length --ylab polyA_occurance_count --title "gene_length vs polyA_occurance_count" --scatter -x length -y count --format svg > 05_polyA/gene.poly_position.csv.length_count.svg

cp src/select_gene.sh $3
grep -v "#"  03_function_annotation/bakta_out/transcripts.gff3 | grep -v region | awk -F "\t" '{if(NF>1){gsub(" ","_",$NF);print "sh select_gene.sh  05_polyA/contig.poly_position.csv "$1" "$4" "$5" \""$NF"\""}}' | sh | grep '^ID' > 05_polyA/anno_contig.polyAAA.position.csv
awk 'BEGIN{print "type\tannotation\tpos"}{if(NF>4){ printf("%s","start\t"$1"\t"$5-$3"\n"); printf("%s","end\t"$1"\t"$4-$(NF-2)"\n");} printf("%s","before\t"$1"\t"$3-$(NF-1)"\n");printf("%s","after\t"$1"\t"$NF-$4"\n"); }' 05_polyA/anno_contig.polyAAA.position.csv |  awk '{if($NF>=0){print $0}}' > 05_polyA/anno_contig.polyAAA.position.csv.start_end_before_after.tsv

csvtk plot box -t  05_polyA/anno_contig.polyAAA.position.csv.start_end_before_after.tsv -g "type" -f "pos" --title "compare" --format svg > 05_polyA/anno_contig.polyAAA.position.csv.start_end_before_after.svg

awk 'BEGIN{print "length,0-10%,10-20%,20-30%,30-40%,40-50%,50-60%,60-70%,70-80%,80-90%,90-100%"}
{leng=$4-$3;printf("%s",leng); for(t=0;t<10;t++){count[t]=0}  for(t=5;t<NF-1;t++){ count[int(10*($t-$3)/leng)]++} for(t=0;t<10;t++){printf("%s", "," count[t])} printf("%s","\n")}' 05_polyA/anno_contig.polyAAA.position.csv > 05_polyA/count_in_gene.csv

awk -F "," 'BEGIN{print "length,all,0-10%,10-20%,20-30%,30-40%,40-50%,50-60%,60-70%,70-80%,80-90%,90-100%,sd"}
{
	if($1!="length"){
	printf("%s",$1);s=0;sd=0;
	for(t=2;t<=NF;t++){s=s+$t;c[t]=$t}
	printf("%s",","s);
	if(s==0){
		for(t=2;t<=NF;t++){printf("%s",",0")}
		printf("%s","\n")
	}
	else{
		for(t=2;t<=NF;t++){printf("%s",","c[t]/s); sd=sd+(c[t]/s-0.1)*(c[t]/s-0.1)  }
		printf("%s",","sqrt(sd)"\n")
	}
	}
}' 05_polyA/count_in_gene.csv | awk -F "," '{if($2>0){print $0}} ' > 05_polyA/per_in_gene.csv


csvtk plot line --scatter 05_polyA/per_in_gene.csv --data-field-x length --data-field-y sd --title "length to sd" --format svg -o 05_polyA/per_in_gene.length_sd.svg
csvtk plot line --scatter 05_polyA/per_in_gene.csv --data-field-x all --data-field-y sd --title "total num to sd" --format svg -o 05_polyA/per_in_gene.polyAcount_sd.svg
grep -v all 05_polyA/per_in_gene.csv | awk -F "," 'BEGIN{print "realtive_position,per"}{for(t=3;t<NF;t++){print t-3","$t}}' | csvtk plot box -g "realtive_position" -f "per" --title "position to percent" --format svg > 05_polyA/per_in_gene.hist.svg


# end info
echo 'Pipeline finished successfully!'


