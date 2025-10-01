#!/usr/bin/env bash
# usage: sh genome_pipeline.sh [/path/to/r1.fq] [/path/to/r2.fq] [output dir] [sample identifier]

cd $3

SAMPLE=$4

# 1. quality control
if [ ! -d "01_filtered_read" ]; then
  mkdir 01_filtered_read
fi

fastp -i $1 -I $2 -o 01_filtered_read/${SAMPLE}.fastp.R1.fq.gz -O 01_filtered_read/${SAMPLE}.fastp.R2.fq.gz --json 01_filtered_read/${SAMPLE}.fastp.json -w 16
rm fastp.html 


# 2. merge paired-reads
pear -f 01_filtered_read/${SAMPLE}.fastp.R1.fq.gz -r 01_filtered_read/${SAMPLE}.fastp.R2.fq.gz -o 01_filtered_read/${SAMPLE} -j 16 2> 01_filtered_read/${SAMPLE}.pear.merge.log
rm 01_filtered_read/${SAMPLE}.pear.merge.log   
rm 01_filtered_read/${SAMPLE}.discarded.fastq    

ls 01_filtered_read/${SAMPLE}.assembled.fastq | awk '{print "pigz -p 16 "$1}' | sh
rm 01_filtered_read/${SAMPLE}.unassembled.forward.fastq 
rm 01_filtered_read/${SAMPLE}.unassembled.reverse.fastq 


# 3. assemble in single cell mode
spades.py -1 01_filtered_read/${SAMPLE}.fastp.R1.fq.gz -2 01_filtered_read/${SAMPLE}.fastp.R2.fq.gz --merged 01_filtered_read/${SAMPLE}.assembled.fastq.gz -o 02_spades --sc -t 16


# 4. taxonomic binning
if [ ! -d "03_taxon_bin" ]; then
  mkdir 03_taxon_bin
fi

DBNAME="/path/to/kraken2_database"  

kraken2 --db $DBNAME 02_spades/scaffolds.fasta \
--report 03_taxon_bin/scaffold_${SAMPLE}.kraken-report.txt --use-mpa-style \
--threads 16 > 03_taxon_bin/scaffold_${SAMPLE}.kraken-out

grep '^C' 03_taxon_bin/scaffold_${SAMPLE}.kraken-out | taxonkit lineage -i 3 | taxonkit reformat -P -i 6 | cut -f 2,7 > 03_taxon_bin/txk_${SAMPLE}.txt


# 5. refinement
python src/taxonkit2binspreader.py 03_taxon_bin/txk_${SAMPLE}.txt 03_taxon_bin/input_binspreader_${SAMPLE}.txt b
python src/get_yaml.py 01_filtered_read/${SAMPLE}.fastp.R1.fq.gz 01_filtered_read/${SAMPLE}.fastp.R2.fq.gz dataset_${SAMPLE}

gfa="02_spades/assembly_graph_with_scaffolds.gfa"
prebin="03_taxon_bin/input_binspreader_${SAMPLE}.txt"
outdir="04_binspreader"
dt="dataset_${SAMPLE}.yaml"

bin-refine $gfa $prebin $outdir -Rprop -t 16 --dataset $dt | tee binspreader-Rprop_${SAMPLE}.log

python src/scaffold_each_bin.py 04_binspreader/binning.tsv binspreader_${SAMPLE} b

for f in $(ls binspreader_${SAMPLE}_bin*txt | cut -f 1 -d '.')
do
	seqkit grep -f ${f}.txt 02_spades/scaffolds.fasta > ${f}.fa
done

rm binspreader_${SAMPLE}_bin*txt
