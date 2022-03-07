cd /home/biostats_share/liucu/Vaping/Scripts

cat > cutadapt_prep.sh 
#!/bin/bash
FASTQC=~/.local/bin/cutadapt
PROJECTDIR=/home/biostats_share/liucu/Vaping
CUTADAPT=~/.local/bin/cutadapt
TRIMMEDFOLDER=$PROJECTDIR/trimmed_reads
# mkdir $TRIMMEDFOLDER
cd $PROJECTDIR/raw_reads

echo 'conda activate py38'
for READ1 in `ls *R1*.fastq.gz`; do

READ2=$(echo | awk -v read1file=$READ1 '{ gsub("R1", "R2", read1file); print read1file }')
OUT1=$(echo | awk -v read1file=$READ1 '{ print "trim_" read1file }')
OUT2=$(echo | awk -v read2file=$READ2 '{ print "trim_" read2file }')

echo $CUTADAPT -u 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--trim-n -j 24 -q 20 --nextseq-trim=20 -m 50 -o $TRIMMEDFOLDER/$OUT1 -p $TRIMMEDFOLDER/$OUT2 $PROJECTDIR/raw_reads/$READ1 $PROJECTDIR/raw_reads/$READ2


done




~/.local/bin/cutadapt -u 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n -j 24 -q 20 --nextseq-trim=20 -m 50 -o /home/biostats_share/liucu/ObeseAsthma_RNAseq/trimmed_reads/trim_OM_06_S19_L003_R1_001.fastq.gz -p /home/biostats_share/liucu/ObeseAsthma_RNAseq/trimmed_reads/trim_OM_06_S19_L003_R2_001.fastq.gz OM_06_S19_L003_R1_001.fastq.gz OM_06_S19_L003_R2_001.fastq.gz


sh cutadapt_prep.sh &> run_cutadapt.sh
sh run_cutadapt.sh > run_cutadapt.log

