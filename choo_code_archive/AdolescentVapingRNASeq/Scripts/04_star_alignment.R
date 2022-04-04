cd /home/biostats_share/liucu/Vaping/Scripts

cat > 04_prep_star_alignment_cmds.sh 
#!/bin/bash
STAR=/home/liucu/STAR/source/STAR
PROJECTDIR=/home/biostats_share/liucu/Vaping
TRIMMEDFOLDER=$PROJECTDIR/trimmed_reads
ALIGNEDFOLDER=$PROJECTDIR/aligned_reads

for R1trimmed in `ls $TRIMMEDFOLDER/*R1*.fastq.gz`; do


R2trimmed=$(echo | awk -v read1file=$READ1 '{ gsub("R1", "R2", read1file); print read1file }')
OUTNAME=$(echo $R1trimmed | awk -F 'trim_|_L003_' '{print $2 }')

echo $STAR --runThreadN 12 --outFilterScoreMinOverLread 0.33 \
  --outFilterMatchNminOverLread 0.33 --outSAMtype BAM Unsorted \
  --quantMode GeneCounts --genomeDir /home/biostats_share/liucu/GENCODE_human_v37 \
  --readFilesCommand zcat --readFilesIn \
  $R1trimmed $R2trimmed --outFileNamePrefix  $ALIGNEDFOLDER/${OUTNAME}_

done

sh 04_prep_star_alignment_cmds.sh > 04_run_star_alignment.sh

mkdir ../aligned_reads
sh 04_run_star_alignment.sh > 04_run_star_alignment.log

cd ../aligned_reads
multiqc . -o multiqc_star_alignment