# my project folder
mkdir /home/biostats_share/liucu/Vaping/Scripts
cd /home/biostats_share/liucu/Vaping/Scripts

# make a script called "runfastqc_raw.sh that does the following commands
# ------------------------------------------------------------------------
cat > runfastqc_raw.sh 
#!/bin/bash
FASTQC=/home/liucu/FastQC/fastqc
TMPFOLDER=/home/biostats_share/liucu/tmp
PROJECTDIR=/home/biostats_share/liucu/Vaping/
mkdir $PROJECTDIR/fastqc_raw

$FASTQC `ls $PROJECTDIR/raw_reads/*.fastq.gz` -o $PROJECTDIR/fastqc_raw/ --noextract -d $TMPFOLDER -t 12


sh runfastqc_raw.sh > runfastqc_raw.log



