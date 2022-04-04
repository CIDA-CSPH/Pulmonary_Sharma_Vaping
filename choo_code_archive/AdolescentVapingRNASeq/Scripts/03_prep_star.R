# # install STAR
# # ------------------------------------------------------------------------------
# git clone https://github.com/alexdobin/STAR.git
# cd STAR/source
# make
# 
# # test working
# /home/liucu/STAR/source/STAR --help 

# get GENCODE "PRI" for human & extract
# ------------------------------------------------------------------------------
mkdir /home/biostats_share/liucu/GENCODE_human_v37
cd /home/biostats_share/liucu/GENCODE_human_v37
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz
gunzip *
  

# setup genome indices
# ------------------------------------------------------------------------------
cd /home/biostats_share/liucu/GENCODE_human_v37
/home/liucu/STAR/source/STAR \
--runThreadN 48 \
--runMode genomeGenerate \
--genomeDir /home/biostats_share/liucu/GENCODE_human_v37 \
--genomeFastaFiles GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile gencode.v37.annotation.gtf \
--sjdbOverhang 150