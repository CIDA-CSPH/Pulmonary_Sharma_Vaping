mkdir /beevol/home/lincuini/R_libs/
  
cat > $HOME/.Renviron
R_LIBS_USER=/beevol/home/lincuini/R_libs/

# Error: C++14 standard requested but CXX14 is not defined

mkdir $HOME/.R/
cat > $HOME/.R/Makevars  
CXX14 = g++
CXX14FLAGS = -g -O2 $(LTO)
CXX14PICFLAGS = -fpic
CXX14STD = -std=gnu++14


# mkdir Datasets
mkdir Analyses
mkdir Analyses/VapingMethylation
mkdir Analyses/VapingMethylation/
mkdir Analyses/VapingMethylation/Data
mkdir Analyses/VapingMethylation/Data/Idats
mkdir Analyses/VapingMethylation/Data/Metadata
mkdir Analyses/VapingMethylation/Scripts

qlogin -R rusage[mem=8]
module load gcc/7.4.0
module load R/4.0.3
R

install.packages("tidyverse")
install.packages("data.table")
install.packages("readxl")
install.packages("BiocManager")

BiocManager::install("minfi")
BiocManager::install("sesame")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("khroma")






mkdir /beevol/home/lincuini/software
cd  /beevol/home/lincuini/software
git clone https://github.com/brentp/combined-pvalues.git

cd combined-pvalues

qlogin -R rusage[mem=8]
module load gcc/7.4.0
module load anaconda/4.7.10
conda create -n py38 python=3.8.5 anaconda
# conda init bash

conda activate py38
pip install setuptools==33.1.1
python setup.py install --install-dir=/beevol/home/lincuini/Python_libs


cd /beevol/home/lincuini/software/combined-pvalues
python -m pip install --upgrade pip --user

python setup.py install --user


# conda setup.py install --install-lib=/beevol/home/lincuini/Python_libs


conda install toolshed
conda config --append channels conda-forge
conda install toolshed --use-local

