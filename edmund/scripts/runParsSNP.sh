#! /usr/bin/zsh

# This fiel runs parsSNP. The install paths are hard-coded for my machine
# so they need to be updated if the script is used anywhere else.
# Annovar is also here:
SNPHOME="/mnt/f/Data/University/Thesis/Tools/ParsSNP"

HERE=$(pwd)

perl $SNPHOME/annovar/convert2annovar.pl $1 --format vcf4 --outfile $SNPHOME/testrun.txt &> /dev/null
perl $SNPHOME/annovar/table_annovar.pl $SNPHOME/testrun.txt $SNPHOME/annovar/humandb/ -buildver hg19 -out $SNPHOME/testrun_anno -remove -protocol refGene,ljb26_all -operation g,f -nastring NA -otherinfo &> /dev/null
# I need to move to this wd as parsSNP is pretty limited
cd $SNPHOME
Rscript ParsSNP_application.r testrun_anno.hg19_multianno.txt &> /dev/null
cd $HERE
