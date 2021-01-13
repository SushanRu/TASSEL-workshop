setwd("C:/Users/sushan/Desktop/workshop/SampleData/File_conversion_R_demo")
library(tidyverse)
library(stringr)
source("File.conversion.source.R")


input.filename <- "example.simple.format.csv"
out.hap.filename <- "example.hmp.txt" # output file name needs to end with hmp.txt

# convert simple format to hapmap
hap.data <- SimpleToHap(input.filename, out.hap.filename)

# convert simple format to vcf
out.vcf.filename <- "example.vcf"
vcf.data <- SimpleToVCF(input.filename, out.vcf.filename)


############################################
######## convert files in terminal #########
############################################

## convert hmp to vcf
# run_pipeline.pl -h example.hmp.txt -export TASSEL.vcf -exportType VCF

## convert vcf to hmp
#run_pipeline.pl -vcf example.vcf -export TASSEL.hmp.txt -exportType Hapmap


## impute with LD KNNi
# run_pipeline.pl -vcf TASSEL.vcf -LDKNNiImputationPlugin -endPlugin -export TASSEL.imputed.vcf -exportType VCF