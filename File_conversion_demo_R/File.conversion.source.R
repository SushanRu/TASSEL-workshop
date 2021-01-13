# functions to convert simple genotype format to hapmp or VCF
# simple format template:
# SNPID   Chrom   Position    REF   ALT   line1   line2   line3   etc.
# gentoypes are coded as AA, AB, BB, based on Ref/alt, A stands for ref, B stands for alt

######################################################################
###                         SimpleToHap                            ###
###   a function to convert simple format to hapmap      format    ###
######################################################################
# input: 
#   infile: the path and file name of the simple format file (.csv)
#   outfile: the path and file name for the vcf file
# output: hapmap file

# simple format template:
# SNPID   Chrom   Position    REF   ALT   line1   line2   line3   etc.


SimpleToHap <- function(infile, outfile){
  input.geno <- as.data.frame(read_csv(infile))
  # function to change AA, AB, BB to A, C, G, T for each row
  recode.marker.hap <- function(geno.temp){
    ref.allele = geno.temp[4]
    alt.allele = geno.temp[5]
    geno2 <- geno.temp[-c(1:5)]
    geno2[geno2 == "AA"] = str_c(ref.allele, ref.allele)
    geno2[geno2 == "BB"] = str_c(alt.allele, alt.allele)
    geno2[geno2 == "AB"] = str_c(ref.allele, alt.allele)
    geno2[is.na(geno2)] = 'NN' # missing values are coded as 'NN'
    return(geno2)
  }
  
  convert.geno <- t(apply(input.geno, 1, recode.marker.hap))
  hap <- cbind(input.geno[,c(1:5)], convert.geno) %>%
    as_tibble() %>%
    rename('rs#' = SNPID, 'chrom' = Chrom, "pos" = Position) %>%
    mutate(alleles = str_c(REF, ALT, sep = "/")) %>% mutate(strand = NA, 'assembly#' = NA, center = NA, protLSID = NA, 'assayLSID' = NA, panelLSID = NA, QCcode = NA) %>%
    select('rs#', alleles, chrom, pos, strand:QCcode, everything()) %>%
    select(-REF, -ALT)
  
  write_tsv(hap, outfile)
  
  return(hap)
}


######################################################################
###                         SimpleToVCF                            ###
###   a function to generate vcf file from simple format           ###
######################################################################
# input: 
#    infile: the path and file name of the simple format file (.csv)
#    outfile: the path and file name for the vcf file
# output: 
# vcf file
# vcf dataframe

SimpleToVCF <- function(infile, outfile){
  input.geno <- as.data.frame(read_csv(infile))
  # function to change AA, AB, BB to 1/1, 1/0, 0/0, change NA to ./. for each row
  recode.marker.to.vcf <- function(geno.temp){
    geno2 <- geno.temp[-c(1:5)]
    geno2[geno2 == "AA"] = '1/1'
    geno2[geno2 == "BB"] = '0/0'
    geno2[geno2 == "AB"] = '1/0'
    geno2[is.na(geno2)] = './.' # missing values are coded as 'NN'
    return(geno2)
  }
  
  convert.geno <- t(apply(input.geno, 1, recode.marker.to.vcf))
  vcf.format <- cbind(input.geno[,c(1:5)], convert.geno) %>%
    as_tibble() %>%
    rename(ID = SNPID, 'CHROM' = Chrom, POS = Position) %>%
    mutate(QUAL = '.', FILTER = 'PASS', INFO = '.', FORMAT = 'GT') %>%
    select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, everything())
  
  
  header <- c("##fileformat=CVFv4.0",
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
              "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">",
              "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">",
              "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">",
              "##FORMAT=<ID=PL,Number=.,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">",
              "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">",
              "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
              "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">")
  write_lines(header, outfile, sep = "\n", append = F)
  line1 <- str_c(colnames(vcf.format), collapse = '\t') %>%
    str_c('#', ., "")
  write_lines(line1, outfile, append = T)
  write_tsv(vcf.format, outfile, append = T)
  return(vcf.format)
}

