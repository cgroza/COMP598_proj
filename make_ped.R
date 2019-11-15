library(tidyverse)
library(snpStats)

phen <- read_csv("simulated_data/phenotype/test_phenotype.csv")
gwas <- read_tsv("simulated_data/GWAS/simulated_trait.sumstats") %>% arrange(SNP_ID)



genotype <- read.snps.long(files = "simulated_data/genotype/long_snps.csv.gz",
                          sample.id = paste(1:5000),
                          snp.id = gwas$SNP_ID,
                          verbose = TRUE, codes = c("0","1", "2"),
                          in.order = FALSE)

write.plink(file.base = "TEST", snp.major = TRUE,
            id = 1:5000,
            father = rep(0, 5000),
            mother = rep(0, 5000),
            sex = rep(1, 5000),
            snps=genotype,
            phenotype=phen$phenotye,
            position = gwas$POS, chromosome = gwas$CHR,
            allele.1 = gwas$REF, allele.2 = gwas$ALT, na.code = 0)
