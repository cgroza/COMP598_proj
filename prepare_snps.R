library(tidyverse)
snps <- read_csv("simulated_data/genotype/test_genotype.csv.gz")
snps.gather <- snps %>%
  gather(key=SNP, value=Value, -X1) %>%
  mutate(Conf=1, X1 = X1+1) %>% arrange(X1, SNP)

write.table(file=gzfile("long_snps.csv.gz"), x=snps.gather, row.names=FALSE, col.names=FALSE, quote = FALSE)

