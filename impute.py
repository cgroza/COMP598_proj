import numpy as np
import pandas as pd
import sys

# load data
tissue = os.path.join("simulated_data", "TWAS", "tissue0")
genotypes = pd.read_csv(os.path.join("simulated_data", "genotype", "mini_test_genotype.csv"))
n_samples = len(genotypes)

genes = {}
# load genes
for gene_csv in os.listdir(tissue):
    print("Loading " + gene_csv)
    eqtls = pd.read_csv(os.path.join(tissue, gene_csv), sep=',', header=0)
    eqtls = eqtls.set_index("Unnamed: 0").T
    genes[os.path.splitext(os.path.basename(gene_csv))[0]] = eqtls

# will contain imputed gene expression for each gene in each sample
imputed_expression = {}
# impute
for gene in genes:
    # baseline expression is 0
    imputed_expression[gene] = [0] * n_samples
    # add the effect of each SNP to the baseline gene expression
    for snp in genes[gene]:
        # effect of eQTL on gene expression
        eqtl_effect = genes[gene][snp]["w_hat"]
        # add the effect of eQTL weighted by genotype for each sample
        for sample in range(n_samples):
            imputed_expression[gene][sample] = imputed_expression[gene][sample] + genotypes[snp][sample] * eqtl_effect


    print("Imputed " + gene + " " + str(imputed_expression[gene]))
