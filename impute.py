import numpy as np
import pandas as pd
import sys

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

imputed_expression = {}
# impute
for gene in genes:
    imputed_expression[gene] = [0] * n_samples
    for snp in genes[gene]:
        eqtl_effect = genes[gene][snp]["w_hat"]
        for sample in range(n_samples):
            imputed_expression[gene][sample] = imputed_expression[gene][sample] + genotypes[snp][sample] * eqtl_effect


    print("Imputed " + gene + " " + str(imputed_expression[gene]))
