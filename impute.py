import os
import numpy as np
import pandas as pd
from sklearn import datasets, linear_model, model_selection
import sys

print("Loading data")
# load data
tissue = os.path.join("simulated_data", "TWAS", "tissue0")
genotypes = pd.read_csv(os.path.join("simulated_data", "genotype", "mini_test_genotype.csv"))
phenotypes = pd.read_csv(os.path.join("simulated_data", "phenotype", "mini_test_phenotype.csv"))
n_samples = len(genotypes)

genes = {}
# load genes
for gene_csv in os.listdir(tissue):
    print("Loading " + gene_csv)
    eqtls = pd.read_csv(os.path.join(tissue, gene_csv), sep=',', header=0)
    eqtls = eqtls.set_index("Unnamed: 0").T
    genes[os.path.splitext(os.path.basename(gene_csv))[0]] = eqtls

print("Imputing gene expression")
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

print("Cross-validating TWAS")
# Gene data to regress on
gene_expr = pd.DataFrame(imputed_expression)
# cross-validate TWAS
twas = linear_model.LinearRegression()
crossval_scores = model_selection.cross_val_score(twas, genotypes, phenotypes['phenotye'], cv=5)
# twas.fit(gene_expr, phenotypes['phenotye'])
print("The cross validation TWAS R2 are " + str(crossval_scores))
