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
    gene_eqtl_effects = genes[gene]
    sample_eqtls = genotypes.loc[:, list(gene_eqtl_effects.columns)]
    imputed_expression[gene] = sample_eqtls.dot(gene_eqtl_effects.T['w_hat'])


print("Cross-validating TWAS")
# Gene data to regress on
gene_expr = pd.DataFrame(imputed_expression)
# cross-validate TWAS
twas = linear_model.LinearRegression()
crossval_scores = model_selection.cross_val_score(twas, genotypes, phenotypes['phenotye'], cv=5)
# twas.fit(gene_expr, phenotypes['phenotye'])
print("The cross validation TWAS R2 are " + str(crossval_scores))
