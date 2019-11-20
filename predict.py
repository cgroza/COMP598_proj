
import os
import numpy as np
import pandas as pd
from sklearn import datasets, linear_model, model_selection
import sys
import pickle

os.environ["OMP_NUM_THREADS"] = "10" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "10" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "10" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "10" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "10" # export NUMEXPR_NUM_THREADS=6

print("Loading data")
# load data
tissue = os.path.join("simulated_data", "eQTL", "tissue1")
genotypes = pd.read_csv(os.path.join("simulated_data", "genotype", "test_genotype.csv.gz"))
phenotypes = pd.read_csv(os.path.join("simulated_data", "phenotype", "test_phenotype.csv"))
n_samples = len(genotypes)

genes = {}
# load genes
for gene_csv in os.listdir(tissue):
    print("Loading " + gene_csv)
    eqtls = pd.read_csv(os.path.join(tissue, gene_csv), sep=',', header=0)
    eqtls = eqtls.set_index("Unnamed: 0").T
    genes[os.path.splitext(os.path.basename(gene_csv))[0]] = eqtls

print("Imputing test gene expression")
# will contain imputed gene expression for each gene in each sample
imputed_expression = {}
# impute
for gene in genes:
    gene_eqtl_effects = genes[gene]
    sample_eqtls = genotypes.loc[:, list(gene_eqtl_effects.columns)]
    imputed_expression[gene] = sample_eqtls.dot(gene_eqtl_effects.T['w_hat'])


print("Predicting TWAS-PRS")
# Gene data to regress on
gene_expr = pd.DataFrame(imputed_expression)
gene_expr.to_csv("test_gene_exp.csv")
# cross-validate TWAS
twas = pickle.load(open("fitted_twas.pickle", 'rb'))
print("The test TWAS R2 is " + str(twas.score(genotypes, phenotypes['phenotye'])))
