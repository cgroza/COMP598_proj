rm LDPRED*
ldpred coord --gf LDREF/1000G.EUR \
       --ssf simulated_data/GWAS/simulated_trait.sumstats \
       --out LDPRED.coordinated \
       --N 10000 \
       --eff_type LINREG \
       --maf 0 \
       --rs SNP_ID \
       --A1 REF --A2 ALT --chr CHR --pos POS --se SE --eff BETA --pval PVAL --reffreq REF_FRQ --ncol N

#ldpred inf --cf LDPRED.coordinated \
#       --h2 0.3 \
#       --ldr 12 \
#       --ldf LDPRED.LD \
#       --out LDPRED_SNP_weights \
#       --use-gw-h2\
#       --N 10000 \
#       # --f 0.0015 \
#       # --n-iter 200

ldpred gibbs --cf LDPRED.coordinated \
       --h2 0.3 \
       --ldr 12 \
       --ldf LDPRED.LD \
       --out LDPRED_SNP_weights \
       --use-gw-h2\
       --N 10000 \
       --f 0.0015 \
       --n-iter 200

ldpred score --gf TEST \
       --rf LDPRED_SNP_weights --rf-format LDPRED \
       --pf TEST.fam --pf-format FAM \
       --f 0.0015 --summary-file TEST_SUMMARY \
       --out TEST_SCORE


