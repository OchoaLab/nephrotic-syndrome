Replication analysis (for NS vs Control and SSNS vs Control)

* extract Bristol individual ID's for disease subtype and age (written into text files)
* extract SNPs at the suggestive significance level from discovery set SAIGE summary stats
* Use those SNPs to subset UKBB controls, for AF test
  * run **allele_freq_discovery.q** to get allele counts for each disease subtype + ancestry from discovery set
* Apply AF test on NS/SSNS control vs UKBB control
  * remove 'null' (false positive) SNPs
  * write SNP id into text file 
* use new subset of SNPs to subset bristol data
    * run **allele_freq_bristol.q** to get allele counts for each disease subtype + ancestry
* Apply AF test on Bristol case vs UKBB control
    * use Bonferroni correction 0.05/n to see how many SNPs are significant
    * extract SNP id's and P-values (for all replication snps) for LD clump - run_clump_plink.q
      * use clump Bristol p-values to get effective number of tests 
      * intersect to find significant independent regions in Bristol

Note: UKBB allele counts are extracted using scripts in the 'UKBB' subfolder. Replication results Manhattan plot figures are available in the 'plots' subfolder



