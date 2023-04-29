Replication analysis (for SSNS vs Control)

* extract Bristol individual ID's for disease subtype and age (written into text files)
* extract SNPs at the suggestive significance level from discovery set GMMAT
* Use those SNPs to subset gnomad controls, for AF test
  * run processing_ssns_ctr.q - outputs gnomad AC for each ancestry
  * with gnomad AC results, write txt file that contains all SNP id's
  * run allele_freq_discovery.q to get allele counts for each ancestry from discovery set
  * run discovery_AC_processing.R to concatenate results - will be merged with gnomad controls for AF tests. 
* Apply AF test on SSNS control vs Gnomad control
  * remove 'null' (false positive) SNPs
  * write SNP id into text file 
* use new subset of SNPs to subset bristol data
    * run allele_freq_bristol.q to get allele counts for each ancestry
    * run bristol_AC_processing.R to concatenate results into one table (written into text file) - this will be merged with gnomad controls for AF tests. 
* Apply AF test on Bristol case vs Gnomad control
    * use Bonferroni correction 0.05/n to see how many SNPs are significant
    * extract SNP id's and P-values (for all replication snps) for LD clump - run_clump_plink.q
      * use clump Bristol p-values to get effective number of tests 
      * intersect to find significant independent regions in Bristol



