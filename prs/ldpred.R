# start interactive shell
# cd /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/prs
# srun --mem 16G -p ochoalab --account ochoalab --pty bash -i
# module load R/4.1.1-rhel8
# R

# based on this tutorial
# https://privefl.github.io/bigsnpr/articles/LDpred2.html

library(bigsnpr)
library(genio)
library(ochoalabtools)
library(ggplot2)
library(dplyr)
library(readr)

# run from here
setwd('/datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/prs')
# name of data to predict on
name <- 'data'
file_phen <- '/datacommons/ochoalab/ssns_gwas/array/patient-data.txt.gz'
# for a better effective sample size estimate
file_phen_base <- '/datacommons/ochoalab/ssns_gwas/imputed/patient-data.txt.gz'
# base summary stats
#file_sumstats <- '/datacommons/ochoalab/ssns_gwas/imputed/ssns_srns/mac20-glmm-score.txt'
file_sumstats <- '/datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20-glmm-score.txt'
# file to subset to clean array SNPs (which passed QC already)
name_array <- '/datacommons/ochoalab/ssns_gwas/array/ssns_tgp_merge_clean'
# output data to reload later
file_betas_filtered <- 'betas-ssns_ctrl-array.txt.gz'

# Tiffany said this PCs file is updated
# '/datacommons/ochoalab/ssns_gwas/replication/bristol_data/GMMAT/ssns_srns/ssns_srns_mac20.eigenvec'

# load validation dataset
# generate this if it doesn't exist
rds <- paste0( name, '.rds' )
if ( !file.exists( rds ) )
    # this generates .bk and .rds files
    snp_readBed( paste0( name, '.bed' ) )
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( rds )

# read phen, only file with actual trait
# load patient data, which excludes TGP info
data <- read_tsv( file_phen, show_col_types = FALSE )

stopifnot( all( obj.bigSNP$fam$sample.ID %in% data$id ) )

# subset and reorder
indexes <- match( obj.bigSNP$fam$sample.ID, data$id )
# reorder phen
data <- data[ indexes, ]
stopifnot( all( obj.bigSNP$fam$sample.ID == data$id ) )

# make numerical version
data$pheno <- NA
data$pheno[ data$diagnosis == 'SRNS' ] <- 1
data$pheno[ data$diagnosis == 'SSNS' ] <- 0

## table( phen$pheno ) # oldest stuff, which is obsolete; it appears to have some bug, all other versions agree with each other but not this one
##   0   1 
## 327 163
## table( phen$diagnosis ) # file Tiffany actually wanted me to use most recently
## NS UNCLASSIFIED             SNS            SRNS            SSNS 
##              70               1             149             364 
## table( data$diagnosis ) # my re-cleaned up data, after subsetting
## NS UNCLASSIFIED            SRNS            SSNS 
##              70             149             365 
## table( data$pheno ) # numerical version of above
##   0   1 
## 365 149 

# transfer to this other object
obj.bigSNP$fam$affection <- data$pheno
# save aligned pheno for later loading
write_lines( data$pheno, 'pheno.txt.gz' )

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection
#NCORES <- nb_cores() # on DCC this is 1 with default settings, which is probably fine because ours is a small dataset
NCORES <- 1 # was needed locally because of automatic blas multithreading leading to two levels of multithreading, something like that
#NCORES <- 20 # for DCC run for precalculating LD matrix
message( NCORES )

# rename columns to what snp_match wants
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))

# load summary statistics (big file because it's from imputed data)
sumstats <- bigreadr::fread2( file_sumstats )
# apply Debo's transformation
sumstats$beta <- sumstats$SCORE / sumstats$VAR
sumstats$beta_se <- 1 / sqrt( sumstats$VAR )
# OLD: exclude "cM", "converged"
#sumstats <- setNames( sumstats[ -c(3, 12) ], c('chr', 'rsid', 'pos', 'a1', 'a0', 'n_eff', 'af', 'beta', 'beta_se', 'p') )
# exclude "MISSRATE" (there's no missingness here), SCORE, VAR
stopifnot( all( sumstats$MISSRATE == 0 ) )
# NOTE: a1, a0 could be reversed!
#sumstats <- setNames( sumstats[ -c(7, 9, 10) ], c('rsid', 'chr', 'pos', 'a1', 'a0', 'n_eff', 'af', 'p', 'beta', 'beta_se') )
sumstats <- setNames( sumstats[ -c(7, 9, 10) ], c('rsid', 'chr', 'pos', 'a0', 'a1', 'n_eff', 'af', 'p', 'beta', 'beta_se') )

# read base gwas data, to calculate better effective sample size
data_base <- read_tsv( file_phen_base, show_col_types = FALSE )
counts <- table( data_base$ssns_ctrl ) # NOTE: for ssns_ctrl!  Adapt if needed for ssns_srns!!!
# overwrite original with new estimate, which is about half as much!
sumstats$n_eff <- 4 / ( 1 / counts[[ '0' ]] + 1 / counts[[ '1' ]] )

# if using ssns-ctrl, reverse signs!  (ssns-srns was fine though)
sumstats$beta <- -sumstats$beta

# load array SNP set
bim <- read_bim( name_array )
# change names to match snp_match
#names( bim ) <- c('chr', 'rsid', 'a0', 'a1')
names( bim )[ names( bim ) == 'id' ] <- 'rsid'
names( bim )[ names( bim ) == 'alt' ] <- 'a1'
names( bim )[ names( bim ) == 'ref' ] <- 'a0'
# remove unneeded column
bim[ names( bim ) == 'posg' ] <- NULL
# convert chr to integer
bim$chr <- as.integer( bim$chr )
# find matching subset
sumstats2 <- snp_match( sumstats, bim )
### ssns-srns
## 12,274,557 variants to be matched.
## 76,584 ambiguous SNPs have been removed.
## 635,789 variants have been matched; 0 were flipped and 0 were reversed.
### ssns-ctrl OLD and NEWEST
## 20,838,869 variants to be matched.
## 82,356 ambiguous SNPs have been removed.
## 672,362 variants have been matched; 0 were flipped and 0 were reversed.
nrow( bim ) # 761,366 # for reference

# before next round, correct some funny business with snp_match's output
sumstats2[ names( sumstats2 ) == 'rsid' ] <- NULL # IDs from BIM, don't want
names( sumstats2 )[ names( sumstats2 ) == 'rsid.ss' ] <- 'rsid' # restore original IDs
sumstats2[ names( sumstats2 ) == '_NUM_ID_.ss' ] <- NULL # new column added, also don't want
sumstats2[ names( sumstats2 ) == '_NUM_ID_' ] <- NULL # new column added, also don't want

# map assuming things are largely aligned already
df_beta <- snp_match( sumstats2, map )
### ssns-srns
## 635,789 variants to be matched.
## 0 ambiguous SNPs have been removed.
## 492,982 variants have been matched; 0 were flipped and 0 were reversed.
### ssns-ctrl OLD and NEWEST
## 672,362 variants to be matched.
## 0 ambiguous SNPs have been removed.
## 495,160 variants have been matched; 0 were flipped and 0 were reversed. # OLD
## 508,082 variants have been matched; 0 were flipped and 0 were reversed. # NEWEST

# here processing suggests QC on sumstats, but no code is provided (there's equations and a massive repo, may have to revisit)

# save df_beta for later use
write_tsv( df_beta, file_betas_filtered )

# corr backing file, to reload if needed
corr_file <- paste0( name, '-corr' )
# NOTE: code blindly adds another .sbk extension if I specify one (leading to double '.sbk.sbk"), so don't include it
# these ones have added extension for use internally
corr_file_rdata <- paste0( corr_file, '.RData' )

if ( file.exists( corr_file_rdata ) ) {
    # load from a different session
    load( corr_file_rdata )
} else {
    # convert basepairs to genetic position this way
    POS2 <- snp_asGeneticPos(CHR, POS, ncores = NCORES)

    # compute correlations (LD matrices) on the fly
    for (chr in 1:22) {
        message(chr)
        ## indices in 'df_beta'
        ind.chr <- which(df_beta$chr == chr)
        ## indices in 'G'
        ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
        # calculate sparse symmetric correlation matrix
        corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                         infos.pos = POS2[ind.chr2], ncores = NCORES)
        # calculate LD scores
        ld0 <- Matrix::colSums(corr0^2)
        # combine data across chromosomes
        if (chr == 1) {
            ld <- ld0
            # set up backing file, first time only
            corr <- as_SFBM(corr0, corr_file, compact = TRUE)
        } else {
            ld <- c(ld, ld0)
            corr$add_columns(corr0, nrow(corr))
        }
    }
    # save object, keep backing file, so we can load it in another session
    save( corr, ld, file = corr_file_rdata )
}

# size of LD data, in GB
file.size(corr$sbk) / 1024^3  # file size in GB
#[1] 4.118546 # ssns-srns
#[1] 4.153004 # ssns-ctrl OLD
#[1] 4.373985 # ssns-ctrl NEWEST

# decide how many individuals are in each subset... example had 503 indivs total, our data has slightly fewer, 490
# only used for LDpred2-grid and lassosum2
set.seed(1)
ind.val <- sample( nrow(G), nrow(G) * 0.7 )   # example 350 (70%) # really "training" set
ind.test <- setdiff( rows_along(G), ind.val ) # example 153 (30%) # testing set
# NEWEST:
length(ind.val)  # [1] 408
length(ind.test) # [1] 176
# save values for later loading
write_lines( ind.val, 'ind-training.txt.gz' )
write_lines( ind.test, 'ind-testing.txt.gz' )
