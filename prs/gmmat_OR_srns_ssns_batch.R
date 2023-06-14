# original downloaded from here:
# /datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/gmmat_OR_srns_ssns.R

library(GMMAT)
library(genio)
library(readr)
library(data.table)

# hardcoded parameters
p_cut_max <- 1e-3
p_cut_min <- 1e-4
n_batches <- 100

# numbers for reference, on score test summary stats table before subsetting
# with p_cut_min <- 0
## nrow( snps )            # [1] 12256040
## sum( snps$PVAL < 1e-4 ) # [1] 2650
## sum( snps$PVAL < 1e-5 ) # [1] 1024
## sum( snps$PVAL < 1e-6 ) # [1] 587
## sum( snps$PVAL < 1e-7 ) # [1] 215
## sum( snps$PVAL < 1e-8 ) # [1] 41
## sum( snps$PVAL < 1e-9 ) # [1] 3

# slice between 1e-3 and 1e-4: 11321

# get batch number from terminal
args <- commandArgs( trailingOnly = TRUE )
batch <- as.numeric( args[1] )
if ( is.na( batch ) )
    stop( 'Usage: <batch> # out of ', n_batches, '!' )
if ( batch < 1 )
    stop( '`batch` must be 1 or more!' )
if ( batch > n_batches )
    stop( '`batch` must be ', n_batches, ' or less!' )

# paths
# overall shared base
base <- '/datacommons/ochoalab/ssns_gwas/'
dir_srns <- paste0( base, 'GMMAT_0418/srns_ctr/srns_ssns/' )
# for genotypes, GRM and PCs
name <- paste0( dir_srns, "srns_ssns_mac20" )
# score test output, needed for prioritizing SNPs
file_gwas_score <- paste0( dir_srns, 'glmm.score.bed.srns_ssns_PC_mac20_nohwe.txt' )
# phenotype and covariates are elsewhere
dir_pheno <- paste0( base, 'nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/' )
file_phen <- paste0( dir_pheno, 'ssns_srns_pheno.txt' ) # need extension if not .phen
file_covar <- paste0( dir_pheno, "ssns_srns_covar.txt" )
# change output path!
file_out <- paste0( base, "GMMAT_0418/PRS/glmm.wald_srns_ssns_batch-", batch, '-', n_batches, ".txt" )

# load and process phenotype data
phen <- read_phen( file_phen )
# this data codes cases/controls as 1/2, but gmmat wants 0/1, just shift down
phen$pheno <- phen$pheno - 1

# covariates don't have standard parser, and not tab-separated either!
message( 'Reading: ', file_covar )
data <- read_table( file_covar, col_names = c('fam', 'id', 'sex', 'race'), col_types = 'cccc' )

# load kinship
grm <- read_grm( name )
# load PCs
eigenvec <- read_eigenvec( name )

# make sure data is aligned
stopifnot( all( data$id == phen$id ) )
stopifnot( all( data$id == grm$fam$id ) )
stopifnot( all( data$id == eigenvec$fam$id ) )

# construct data for model
data$trait <- phen$pheno
data$PCs <- eigenvec$eigenvec

# decide which SNPs to run based on score tests (i.e. p-value thresholds)
# this file is huge, select subset to save memory
message( 'Reading: ', file_gwas_score )
snps <- fread(
    file_gwas_score,
    select = c('SNP', 'PVAL'),
    showProgress = FALSE
)

# filter data by p-value threshold, after that we only need SNPs
snps <- snps$SNP[ p_cut_min < snps$PVAL & snps$PVAL <= p_cut_max ]
# report number of total SNPs
message( 'Total SNPs (all batches): ', length( snps ) )
# get batch to calculate in this particular run
# https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
snp_batches <- chunk( snps, n_batches )
# consider current batch now:
snps <- snp_batches[[ batch ]]
message( 'Batch number: ', batch , '/', n_batches )
message( 'Total SNPs (this batch): ', length( snps ) )

message( 'Running glmm.wald!' )
x <- glmm.wald(
    trait ~ sex + race + PCs,
    data = data,
    kins = grm$kinship,
    id = "id",
    family = binomial(link = "logit"),
    infile = name,
    snps = snps
)

# save output
message( 'Writing: ', file_out )
write_tsv( x, file_out )
