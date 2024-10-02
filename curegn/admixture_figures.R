# CureGN admixture results

library(tidyverse)
library(popkin)
library(genio)
library(ochoalabtools)

# admixture plot for manuscript
# covar file contains all individuals before excluding KO0939 and KO1059
# need all samples to match with ADMIXTURE output Q

setwd( '/datacommons/ochoalab/curegn/merge_tgp/admixture' )

name_Q <- "curegn_tgp_merge_maf10.5.Q"

# load CureGN info
data_curegn <- read_tsv("curegn_covar_diseasesubtype.txt", col_types = 'cicclll') %>%
    mutate( dataset = 'CureGN', ssns_srns = ssns | srns ) %>%
    select( -sex, -race, -ssns, -srns )

# load TGP too
data_tgp <- read_tsv( '/datacommons/ochoalab/ssns_gwas/imputed/patient-data.txt.gz', col_types = 'cccccdciiii' ) %>% 
    filter( dataset == "tgp" ) %>% 
    select( id, ancestry ) %>% 
    mutate( ancestry = toupper( ancestry ), dataset = 'TGP', ns = FALSE, ssns_srns = FALSE )

# merge CureGN with TGP
data <- bind_rows( data_curegn, data_tgp )

# read admixture Q matrix:
Q <- read_matrix( name_Q )
# Q and data are aligned, use that to label ancestries systematically
colnames(Q) <- admix_label_cols( Q, data$ancestry )

# these are already aligned, merge Q and data!
Q <- bind_cols(Q, data)
# reorder columns with distance from Africa, for plot to look nicer
Q <- Q %>% relocate( AFR, EUR, SAS, EAS, AMR ) %>% 
    # remove individuals that no longer consent
    filter(!id %in% c("KO0939", "KO1059"))

# Figure (A) CureGN + TGP
# split groups by ancestry and dataset (curegn/tgp)
# reuse previously-calculated major ancestry filters of 80%, except SAS by 50%
new_label_order <- c("AFR\nCureGN", "AFR_admix\nCureGN", "AFR\nTGP", "AFR_admix\nTGP",
                     "EUR\nCureGN", "EUR_admix\nCureGN", "EUR\nTGP",
                     "SAS\nCureGN", "Asian_admix\nCureGN", "SAS\nTGP", "Asian_admix\nTGP",
                     "EAS\nCureGN", "EAS\nTGP", 
                     "Other\nCureGN", "AMR\nTGP")
# arrange individuals by each ancestry subgroup separately
Q = Q %>%
    mutate( new_label = factor( paste0(ancestry,"\n",dataset), levels = new_label_order ) ) %>% 
    arrange(
        new_label,
        ifelse( ancestry %in% c('AFR', 'AFR_admix'), -AFR,
        ifelse( ancestry %in% c('EUR', 'EUR_admix'), -EUR,
        ifelse( ancestry %in% c('SAS', 'Asian_admix'), -SAS,
        ifelse( ancestry == 'EAS', -EAS, -AMR )
        )))
    )

# just for plot legend, rename ancestry proportion vars
Q = Q %>% rename(African = AFR, European = EUR, 'South Asian' = SAS, 'East Asian' = EAS, 'Native American' = AMR) 

# Figure (B)
# cleaned up admixture plot after 80% filter for AFR and EUR and 50% filter for SAS
# remove TGP
# filter for NS only
Q2 = Q %>% 
    filter(dataset == "CureGN") %>%  
    filter( id %in% data$id[ data$ns ] )

# Figure (C)
# like B but filter further for SSNS & SRNS only
Q3 = Q2 %>% 
    filter(id %in% data$id[ data$ssns_srns ] )

wh <- fig_scale( 2/1.5 )
panel_letters_adj <- -0.065
fig_start( name_Q, width = wh[1], height = wh[2], mar_t = 2 )
par( oma = c(2, 0, 0, 0) )
# unfortunately need to define a layout
layout(
    matrix(
        c(1:3, 2, 4, 2),
        nrow = 3,
        byrow = TRUE
    ),
    widths = c(1, 0.1)
)
plot_admix(
    Q[,1:5],
    labs = as.character( Q$new_label ),
    labs_cex = 0.6,
    labs_line = 3,
    labs_even = TRUE,
    labs_even_line = 2,
    xlab = NA,
    panel_letters = 'A',
    panel_letters_adj = panel_letters_adj,
    layout_add = FALSE,
    leg_width = 0.1,
    leg_cex = 0.8,
    leg_mar = c(8, 0, 8, 3)
)
plot_admix(
    Q2[,1:5],
    labs = Q2$ancestry,
    labs_cex = 0.7,
    labs_line = 2,
    labs_even = TRUE,
    xlab = NA,
    panel_letters = 'B',
    panel_letters_adj = panel_letters_adj,
    leg_omit = TRUE
)
plot_admix(
    Q3[,1:5],
    labs = Q3$ancestry,
    labs_cex = 0.7,
    labs_line = 2,
    labs_even = TRUE,
    xlab_line = 3.7,
    panel_letters = 'C',
    panel_letters_adj = panel_letters_adj,
    leg_omit = TRUE
)
fig_end()

# report sample sizes for legend
message( 'A) n = ', nrow( Q ) )
table( Q$dataset )
message( 'B) n = ', nrow( Q2 ) )
table( Q2$ancestry )
message( 'C) n = ', nrow( Q3 ) )
table( Q3$ancestry )
