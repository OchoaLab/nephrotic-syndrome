library(genpwr)

# total sample size, which we might want to increase
N <- 250
# prevalence of SRNS
p <- 0.2
# actual numbers of cases and controls
n_case <- N * p
n_ctrl <- N * ( 1 - p )
# expected risk haplotype proportions in the population
af_case <- 0.29
af_ctrl <- 0.12
# derive OR from here
or <- af_case * ( 1 - af_ctrl ) / ( af_ctrl * ( 1 - af_case ) )
# and MAF in sample
af <- ( af_case * n_case + af_ctrl * n_ctrl ) / N

pw <- genpwr.calc(
    calc = "power",
    model = "logistic",
    N = N,
    Case.Rate = p,
    MAF = af, # seq(0.05, 0.3, by = 0.05),
    OR = or, # seq(1.1, 3.0, by = 0.1),
    Alpha = 0.01, # 5e-8,
    True.Model = "Additive", 
    Test.Model = "Additive"
)
