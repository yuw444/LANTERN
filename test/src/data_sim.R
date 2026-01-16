library(SeqArray)
library(future)
library(future.apply)
library(MASS)
library(data.table)
library(dplyr)
library(GMMAT)

gds <- seqOpen("/scratch/g/pauer/Yu/Tractor-RVA/test/data/aa_ds.gds")
gds_ids <- seqGetData(gds, "sample.id")
seqClose(gds)

set.seed(926)
ids_sub <- sample(gds_ids, 3000)
n <- length(ids_sub)

###################################
####### simulate data file #######
###################################
df_pheno <- data.frame(
    id = ids_sub,
    y = rpois(n, lambda = 5),
    age = rnorm(n, mean = 50, sd = 10),
    sex = sample(c("M", "F"), n, replace = TRUE),
    pc1 = rnorm(n),
    pc2 = rnorm(n)
)

write.table(
    df_pheno,
    file = "/scratch/g/pauer/Yu/Tractor-RVA/test/output/data_file.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

###################################
##### simulate kinship matrix #####
###################################


ids_sub2 <- sample(gds_ids, 3000)
n2 <- length(ids_sub2)
m <- 1000

# simulate SNP allele frequencies and genotypes (0,1,2)
p <- runif(m, 0.05, 0.5)
G <- matrix(rbinom(n2 * m, 2, rep(p, each = n2)), nrow = n2, ncol = m)

# center and scale columns by expected SD sqrt(2p(1-p))
Gc <- sweep(G, 2, 2 * p, "-")
sdvec <- sqrt(2 * p * (1 - p))
Z <- sweep(Gc, 2, sdvec, "/")

# genomic relationship matrix (GRM) and kinship matrix = GRM / 2
GRM <- tcrossprod(Z) / m
kinship <- GRM / 2

rownames(kinship) <- colnames(kinship) <- ids_sub2

# save results
saveRDS(kinship, "/scratch/g/pauer/Yu/Tractor-RVA/test/output/kinship.rds")