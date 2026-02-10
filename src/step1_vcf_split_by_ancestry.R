library(optparse)

# write option for chr_id
option_list <- list(
    make_option(
        c("--bed"),
        type = "character",
        default = "/scratch/g/pauer/Yu/smmat/rawdata/rfmix_merged.bed",
        help = "Path to PLINK .bed file [default %default]"
    ),
    make_option(
        c("--bim"),
        type = "character",
        default = "/scratch/g/pauer/Yu/smmat/rawdata/rfmix_merged.bim",
        help = "Path to PLINK .bim file [default %default]"
    ),
    make_option(
        c("--fam"),
        type = "character",
        default = "/scratch/g/pauer/Yu/smmat/rawdata/rfmix_merged.fam",
        help = "Path to PLINK .fam file [default %default]"
    ),
    make_option(
        c("-i", "--vcf_path"),
        type = "character",
        default = "/scratch/g/pauer/Yu/smmat/src/python_split/output/chr15.maf0.01.intersected.vcf.gz",
        help = "Chromosome ID to process [default %default]"
    ),
    make_option(
        c("-o", "--out_path"),
        type = "character",
        default = "/scratch/g/pauer/Yu/smmat/src/simulation4/meta/pure_ancestry_bcf/",
        help = "Output directory [default %default]"
    ),
    make_option(
        c("-c", "--chr_id"),
        type = "character",
        default = "15",
        help = "Chromosome ID to process [default %default]"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

library(data.table)
library(dplyr)
library(tidyr)
library(doParallel)
library(vcfR)
library(stringr)
## convert vcf to gds
library(SeqArray)
library(SeqVarTools)

# Load the bed file for ancestry
library(snpStats)
df_geno <- read.plink(
    bed = opt$bed,
    bim = opt$bim,
    fam = opt$fam
)

sample_id_ancestry <- df_geno$fam$member

df_geno2 <- df_geno$genotypes@.Data

format(object.size(df_geno2), units = "Mb")

df_geno3 <- as.data.table(t(df_geno2))
regions <- colnames(df_geno2)

df_geno4 <- df_geno3 %>%
    mutate(regions = regions) %>%
    separate(regions, into = c("chr", "pos1"), sep = ":") %>%
    separate(pos1, into = c("pos1", "pos2"), sep = "-") %>%
    mutate(
        pos1 = as.numeric(pos1),
        pos2 = as.numeric(pos2)
    )

### Perform the ancestry splitting for each chromosome
# Set up parallel backend
registerDoParallel(cores = 5)
vcf_path <- opt$vcf_path
out_path <- opt$out_path

dir.exists(paste0(out_path, "/cache/")) ||
    dir.create(paste0(out_path, "/cache/"), recursive = TRUE)

df_geno4_chr <- df_geno4 %>%
    filter(chr == opt$chr_id) %>%
    dplyr::select(-chr) %>%
    mutate(across(where(is.raw), as.integer)) %>%
    as.data.table()

colnames(df_geno4_chr) <- c(
    paste(
        colnames(df_geno4_chr)[seq_len(ncol(df_geno4_chr) - 2)],
        "a",
        sep = "_"
    ),
    "pos1",
    "pos2"
)

setDT(df_geno4_chr)
setkey(df_geno4_chr, pos1, pos2)

cmd_index <- paste0(
    "bcftools index ",
    vcf_path
)

cmd <- paste0(
    'bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n" ',
    vcf_path
)

cmd_id <- paste0(
    "bcftools query -l ",
    vcf_path
)

system(cmd_index)

sample_id_vcf <- fread(cmd = cmd_id, header = FALSE)[[1]]

colnames <- c("CHROM", "POS", "ID", "REF", "ALT", sample_id_vcf)

df_vcf <- setNames(data.table(fread(cmd = cmd)), colnames)

df_vcf_gt <- df_vcf[, -c("CHROM", "POS", "ID", "REF", "ALT")]

format(object.size(df_vcf_gt), units = "Mb")

n_obs <- ncol(df_vcf_gt)
nsnps <- nrow(df_vcf_gt)

lst_gt_convert <- foreach(i = seq_len(ncol(df_vcf_gt))) %dopar%
    {
        gt <- unlist(df_vcf_gt[, ..i])
        gt_new <- sapply(strsplit(as.character(gt), "[/|]"), function(y) {
            sum(as.numeric(y), na.rm = TRUE)
        })
        return(gt_new)
    }

df_gt_convert <- do.call(cbind, lst_gt_convert)

colnames(df_gt_convert) <- sample_id_vcf

format(object.size(df_gt_convert), units = "Mb")

df_gt_convert <- as.data.table(df_gt_convert) %>%
    mutate(
        pos1 = df_vcf$POS,
        pos2 = df_vcf$POS
    )

setDT(df_gt_convert)
setkey(df_gt_convert, pos1, pos2)

df_comb_chr <- foverlaps(df_gt_convert, df_geno4_chr, nomatch = NA)

intersect_id_vcf <- intersect(sample_id_vcf, sample_id_ancestry)

df_geno5_chr <- df_comb_chr %>%
    dplyr::select(starts_with(intersect_id_vcf), pos1, pos2)

if (
    sum(
        colnames(df_geno5_chr)[seq_along(intersect_id_vcf) * 2 - 1] ==
            paste(
                colnames(df_geno5_chr)[seq_along(intersect_id_vcf) * 2],
                "_a",
                sep = ""
            )
    ) !=
        length(intersect_id_vcf)
) {
    stop("The result sample_ids does not match in two matrices")
}

df_ee_chr <- df_aa_chr <- df_gt_convert

# Table representation
# ---------------------
# Group         GT                  Frequency     Condition                 AA      EE
# AA(03)        1/1(2)              N1            (pt == 3 & gt == 2)       2       0
#               1/0(1), 0/1(1)      N2            (pt == 3 & gt == 1)       1       0
#               0/0(0)              N3            (pt == 3 & gt == 0)       0       0
#
# AE(02)        1/1(2)              N4            (pt == 2 & gt == 2)       1       1
#               1/0(1), 0/1(1)      N5            (pt == 2 & gt == 1)       p1      p2
#               0/0(0)              N6            (pt == 2 & gt == 0)       0       0
#
# EE(01)        1/1(2)              N7            (pt == 1 & gt == 2)       0       2
#               1/0(1), 0/1(1)      N8            (pt == 1 & gt == 1)       0       1
#               0/0(0)              N9            (pt == 1 & gt == 0)       0       0
#
# Formulas
# --------
# p1 = (2N1 + N2 + N4)/(2N1 + N2 + 2N4 + 2N7 + N8)
#
# p2 = (N4 + 2N7 + N8)/(2N1 + N2 + 2N4 + 2N7 + N8)
#
# The key is to find N5, sum(pt == 2 & gt == 1)
# the denominator is `D = sum(gt) - sum(pt == 2 & gt == 1)`
# First ignore N3, N6, N9, which is (gt !=0)
# p1 = sum(gt[(pt == 3 & gt >= 1)] + sum(pt == 2 & gt == 2) = D
# p2 = sum(gt[(pt == 1 & gt >= 1)] + sum(pt == 2 & gt == 2) = D

ee_f <- function(pt, gt) {
    ee <- rep(0, length(pt))
    ## singletons case
    if (sum(gt == 1) == 1 && sum(pt == 2 & gt == 1) == 1) {
        ee[gt == 1 & pt == 2] <- 0.5
        return(ee)
    }

    N5 <- sum(gt[pt == 2 & gt == 1])
    if (is.na(N5)) {
        N5 <- 0
    }
    if (N5 == sum(gt)) {
        p2 <- 0.5
    } else {
        p2 <- (sum(gt[(pt == 1 & gt >= 1)]) + sum(pt == 2 & gt == 2)) /
            (sum(gt) - N5) *
            0.9999
    }
    ee[gt == 1 & pt == 1] <- 1
    ee[gt == 1 & pt == 2] <- p2
    ee[gt == 2 & pt == 1] <- 2
    ee[gt == 2 & pt == 2] <- 1
    return(ee)
}

aa_f <- function(pt, gt) {
    aa <- rep(0, length(pt))
    ## singleton case
    if (sum(gt == 1) == 1 && sum(pt == 2 & gt == 1) == 1) {
        aa[gt == 1 & pt == 2] <- 0.5
        return(aa)
    }
    N5 <- sum(gt[pt == 2 & gt == 1])

    if (is.na(N5)) {
        N5 <- 0
    }
    if (N5 == sum(gt)) {
        p1 <- 0.5
    } else {
        p1 <- (sum(gt[(pt == 3 & gt >= 1)]) + sum(pt == 2 & gt == 2)) /
            (sum(gt) - N5) *
            0.9999
    }
    aa[gt == 1 & pt == 2] <- p1
    aa[gt == 1 & pt == 3] <- 1
    aa[gt == 2 & pt == 2] <- 1
    aa[gt == 2 & pt == 3] <- 2
    return(aa)
}

format(object.size(df_geno5_chr), units = "Mb")
format(object.size(df_gt_convert), units = "Mb")

df_pt6 <- df_geno5_chr %>% dplyr::select(ends_with("_a"))
df_gt6 <- df_geno5_chr %>% dplyr::select(!ends_with("_a"), -c(pos1, pos2))
df_pos6 <- df_gt_convert %>% dplyr::select(pos1)

df_pt6 %>%
    mutate(
        pos1 = df_pos6$pos1
    ) %>%
    write.table(
        file = paste0(out_path, "/cache/pt_matrix_chr", opt$chr_id, ".tsv"),
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
    )

stopifnot(
    all(colnames(df_gt6) == intersect_id_vcf)
)

format(object.size(df_pt6), units = "Mb")
format(object.size(df_gt6), units = "Mb")
format(object.size(df_pos6), units = "Mb")

rst_chr <- foreach(j = seq_len(nsnps)) %dopar%
    {
        pt <- unlist(df_pt6[j, ])
        gt <- unlist(df_gt6[j, ])
        pos <- unlist(df_pos6[j, ])

        if (any(is.na(pt)) || sum(gt) == 0) {
            return(FALSE)
        } else {
            return(list(
                pos = pos,
                ee = ee_f(pt, gt),
                aa = aa_f(pt, gt)
            ))
        }
    }

vcf_pos_idx_rm <- which(unlist(lapply(rst_chr, function(x) {
    identical(x, FALSE)
})))
chr_id <- opt$chr_id

pos_to_rm <- df_pos6[vcf_pos_idx_rm, ]
pos_to_rm_vec <- paste(
    "chr",
    chr_id,
    "\t",
    pos_to_rm$pos1,
    sep = ""
)

write.table(
    pos_to_rm_vec,
    file = paste0(out_path, "/cache/pos_to_rm.txt"),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
)

sample_idx_keep <- which(sample_id_vcf %in% intersect_id_vcf)
sample_id_keep <- sample_id_vcf[sample_idx_keep]
writeLines(
    sample_id_keep,
    paste0(out_path, "/cache/sample_id_keep.txt")
)

## Check if the GT sample match with PT sample
stopifnot(
    all(sample_id_keep == intersect_id_vcf)
)

## use bcftools to remove the variant that is not in the track and subset the sample
nrow_pos <- as.numeric(strsplit(
    system(paste0("wc -l ", out_path, "/cache/pos_to_rm.txt"), intern = TRUE),
    " "
)[[1]][1])
nrow_rm_ids <- as.numeric(strsplit(
    system(
        paste0("wc -l ", out_path, "/cache/sample_id_keep.txt"),
        intern = TRUE
    ),
    " "
)[[1]][1])

if (nrow_pos == 1 && nrow_rm_ids == length(sample_id_vcf)) {
    system(paste0(
        "cp ",
        vcf_path, ## subset the sample
        " ",
        out_path,
        "/cache/subset.vcf"
    ))
} else {
    system(paste0(
        "bcftools view ",
        ifelse(
            nrow_pos > 1,
            paste0("-T ^", out_path, "/cache/pos_to_rm.txt "),
            " "
        ), ## remove the variants
        ifelse(
            nrow_rm_ids != length(sample_id_vcf),
            paste0("-S ", out_path, "/cache/sample_id_keep.txt "),
            " "
        ),
        vcf_path,
        " > ", ## subset the sample
        out_path,
        "/cache/subset.vcf"
    ))
}

ee_chr <- do.call(
    rbind,
    lapply(rst_chr, function(x) {
        if (!identical(x, FALSE)) {
            return(c(chr_id, x$pos, x$ee))
        } else {
            return(NULL)
        }
    })
)
aa_chr <- do.call(
    rbind,
    lapply(rst_chr, function(x) {
        if (!identical(x, FALSE)) {
            return(c(chr_id, x$pos, x$aa))
        } else {
            return(NULL)
        }
    })
)

colnames(ee_chr) <- c("chr", "pos", intersect_id_vcf)
colnames(aa_chr) <- c("chr", "pos", intersect_id_vcf)

ee_chr <- as.data.table(ee_chr)
aa_chr <- as.data.table(aa_chr)

ee_chr_reorder <- ee_chr %>%
    dplyr::mutate(
        chr = paste("chr", chr_id, sep = ""),
        pos = format(as.numeric(pos), scientific = FALSE, nsmall = 0)
    ) %>%
    arrange(chr, pos)
aa_chr_reorder <- aa_chr %>%
    dplyr::mutate(
        chr = paste("chr", chr_id, sep = ""),
        pos = format(as.numeric(pos), scientific = FALSE, nsmall = 0)
    ) %>%
    arrange(chr, pos)

## use bcftools to annotate the DS field
## Let's use data.table to write the tsv file instead
ee_ds_file <- paste0(out_path, "/cache/ee_ds.tsv")
aa_ds_file <- paste0(out_path, "/cache/aa_ds.tsv")
data.table::fwrite(
    ee_chr_reorder,
    file = ee_ds_file,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
)

data.table::fwrite(
    aa_chr_reorder,
    file = aa_ds_file,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
)

system(
    paste0("bgzip -c ", ee_ds_file, " > ", ee_ds_file, ".gz")
)
system(
    paste0("bgzip -c ", aa_ds_file, " > ", aa_ds_file, ".gz")
)

system(
    paste0("tabix -s 1 -b 2 -e 2 ", ee_ds_file, ".gz")
)
system(
    paste0("tabix -s 1 -b 2 -e 2 ", aa_ds_file, ".gz")
)

writeLines(
    c(
        "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage of the alternate allele\">"
    ),
    paste0(out_path, "/cache/header.txt")
)

system(paste0(
    "bcftools annotate -a ",
    out_path,
    "/cache/aa_ds.tsv.gz -c CHROM,POS,FORMAT/DS ",
    "-h ",
    out_path,
    "/cache/header.txt ",
    out_path,
    "/cache/subset.vcf > ",
    out_path,
    "/aa_ds.vcf"
))

system(paste0(
    "bcftools annotate -a ",
    out_path,
    "/cache/ee_ds.tsv.gz -c CHROM,POS,FORMAT/DS ",
    "-h ",
    out_path,
    "/cache/header.txt ",
    out_path,
    "/cache/subset.vcf > ",
    out_path,
    "/ee_ds.vcf"
))

seqVCF2GDS(
    paste0(out_path, "/aa_ds.vcf"),
    paste0(out_path, "/aa_ds.gds")
)
seqVCF2GDS(
    paste0(out_path, "/ee_ds.vcf"),
    paste0(out_path, "/ee_ds.gds")
)
