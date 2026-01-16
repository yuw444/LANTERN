library(optparse)
option_list <- list(
  make_option(
    c("-p", "--pt"),
    type = "character",
    default = "./cache/pt_matrix_chr15.tsv",
    help = "Path to PT matrix TSV file [default %default]",
    metavar = "file"
  ),
  make_option(
    c("-g", "--gene_group"),
    type = "character",
    default = "gene_group.tsv",
    help = "Path to gene group TSV file [default %default]",
    metavar = "file"
  ),
  make_option(
    c("-o", "--out_file"),
    type = "character",
    default = "output.tsv",
    help = "Output TSV path [default %default]",
    metavar = "file"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

pt_path <- opt$pt
gene_group_path <- opt$gene
out_file <- opt$out_file

library(data.table)
library(dplyr)
df_pt <- data.table::fread(
  pt_path,
  header = TRUE
) %>%
  as.data.frame()

df_pt1 <- df_pt %>%
  dplyr::select(-pos1)

df_gene_group <- setNames(
  fread(
    gene_group_path
  ),
  c("gene", "chr", "pos", "ref", "alt", "weight")
)

# df_gene_group %>%
#   group_by(gene) %>%
#   summarise(n = n()) %>%
#   pull(n) %>%
#   summary()

# str(df_gene_group)

genes_unique <- unique(df_gene_group$gene)

library(doParallel)
registerDoParallel(cores = 10)
df_weight <- foreach(genes = genes_unique, .combine = rbind.data.frame) %dopar%
  {
    vec_pos <- df_gene_group %>%
      filter(gene == genes) %>%
      pull(pos)

    pos_idx <- match(vec_pos, df_pt$pos1)

    a <- apply(df_pt1[pos_idx, ] == 03, 1, sum) ## african ancestry
    e <- apply(df_pt1[pos_idx, ] == 01, 1, sum) ## european ancestry

    return(data.frame(
      gene = genes,
      a = a,
      e = e
    ))
  }


df_weight_final <- df_weight %>%
  group_by(gene) %>%
  summarise(
    a = median(a),
    e = median(e)
  )

write.table(
  df_weight_final,
  file = out_file,
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)
