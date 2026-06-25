genelist<-AnnotationDbi::select(org.Hs.eg.db, keys=data,keytype="SYMBOL", columns = "ENTREZID", drop = T)
# ============================================================================
# K to R Enrichment Analysis
# ============================================================================
# Purpose: Identify proteins associated with the K to R substitution and
#          perform GO, KEGG, and gene set enrichment analyses.
# ============================================================================

rm(list = ls())

library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
input_file <- 'subs_normexpr_filt_blacklist.csv'
target_substitution <- 'K to R'
minimum_protein_count <- 2

output_protein_count <- 'KtoR_proteincount.txt'
output_kegg <- 'enrichprotein_KtoR_KEGG.txt'
output_go <- 'enrichprotein_KtoR_GO.txt'
output_hallmark <- 'enrichprotein_KtoR_hallmark.txt'
output_c6 <- 'enrichprotein_KtoR_C6.txt'

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
write_enrichment_table <- function(result, file_path) {
  write.table(
    as.data.frame(result),
    file_path,
    sep = '\t',
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
}

run_msigdb_enrichment <- function(genes, category, output_file) {
  term_to_gene <- msigdbr(species = 'Homo sapiens', category = category) %>%
    dplyr::select(gs_name, gene_symbol)

  enrichment <- enricher(
    gene = genes,
    TERM2GENE = term_to_gene,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    pAdjustMethod = 'BH'
  )

  write_enrichment_table(enrichment, output_file)
}

# ---------------------------------------------------------------------------
# Protein frequency summary for the target substitution
# ---------------------------------------------------------------------------
subs <- read.csv(input_file, row.names = 1, sep = ',', header = TRUE)
target_data <- subs %>% filter(substitution == target_substitution)

protein_count <- sort(table(target_data$protein), decreasing = TRUE)
write.table(
  as.data.frame(protein_count),
  output_protein_count,
  sep = '\t',
  col.names = c('protein', 'count'),
  row.names = FALSE,
  quote = FALSE
)

proteins_for_enrichment <- names(protein_count)[protein_count >= minimum_protein_count]

# ---------------------------------------------------------------------------
# Map proteins to gene symbols and Entrez IDs
# ---------------------------------------------------------------------------
gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = proteins_for_enrichment,
  keytype = 'SYMBOL',
  columns = 'ENTREZID',
  drop = TRUE
)

# Keep only rows with valid Entrez IDs for downstream enrichment
gene_mapping <- gene_mapping %>% filter(!is.na(ENTREZID))

# ---------------------------------------------------------------------------
# Enrichment analyses
# ---------------------------------------------------------------------------
kegg <- enrichKEGG(
  gene_mapping$ENTREZID,
  organism = 'hsa',
  keyType = 'kegg',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  use_internal_data = FALSE
)
write_enrichment_table(kegg, output_kegg)

go <- enrichGO(
  gene_mapping$SYMBOL,
  OrgDb = org.Hs.eg.db,
  ont = 'ALL',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  keyType = 'SYMBOL'
)
write_enrichment_table(go, output_go)

run_msigdb_enrichment(gene_mapping$SYMBOL, 'H', output_hallmark)
run_msigdb_enrichment(gene_mapping$SYMBOL, 'C6', output_c6)

cat('K to R enrichment analysis completed.\n')
