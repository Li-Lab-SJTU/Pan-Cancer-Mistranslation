# ============================================================================
# Cell Line (HEK293T) Mistranslation Experiment Analysis
# ============================================================================
# Purpose: Analyze K-to-R mistranslation in HEK293T cell line experiments,
#          including enrichment of K-to-R substitution proteins, differential
#          expression between K-R and mock conditions, and gene set
#          enrichment analysis (GSEA) of the differential expression results.
# ============================================================================

rm(list = ls())

library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(ReactomePA)
library(data.table)
library(biomaRt)

# ============================================================================
# Enrichment Analysis of K-to-R Substitution Proteins
# ============================================================================
# Perform GO and Hallmark gene set enrichment on the proteins bearing the
# K-to-R substitution identified in the HEK293T experiment.
# ============================================================================

df <- read.csv(
  'subs_normexpr_filt_HEK293T_KtoR.csv',
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

genelist <- df$protein

# ---------------------------------------------------------------------------
# Gene Ontology (GO) Enrichment
# ---------------------------------------------------------------------------
go <- enrichGO(
  genelist,
  OrgDb = org.Hs.eg.db,
  ont = 'ALL',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  keyType = 'SYMBOL'
)
go <- as.data.frame(go)
write.table(
  go,
  'go_KtoR.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------------------------------------------------------
# Hallmark Gene Set Enrichment
# ---------------------------------------------------------------------------
m_t2g <- msigdbr(species = 'Homo sapiens', category = 'H') %>%
  dplyr::select(gs_name, gene_symbol)

Hallmark <- enricher(
  gene = genelist,
  TERM2GENE = m_t2g,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = 'BH'
)
Hallmark <- as.data.frame(Hallmark)
write.table(
  Hallmark,
  'Hallmark_KtoR.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# K-R vs Mock Differential Expression Analysis
# ============================================================================
# Identify differentially expressed proteins between K-R and mock conditions
# in HEK293T cells using a Wilcoxon rank-sum test on median-normalized,
# log2-transformed expression values.
# ============================================================================

rm(list = ls())

expr <- read.table(
  'all_expression_protein_HEK293T.tsv',
  sep = '\t',
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Replace zero values with NA, then median-normalize and log2-transform
expr[expr == 0] <- NA
col_medians <- apply(expr, 2, median, na.rm = TRUE)
expr_norm <- sweep(expr, 2, col_medians, '/')
expr_log2 <- log2(expr_norm)

# Keep only proteins with complete (non-NA) observations across samples
expr_filt <- expr_log2[complete.cases(expr_log2), , drop = FALSE]
expr_norm_filt <- expr_norm[rownames(expr_filt), , drop = FALSE]

# Locate K-R (experimental) and mock (control) sample columns
exp_cols <- grep('K-R', colnames(expr_filt))
ctrl_cols <- grep('mock', colnames(expr_filt))

if (length(exp_cols) == 0 || length(ctrl_cols) == 0) {
  stop('There are no K-R or mock samples in the expression data. Please check the column names.')
}

# Wilcoxon rank-sum test for each protein
p_value <- apply(expr_filt, 1, function(x) {
  wilcox.test(
    as.numeric(x[exp_cols]),
    as.numeric(x[ctrl_cols]),
    paired = FALSE
  )$p.value
})

# Fold change as the ratio of mean expression between groups
foldchange <- apply(expr_norm_filt, 1, function(x) {
  mean(as.numeric(x[exp_cols]), na.rm = TRUE) /
    mean(as.numeric(x[ctrl_cols]), na.rm = TRUE)
})

# Compile differential expression results with FDR correction
DEP <- data.frame(
  protein = rownames(expr_filt),
  foldchange = foldchange,
  p.value = p_value,
  p.adjust = p.adjust(p_value, method = 'fdr'),
  row.names = NULL,
  check.names = FALSE
)

# Parse FASTA headers to map UniProt IDs to gene names
fasta_lines <- readLines('uniprotkb_reviewed_true_AND_model_organ.fasta')
gene_names <- tibble(line = fasta_lines) %>%
  filter(str_starts(line, '^>sp\\|')) %>%
  mutate(
    uniprot_id = str_extract(line, '(?<=sp\\|)[A-Z0-9]+(?=\\|)'),
    gene_name = str_extract(line, '(?<=GN=)[^ ]+')
  ) %>%
  dplyr::select(uniprot_id, gene_name)

DEP <- DEP %>%
  mutate(uniprot_id = str_extract(protein, '(?<=sp\\|)[A-Z0-9]+(?=\\|)')) %>%
  left_join(gene_names, by = 'uniprot_id')

write.table(
  DEP,
  'DEP_KtoR.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# Gene Set Enrichment Analysis (GSEA) of Differential Expression Results
# ============================================================================
# Perform GSEA using GO Biological Process, KEGG, and Hallmark gene sets on
# the ranked list of differentially expressed proteins.
# ============================================================================

# Build the ranked gene list from fold change
genelist_input <- DEP %>%
  mutate(log2fc = log2(foldchange)) %>%
  dplyr::select(gene_name, log2fc)

# Map gene symbols to Entrez IDs
genename <- as.character(genelist_input[, 1])
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genename,
  keytype = 'SYMBOL',
  columns = 'ENTREZID'
)

# Remove duplicated symbols
non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)
gene_map <- gene_map[non_duplicates_idx, ]
colnames(gene_map)[1] <- 'gene_name'

# Join with fold change values and build the named ranked list
temp <- inner_join(gene_map, genelist_input, by = 'gene_name')
temp <- temp[, -1]
temp <- na.omit(temp)
temp$log2fc <- sort(temp$log2fc, decreasing = TRUE)
geneList <- temp[, 2]
names(geneList) <- as.character(temp[, 1])

# ---------------------------------------------------------------------------
# GSEA: Gene Ontology (Biological Process)
# ---------------------------------------------------------------------------
Go_gseresult <- gseGO(
  geneList,
  'org.Hs.eg.db',
  keyType = 'ENTREZID',
  ont = 'BP',
  pvalueCutoff = 0.05
)
go_results <- as.data.frame(Go_gseresult)
write.table(
  go_results,
  'GSEA_GO_KtoR.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------------------------------------------------------
# GSEA: KEGG Pathways
# ---------------------------------------------------------------------------
KEGG_gseresult <- gseKEGG(geneList, pvalueCutoff = 0.05)
kegg_results <- as.data.frame(KEGG_gseresult)
write.table(
  kegg_results,
  'GSEA_KEGG_KtoR.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------------------------------------------------------
# GSEA: Hallmark Gene Sets
# ---------------------------------------------------------------------------
msig_h_df <- msigdbr(species = 'Homo sapiens', category = 'H')
msig_h <- msig_h_df[, c('gs_name', 'gene_symbol')]

genelist_input$log2fc <- sort(genelist_input$log2fc, decreasing = TRUE)
geneList <- genelist_input[, 2]
names(geneList) <- as.character(genelist_input[, 1])

Hallmark_gseresult <- GSEA(
  geneList = geneList,
  TERM2GENE = msig_h,
  minGSSize = 10,
  pvalueCutoff = 0.05
)
Hallmark_results <- as.data.frame(Hallmark_gseresult)
write.table(
  Hallmark_results,
  'GSEA_Hallmark_KtoR.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

cat('Cell line (HEK293T) mistranslation analysis completed.\n')
