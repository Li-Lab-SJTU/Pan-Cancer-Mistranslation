# ============================================================================
# Correlation Analysis of Mistranslation and Protein Expression
# ============================================================================
# Purpose: Analyze correlations between mistranslation and protein
#          expression levels, followed by functional enrichment analysis.
# ============================================================================

rm(list = ls())
library(tidyverse)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(msigdbr)
library(AnnotationDbi)

# Load mistranslation data
subs <- read.csv('subs_normexpr_filt_blacklist.csv', row.names = 1, sep = ',', header = TRUE)

# Select analysis type: modify to switch analysis ("aas" or "error")
analysis_type <- "error"

# ============================================================================
# Data Preparation - Amino Acid Substitution Analysis
# ============================================================================
# Aggregate expression by amino acid substitution type
if (analysis_type == "aas") {
  subs <- subs %>%
    select(substitution, matches("Tumor")) %>%
    mutate(across(where(is.numeric), ~ na_if(., 0))) %>%
    group_by(substitution) %>%
    summarise(
      across(everything(),
             ~ if (all(is.na(.))) NA else mean(., na.rm = TRUE)
      )
    ) %>%
    rownames_to_column(var = "rowname") %>%
    mutate(rowname = substitution) %>%
    column_to_rownames(var = "rowname") %>%
    select(-substitution)
}

# ============================================================================
# Data Preparation - Error-level Analysis
# ============================================================================
# Use error-level expression data (protein_position_substitution)
if (analysis_type == "error") {
  subs <- subs %>%
    select(matches("Tumor")) %>%
    mutate(across(where(is.numeric), ~ na_if(., 0)))
}

# ============================================================================
# Data Filtering and Preprocessing (Shared for all analysis types)
# ============================================================================

# Filter out low-quality samples
subs <- subs %>%
  select(!(
    contains("Disqualified") | contains("KoreanReference") |
    contains("WU.PDA1") | contains("WU.pooled")
  ))

# Load and process protein expression data
expression <- read.csv('protein_expression.tsv', sep = '\t', header = TRUE, row.names = 1)
expression <- expression %>% rename_with(~ str_remove_all(., "Primary."))
expression <- expression %>%
  select(matches("Tumor")) %>%
  select(!(
    contains("Disqualified") | contains("KoreanReference") |
    contains("WU.PDA1") | contains("WU.pooled") | contains("Not")
  ))

# Calculate non-NA count and filter proteins with insufficient data
NA_rate_retain <- 0.1
calculate_non_na_count <- function(x, df) {
  return(ncol(df) - sum(is.na(x)))
}
len_non_na <- apply(expression, 1, calculate_non_na_count, expression)
retain_rows <- which(len_non_na >= ncol(expression) * NA_rate_retain)
expression_filt <- expression[retain_rows, ]

# ============================================================================
# Correlation Analysis
# ============================================================================

# Align sample columns between mistranslation and protein expression data
common_samples <- intersect(colnames(subs), colnames(expression))
subs <- subs %>% select(all_of(common_samples))
expression_filt <- expression_filt %>% select(all_of(common_samples))

# Select target mistranslation for analysis
subs_filt <- subs['K to R', ]

# Initialize correlation results data frame
correlation <- data.frame(
  protein = character(),
  cor = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Calculate Pearson correlation between mistranslation and each protein
for (protein in rownames(expression_filt)) {
  cor_test_result <- try({
    cor_test <- cor.test(
      as.numeric(subs_filt[1, ]),
      as.numeric(expression_filt[protein, ]),
      method = "pearson"
    )
    correlation <- rbind(
      correlation,
      data.frame(
        protein = protein,
        cor = cor_test$estimate,
        p_value = cor_test$p.value
      )
    )
  }, silent = TRUE)
  
  if (inherits(cor_test_result, "try-error")) {
    cat("Warning: Error calculating correlation for", protein, "- skipping.\n")
  }
}

# Retrieve gene symbols from UniProt database
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

correlation <- correlation %>%
  mutate(uniprot_id = str_extract(protein, "(?<=\\|)[^|]+")) %>%
  left_join(
    getBM(
      attributes = c("uniprotswissprot", "hgnc_symbol"),
      filters = "uniprotswissprot",
      values = unique(.$uniprot_id),
      mart = mart
    ) %>%
      rename(
        uniprot_id = uniprotswissprot,
        gene_name = hgnc_symbol
      ),
    by = "uniprot_id"
  ) %>%
  select(-uniprot_id)

# Export correlation results
write.table(
  correlation,
  'correlation_KtoR_pearson.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# Functional Enrichment Analysis
# ============================================================================

# Load correlation results and filter for significant correlations
correlation <- read.table(
  "correlation_KtoR_pearson.txt",
  header = TRUE,
  fill = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Apply filtering criteria: p-value < 0.05 and correlation > 0.25
correlation_filt <- correlation %>%
  filter(p_value < 0.05 & cor > 0.25)

protein <- correlation_filt$protein

# Parse FASTA file to extract gene names and UniProt IDs
fasta_lines <- readLines("uniprotkb_reviewed_true_AND_model_organ.fasta")
gene_names <- tibble(line = fasta_lines) %>%
  filter(str_starts(line, "^>sp\\|")) %>%
  mutate(
    uniprot_id = str_extract(line, "(?<=sp\\|)[A-Z0-9]+(?=\\|)"),
    gene_name = str_extract(line, "(?<=GN=)[^ ]+")
  ) %>%
  select(uniprot_id, gene_name)

# Map protein IDs to gene names
names_df <- tibble(original_protein = protein) %>%
  mutate(uniprot_id = str_extract(original_protein, "(?<=sp\\|)[A-Z0-9]+(?=\\|)")) %>%
  left_join(gene_names, by = "uniprot_id")

# Retrieve Entrez IDs for enrichment analysis
genelist <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = names_df$gene_name,
  keytype = "SYMBOL",
  columns = "ENTREZID",
  drop = TRUE
)

# ============================================================================
# Gene Ontology (GO) Enrichment Analysis
# ============================================================================
go <- enrichGO(
  genelist$SYMBOL,
  OrgDb = org.Hs.eg.db,
  ont = 'ALL',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  keyType = "SYMBOL"
)
go <- as.data.frame(go)
write.table(
  go,
  'correlation_KtoR_pearson_go_0.25.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# KEGG Pathway Enrichment Analysis
# ============================================================================
kegg <- enrichKEGG(
  genelist$ENTREZID,
  organism = 'hsa',
  keyType = 'kegg',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  use_internal_data = FALSE
)
kegg <- as.data.frame(kegg)
write.table(
  kegg,
  'correlation_KtoR_pearson_kegg_0.25.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# Reactome Pathway Enrichment Analysis
# ============================================================================
Reactome <- enrichPathway(
  genelist$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
Reactome <- as.data.frame(Reactome)
write.table(
  Reactome,
  'correlation_KtoR_pearson_reactome_0.25.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# Hallmark Gene Set Enrichment Analysis
# ============================================================================
m_t2g_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
Hallmark <- enricher(
  gene = genelist$SYMBOL,
  TERM2GENE = m_t2g_hallmark,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH"
)
Hallmark <- as.data.frame(Hallmark)
write.table(
  Hallmark,
  'correlation_KtoR_pearson_hallmark_0.25.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# C6 Oncogenic Signature Enrichment Analysis
# ============================================================================
m_t2g_c6 <- msigdbr(species = "Homo sapiens", category = "C6") %>%
  dplyr::select(gs_name, gene_symbol)
C6 <- enricher(
  gene = genelist$SYMBOL,
  TERM2GENE = m_t2g_c6,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH"
)
C6 <- as.data.frame(C6)
write.table(
  C6,
  'correlation_KtoR_pearson_C6_0.25.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# Analysis Complete
# ============================================================================
cat("\nAnalysis completed successfully!\n")
cat("Output files generated:\n")
cat("  - correlation_KtoR_pearson.txt\n")
cat("  - correlation_KtoR_pearson_go_0.25.txt\n")
cat("  - correlation_KtoR_pearson_kegg_0.25.txt\n")
cat("  - correlation_KtoR_pearson_reactome_0.25.txt\n")
cat("  - correlation_KtoR_pearson_hallmark_0.25.txt\n")
cat("  - correlation_KtoR_pearson_C6_0.25.txt\n")
