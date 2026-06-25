# ============================================================================
# Differential Expression Analysis of Mistranslation Events
# ============================================================================
# Purpose: Perform differential expression analysis between tumor and normal
#          samples for mistranslation events, identify pan-cancer errors,
#          and perform functional enrichment analysis.
# ============================================================================

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)

# ============================================================================
# Data Loading and Preprocessing
# ============================================================================

# Load data containing mistranslation events and their expression levels
subs <- read.csv('subs_normexpr_filt_blacklist.csv', row.names = 1, sep = ',', header = TRUE)

# Select Tumor and Normal samples and replace 0 values with NA
data <- subs %>%
  dplyr::select(matches("Tumor|Normal")) %>%
  mutate(across(where(is.numeric), ~ na_if(., 0)))

data <- data %>%
  dplyr::select(!(\
    contains("Disqualified") | contains("KoreanReference") |\
    contains("WU.PDA1") | contains("WU.pooled")\
  ))

# ============================================================================
# Data Filtering Based on Non-NA Counts
# ============================================================================

# Retain rows where at least 20% of samples have non-NA values in each group
NA_rate_retain <- 0.2

# Separate samples by group
control_sample <- colnames(data)[grepl('Normal', colnames(data))]
case_sample <- colnames(data)[grepl('Tumor', colnames(data))]

# Calculate number of non-NA values per row for each group
calculate_non_na_count <- function(x, df) {
  return(ncol(df) - sum(is.na(x)))
}

len_control <- apply(data[, control_sample], 1, calculate_non_na_count, data[, control_sample])
len_case <- apply(data[, case_sample], 1, calculate_non_na_count, data[, case_sample])

# Identify rows meeting the retention criteria
retain_rows <- which(
  len_control >= length(control_sample) * NA_rate_retain &
  len_case >= length(case_sample) * NA_rate_retain
)
data_filt <- data[retain_rows, ]

# ============================================================================
# Differential Expression Analysis Using Wilcoxon Test
# ============================================================================
# Compare log2-transformed expression levels between tumor and normal samples
# Calculate fold change (FC) as the ratio of mean expression between groups

p_value <- c()
fold_change <- c()

for (j in 1:nrow(data_filt)) {
  # Perform Wilcoxon rank-sum test (non-parametric alternative to t-test)
  test <- wilcox.test(
    log2(unlist(data_filt[j, case_sample])),
    log2(unlist(data_filt[j, control_sample])),
    paired = FALSE
  )
  p_value <- append(p_value, test$p.value)
  
  # Calculate fold change as mean ratio
  fold_change <- append(
    fold_change,
    mean(unlist(data_filt[j, case_sample]), na.rm = TRUE) /
    mean(unlist(data_filt[j, control_sample]), na.rm = TRUE)
  )
}

# Apply FDR correction
p_adjust <- p.adjust(p_value, method = "fdr")

# Compile results
DEE <- data.frame(
  error = rownames(data_filt),
  p.value = p_value,
  p.adjust = p_adjust,
  fc = fold_change
)

DEE <- DEE %>% arrange(p.adjust)

# Export differential expression results
write.table(
  DEE,
  'differential_normexpr_error.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf("Differential expression analysis completed: %d errors analyzed\n", nrow(data_filt)))

# ============================================================================
# Pan-Cancer Events Frequency Analysis
# ============================================================================
# Identify events that are significantly differentially expressed across
# multiple cancer types, indicating pan-cancer significance

# Clear workspace for pan-cancer analysis
rm(list = ls())

# Define cancer types to analyze
cancer_types <- c(
  'BRCA', 'CCRCC', 'CRC', 'EOGC', 'GBM', 'HCC', 'HNSCC', 'LSCC',
  'LUAD', 'LUAD_cnhpp', 'OV', 'PDAC', 'UCEC'
)

# Aggregate significant events across cancer types
error_frequency <- cancer_types %>%
  map_dfr(~ {
    # Read differential expression results from each cancer type folder
    list.files(.x, pattern = "differential_normexpr_error.txt", full.names = TRUE) %>%
      read_delim(delim = "\t", col_names = TRUE, quote = "") %>%
      # Apply filtering criteria: FDR-adjusted p-value < 0.05 and fold change > 1.5
      filter(p.adjust < 0.05 & fc > 1.5) %>%
      pull(error) %>%
      tibble(error = .)
  }) %>%
  # Count occurrences of each event across cancer types
  count(error) %>%
  arrange(desc(n))

# Export pan-cancer results
write.table(
  error_frequency,
  'pancancer_differential_normexpr_error.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf("Pan-cancer analysis completed: %d significant errors identified\n", nrow(error_frequency)))

# ============================================================================
# Functional Enrichment Analysis
# ============================================================================
# Clear workspace and reload required libraries
rm(list = ls())
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)

# Select analysis type: modify to switch between "cancer_specific" and "pancancer"
analysis_type <- "cancer_specific"

# ============================================================================
# Cancer-Specific Enrichment
# ============================================================================
if (analysis_type == "cancer_specific") {
  # Load differential expression results
  DEE <- read.table(
    "differential_normexpr_error.txt",
    header = TRUE,
    fill = TRUE,
    sep = "\t",
    check.names = FALSE
  )
  
  # Classify events as Up/Down regulated or Not significant
  DEE <- DEE %>%
    mutate(
      change = case_when(
        p.adjust < 0.05 & fc > 1.5 ~ "Up",
        p.adjust < 0.05 & fc < 2^(-log2(1.5)) ~ "Down",
        TRUE ~ "Not significant"
      )
    )
  
  # Select top 30 up-regulated events
  DEE_up <- subset(DEE, change == 'Up')
  DEE_up <- DEE_up %>% arrange(desc(fc)) %>% slice_head(n = 30)
  analysis_name <- "cancer_specific"
}

# ============================================================================
# Pan-Cancer Enrichment
# ============================================================================
if (analysis_type == "pancancer") {
  # Load pan-cancer results and filter for events found in >4 cancer types
  panDE <- read.table(
    "pancancer_differential_normexpr_error.txt",
    header = TRUE,
    fill = TRUE,
    sep = "\t",
    check.names = FALSE
  )
  DEE_up <- panDE %>% filter(n > 4)
  analysis_name <- "pancancer"
}

# ============================================================================
# Gene List Preparation and Enrichment Analysis (Shared)
# ============================================================================

# Extract gene names from event identifiers (format: gene_position_substitution)
genelist <- sapply(strsplit(DEE_up$error, "_"), function(x) x[1])

# ============================================================================
# Gene Ontology (GO) Enrichment Analysis
# ============================================================================
go_results <- enrichGO(
  genelist,
  OrgDb = org.Hs.eg.db,
  ont = 'ALL',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  keyType = "SYMBOL"
)

write.table(
  go_results,
  'go_up.txt',
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
  gene = genelist,
  TERM2GENE = m_t2g_hallmark,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH"
)
Hallmark <- as.data.frame(Hallmark)

write.table(
  Hallmark,
  'hallmark_up.txt',
  sep = '\t',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ============================================================================
# Analysis Complete
# ============================================================================
cat(sprintf("\nEnrichment analysis completed: %s\n", analysis_name))
cat(sprintf("Output files: go_up.txt, hallmark_up.txt\n"))
