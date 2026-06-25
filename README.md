# Mistranslation Analysis Pipeline

This repository contains a analysis pipeline for identifying mistranslation events from proteomics data and evaluating their biological relevance through differential expression, correlation, survival, and enrichment analyses.

## Project Summary

1. Generate mutation FASTA files from MSGF+ search results.
2. Remove identified spectra from MGF files.
3. Detect and filter mistranslation events.
4. Quantify differential expression of mistranslation events.
5. Correlate mistranslation events with protein expression.
6. Perform survival analysis.
7. Run functional enrichment analyses for mistranslation-associated proteins.
8. Analyze K to R mistranslation enrichment
9. Analyze cell line experiments.

## Repository Structure

```text
code/
├── identification&filter.py
├── mutatefasta_generation.py
├── unidentified_spec_mgf.py
├── differential_expression.R
├── correlation.R
├── survival.R
├── KtoR_enrichment.R
├── cellline_experiment.R
├── README.md
├── LICENSE
├── CITATION.cff
├── subs_normexpr_filt_blacklist.csv
├── snp_uniprot.txt
├── msgf_out/
└── expression/
```

## Requirements

### Python

- Python 3
- pandas
- numpy
- biopython

### R

- R 4.0 or later recommended
- tidyverse
- survival
- survminer
- clusterProfiler
- enrichplot
- org.Hs.eg.db
- msigdbr
- ReactomePA
- biomaRt
- AnnotationDbi

## Workflow

### 1. Generate mutation FASTA files

Script: `mutatefasta_generation.py`

This script extracts peptide sequences from MSGF+ output tables and generates all single amino acid substitution variants in FASTA format.

Main output:

- `mutate_peptides_<experiment>.fasta`

### 2. Remove unidentified spectra from MGF files

Script: `unidentified_spec_mgf.py`

This script reads MSGF+ results, extracts spectrum identifiers, and removes matched spectra from the corresponding MGF file.

Main output:

- `<fraction>_filtered.mgf`

### 3. Identify and filter Mistranslation events

Script: `identification&filter.py`

This script maps peptide-level events to proteins, identifies amino acid substitutions, normalizes expression values, removes known SNP-related events, and writes the final filtered table.

Main outputs:

- `subs_twostep.csv`
- `subs_normexpr_filt.csv`

### 4. Differential expression analysis

Script: `differential_expression.R`

This script performs tumor-versus-normal differential expression analysis using a Wilcoxon test, applies a pancancer frequency filter, and performs enrichment analysis.

Main outputs:

- `differential_normexpr_error.txt`
- `pancancer_differential_normexpr_error.txt`
- `go_up.txt`
- `hallmark_up.txt`

### 5. Correlation analysis

Script: `correlation.R`

This script computes Pearson correlations between a selected mistranslation event and protein expression, then performs GO, KEGG, Reactome, Hallmark, and C6 enrichment analyses on the correlated proteins.

Main outputs:

- `correlation_KtoR_pearson.txt`
- `correlation_KtoR_pearson_go_0.25.txt`
- `correlation_KtoR_pearson_kegg_0.25.txt`
- `correlation_KtoR_pearson_reactome_0.25.txt`
- `correlation_KtoR_pearson_hallmark_0.25.txt`
- `correlation_KtoR_pearson_C6_0.25.txt`

### 6. Survival analysis

Script: `survival.R`

This script converts event-level expression into survival-ready input, identifies optimal cutpoints, and performs Cox regression and log-rank testing.

Main outputs:

- `survival_normexpr_error_cox.txt`
- `survival_normexpr_error.txt`

### 7. K to R enrichment analysis

Script: `KtoR_enrichment.R`

This script summarizes proteins associated with the `K to R` substitution and performs GO, KEGG, Hallmark, and C6 enrichment analyses.

Main outputs:

- `KtoR_proteincount.txt`
- `enrichprotein_KtoR_KEGG.txt`
- `enrichprotein_KtoR_GO.txt`
- `enrichprotein_KtoR_hallmark.txt`
- `enrichprotein_KtoR_C6.txt`

### 8. Cell line experiment analysis

Script: `cellline_experiment.R`

This script analyzes cell lines data for the K to R substitution, performs differential expression analysis, and runs GSEA-style enrichment analyses.

Main outputs:

- `go_KtoR.txt`
- `Hallmark_KtoR.txt`
- `DEP_KtoR.txt`
- `GSEA_GO_KtoR.txt`
- `GSEA_KEGG_KtoR.txt`
- `GSEA_Hallmark_KtoR.txt`

## Running the Scripts

Example commands for the Python scripts:

```bash
python identification&filter.py
python mutatefasta_generation.py
python unidentified_spec_mgf.py
```

Example commands for the R scripts:

```r
source('differential_expression.R')
source('correlation.R')
source('survival.R')
source('KtoR_enrichment.R')
source('cellline_experiment.R')
```
## License

This project is distributed under the MIT License. See [LICENSE](LICENSE) for details.
