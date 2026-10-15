# ==============================================================================
# Title:       OMIM Disease–DrugBank Integration and Drug Target Analysis
# Project:     FP7 NGS Rare Variants
# Description: Integrates pain-related OMIM disease phenotypes with DrugBank
#              drug and drug-target information to identify approved and
#              investigational compounds targeting genes implicated in rare
#              pain-related Mendelian diseases. Produces an integrated
#              disease–gene–drug dataset with ATC pharmacological
#              classification and pain-use annotation, and exports summary
#              tables for reporting.
#
# Authors:     Jorn Lotsch
# Date:        2026-09-05
# Revised:     2026-10-11
#
# Input:
#   - FILE_DRUGBANK_XML ("full database.xml")
#                              — DrugBank full database XML release
#                                (version 5.1.14, released 2026-01-04).
#                                NOT PROVIDED WITH THIS CODE: the DrugBank
#                                database is copyrighted and must not be
#                                redistributed. Users have to obtain the
#                                file themselves from https://go.drugbank.com
#                                under their own DrugBank license and set
#                                FILE_DRUGBANK_XML to its local path.
#   - FILE_OMIM_TSV ("OMIM-Entry-Search.tsv")
#                              — OMIM text query export (tab-separated,
#                                4-line header; obtain via https://omim.org)
#   - FILE_WORKSPACE ("DrugbankAuslesen.img")
#                              — reserved path for an optional local cache.
#                                The workspace save/load calls below are
#                                commented out, so this file is not created
#                                or loaded automatically.
#
# Switch read_DrugBank:
#   Set TRUE to read and parse FILE_DRUGBANK_XML with dbparser. DrugBank XML
#   is copyrighted and is not provided with this code; users must obtain it
#   from DrugBank under their own license and set its local path above.
#   - read_DrugBank <- TRUE   Read and parse the XML for this run. Required
#                             when DrugBank objects are not already loaded.
#   - read_DrugBank <- FALSE  Skip XML parsing; use only when the parsed
#                             DrugBank objects are already available in the
#                             current R session.
#
# Output (tables):
#   - Table_OMIM_Diseases.csv          — OMIM disease phenotype entries
#   - Table_OMIM_With_Genes.csv        — OMIM entries with Entrez gene mapping
#   - Table_Drug_Targets_Filtered.csv  — DrugBank targets overlapping OMIM genes
#   - Table_Drugs_Filtered.csv         — approved/investigational drugs for
#                                        OMIM disease gene targets
#   - Table_OMIM_Drug_Integrated.csv   — full integrated disease–gene–drug
#                                        table with ATC and pain annotation
#   - Table_OMIM_Drug_Integrated_short_unique.csv
#                                      — one row per drug with pain annotation
#
# Output (figures):
#   - OMIM_Drug_Sankey_ggplot.pdf/.png — ggplot2 Sankey diagram
#   - OMIM_Drug_Sankey_plotly.html     — interactive Plotly Sankey diagram
#   - OMIM_Drug_Sankey_plotly.svg      — Plotly Sankey (static export)
#   - OMIM_Drug_Chord.pdf/.svg         — circlize chord diagram
#
# R package dependencies:
#   dbparser (version 1.2.0 required; other versions did not work for us),
#   dplyr, stringr, AnnotationDbi, org.Hs.eg.db,
#   ggplot2, ggsankey, plotly, htmlwidgets, reticulate, circlize, viridisLite
#
# Usage:
#   1. Obtain the DrugBank XML and the OMIM export and adapt BASE_PATH and the
#      FILE_* paths in the configuration section.
#   2. First run: set read_DrugBank <- TRUE.
#   3. Set read_DrugBank <- FALSE only if parsed DrugBank objects are already
#      loaded in the current R session.
#   Rscript omim_drugbank.R
# ==============================================================================

#################################### Configuration ########################################################################

# Base paths - centralized configuration
BASE_PATH <- "/home/joern/Aktuell/FP7NGS_RareVariants"
PATH_DRUGBANK <- file.path(BASE_PATH, "11DrugBank")
PATH_R <- file.path(BASE_PATH, "08AnalyseProgramme", "R")
PATH_OMIM <- file.path(BASE_PATH, "10OMIM")

# File paths
FILE_DRUGBANK_XML <- file.path(PATH_DRUGBANK, "full database.xml")
FILE_OMIM_TSV <- file.path(PATH_OMIM, "OMIM-Entry-Search.tsv")
FILE_WORKSPACE <- file.path(PATH_R, "DrugbankAuslesen.img")

read_DrugBank <- FALSE  # To skip reading if run accidentally

#################################### Libraries ########################################################################
# Suppress package startup messages
suppressPackageStartupMessages({
  library(dbparser)
  library(dplyr)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# Only load plotting libraries if needed
# library(ggplot2)
# library(ComplexHeatmap)
# library(viridis)
# library(ggthemes)
# library(ggwordcloud)
# library(ggridges)

#################################### Helper Functions ########################################################################

#' Check if file exists and stop with informative message if not
check_file_exists <- function(filepath, description = "File") {
  if (!file.exists(filepath)) {
    stop(sprintf("%s not found at: %s", description, filepath))
  }
  invisible(TRUE)
}

#' Report join statistics
report_join_stats <- function(before, after, join_type = "join") {
  message(sprintf("%s: %d rows before, %d rows after",
                  join_type, nrow(before), nrow(after)))
}

######################### Load DrugBank Data ######################################################
if (read_DrugBank) {
  # Reference: https://cran.r-project.org/web/packages/dbparser/vignettes/dbparser.html
  
  # Validate file existence before parsing
  check_file_exists(FILE_DRUGBANK_XML, "DrugBank XML")
  
  # Parse XML database (this loads data into dbparser's internal state)
  message("Parsing DrugBank XML database...")
  read_drugbank_xml_db(FILE_DRUGBANK_XML)
  
  # Extract all relevant DrugBank tables
  message("Loading DrugBank tables...")
  drugs <- drugs()
  drug_groups <- drug_groups()
  
  # Drug targets and related tables
  drug_targets <- targets()
  message(sprintf("  - Unique drug targets: %d", length(unique(drug_targets$id))))
  
  drug_targets_actions <- targets_actions()
  drug_targets_polypeptides <- targets_polypeptides()
  message(sprintf("  - Unique target polypeptides (parent_id): %d",
                  length(unique(drug_targets_polypeptides$parent_id))))
  message(sprintf("  - Unique target polypeptides (id): %d",
                  length(unique(drug_targets_polypeptides$id))))
  
  # Additional target annotations (loaded but not actively used)
  drug_targets_polypep_ex_ident <- targets_polypep_ex_ident()
  drug_targets_polypeptides_go <- targets_polypeptides_go()
  drug_targets_polypeptides_pfams <- targets_polypeptides_pfams()
  drug_targets_polypeptides_syn <- targets_polypeptides_syn()
  
  # Enzymes and transporters (loaded but not actively used)
  drug_enzymes <- enzymes()
  drug_enzymes_actions <- enzymes_actions()
  drug_transporters <- transporters()
  drug_transporters_actions <- transporters_actions()
  
  # Optional: Save/load workspace for faster subsequent runs
  # save.image(file = FILE_WORKSPACE)
}
# load(file = FILE_WORKSPACE)

######################### Load and Process OMIM Data ######################################################

# Load OMIM pain-related entries
check_file_exists(FILE_OMIM_TSV, "OMIM TSV")
message("Loading OMIM data...")

omim_raw <- read.csv(
  FILE_OMIM_TSV,
  sep = "\t",
  skip = 4,
  stringsAsFactors = FALSE
)

message(sprintf("  - Total hits: %d", nrow(omim_raw)))

# Identify phenotype entries with known molecular basis (marked with #)
# CHANGED: corrected prefix meanings
# OMIM prefix meanings: # = phenotype with known molecular basis, * = gene,
#                       % = phenotype with unknown molecular basis,
#                       + = gene and phenotype combined
omim_diseases <- omim_raw %>%
  dplyr::mutate(
    is_disease = grepl("^#", MIM.Number),
    OMIM_clean = str_replace_all(MIM.Number, "^[#*%+]+", "")
  ) %>%
  dplyr::filter(is_disease)

message(sprintf("  - Disease entries identified: %d", nrow(omim_diseases)))

# Map OMIM IDs to ENTREZ gene IDs
message("Mapping OMIM IDs to ENTREZ gene IDs...")
omim_keys <- unique(na.omit(omim_diseases$OMIM_clean))

# Suppress expected 1:many mapping warnings
gene2entrez <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = omim_keys,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "OMIM"
  )
)

# Join gene information to OMIM data
omim_with_genes <- omim_diseases %>%
  dplyr::left_join(gene2entrez, by = c("OMIM_clean" = "OMIM"))

# CHANGED: report mapped diseases, disease-gene pairs and unique genes separately
# (the former "92/80 (115.0%)" divided disease-gene pairs by diseases)
n_pairs <- sum(!is.na(omim_with_genes$ENTREZID))
n_diseases_mapped <- omim_with_genes %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::pull(OMIM_clean) %>%
  dplyr::n_distinct()
message(sprintf("  - OMIM entries with at least one gene: %d/%d",
                n_diseases_mapped, nrow(omim_diseases)))
message(sprintf("  - Disease-gene pairs: %d; unique genes: %d",
                n_pairs, dplyr::n_distinct(omim_with_genes$ENTREZID, na.rm = TRUE)))

######################### Map Drug Targets to ENTREZ IDs ######################################################

# Add ENTREZ IDs to drug target polypeptides based on gene symbols
message("Mapping drug target genes to ENTREZ IDs...")

drug_targets_with_entrez <- drug_targets_polypeptides %>%
  dplyr::mutate(
    ENTREZ_ID = suppressMessages(
      as.character(
        mapIds(
          org.Hs.eg.db,
          keys = gene_name,
          column = "ENTREZID",
          keytype = "SYMBOL",
          multiVals = "first"  # Takes first match when multiple exist
        )
      )
    )
  )

n_drug_targets_mapped <- sum(!is.na(drug_targets_with_entrez$ENTREZ_ID))
message(sprintf("  - Drug targets with ENTREZ mapping: %d/%d (%.1f%%)",
                n_drug_targets_mapped, nrow(drug_targets_polypeptides),
                100 * n_drug_targets_mapped / nrow(drug_targets_polypeptides)))

######################### Filter to OMIM-Relevant Drugs ######################################################

# Get list of ENTREZ IDs from OMIM diseases
omim_entrez_ids <- unique(na.omit(as.character(omim_with_genes$ENTREZID)))
message(sprintf("Filtering drugs targeting OMIM disease genes (%d genes)...",
                length(omim_entrez_ids)))

# Filter drug targets to those matching OMIM genes
drug_targets_in_omim <- drug_targets_with_entrez %>%
  dplyr::filter(ENTREZ_ID %in% omim_entrez_ids)

message(sprintf("  - Drug target polypeptides matching OMIM: %d", nrow(drug_targets_in_omim)))
# CHANGED: added count of OMIM genes that are drug targets at all
message(sprintf("  - OMIM genes that are drug targets (before drug filter): %d/%d",
                dplyr::n_distinct(drug_targets_in_omim$ENTREZ_ID), length(omim_entrez_ids)))

# Get corresponding drug target records
drug_targets_filtered <- drug_targets %>%
  dplyr::filter(id %in% drug_targets_in_omim$parent_id)

message(sprintf("  - Drug targets matching OMIM: %d", nrow(drug_targets_filtered)))

# Get drug information for these targets
drugs_filtered_initial <- drugs$general_information %>%
  dplyr::filter(primary_key %in% drug_targets_filtered$parent_key)

message(sprintf("  - Initial drugs targeting OMIM genes: %d", nrow(drugs_filtered_initial)))

# Filter out drugs with unwanted groups (withdrawn, vet_approved, nutraceutical, illicit, experimental)
# Note: a drug is excluded if it carries ANY of these labels, also when it is
# additionally labelled approved or investigational.
excluded_groups <- c("withdrawn", "vet_approved", "nutraceutical", "illicit", "experimental")

# Get drug IDs that have any of the excluded groups
drugs_to_exclude <- drugs$groups %>%
  dplyr::filter(group %in% excluded_groups) %>%
  dplyr::pull(`drugbank-id`) %>%
  unique()

# CHANGED: message label (count refers to the complete DrugBank database)
message(sprintf("  - Drugs with an excluded group label in complete DrugBank: %d",
                length(drugs_to_exclude)))

# Filter out excluded drugs
drugs_filtered <- drugs_filtered_initial %>%
  dplyr::filter(!primary_key %in% drugs_to_exclude)

message(sprintf("  - Filtered drugs targeting OMIM genes: %d", nrow(drugs_filtered)))

######################### Create Integrated OMIM-Drug Dataset ######################################################

message("Building integrated OMIM-Drug dataset...")

# Join OMIM diseases with drug information through ENTREZ IDs
omim_drug_integrated <- omim_with_genes %>%
  # Link to drug target polypeptides via ENTREZ ID
  dplyr::left_join(
    drug_targets_in_omim %>%
      dplyr::select(ENTREZ_ID, parent_id),
    by = c("ENTREZID" = "ENTREZ_ID"),
    relationship = "many-to-many"
  ) %>%
  # Link to drug targets via parent_id
  # CHANGED: retained drugs only. Previously, excluded drugs (withdrawn,
  # experimental, ...) stayed in the table with drug_id and ATC code and were
  # counted in the associations, linked diseases and ATC summaries.
  dplyr::left_join(
    drug_targets_filtered %>%
      dplyr::filter(parent_key %in% drugs_filtered$primary_key) %>%
      dplyr::select(id, parent_key),
    by = c("parent_id" = "id"),
    relationship = "many-to-many"
  ) %>%
  # Link to drug information via parent_key
  dplyr::left_join(
    drugs_filtered %>%
      dplyr::select(primary_key, drug_name = name, drug_description = description),
    by = c("parent_key" = "primary_key"),
    relationship = "many-to-many"
  ) %>%
  # Add ATC codes and classifications
  dplyr::left_join(
    drugs$atc_codes %>%
      dplyr::select(drugbank_id = `drugbank-id`, atc_code, atc_level_1 = level_1,
                    atc_level_2 = level_2, atc_level_3 = level_3),
    by = c("parent_key" = "drugbank_id"),
    relationship = "many-to-many"
  ) %>%
  # Clean up intermediate columns
  dplyr::select(-parent_id) %>%
  # Rename for clarity
  dplyr::rename(drug_id = parent_key) %>%
  # CHANGED: the same polypeptide can be listed under several target records
  dplyr::distinct()

######################### Annotation: ATC System and Pain Use ######################################################

# Add ATC system classification (anatomical main group, first letter of the code)
ATC_SYSTEM_LOOKUP <- c(
  A = "Alimentary tract and metabolism",
  B = "Blood and blood forming organs",
  C = "Cardiovascular system",
  D = "Dermatologicals",
  G = "Genito-urinary system and sex hormones",
  H = "Systemic hormonal preparations, excluding sex hormones and insulins",
  J = "Antiinfectives for systemic use",
  L = "Antineoplastic and immunomodulating agents",
  M = "Musculo-skeletal system",
  N = "Nervous system",
  P = "Antiparasitic products, insecticides and repellents",
  R = "Respiratory system",
  S = "Sensory organs",
  V = "Various"
)

omim_drug_integrated <- omim_drug_integrated %>%
  dplyr::mutate(
    atc_system = dplyr::case_when(
      !is.na(atc_code) & nchar(atc_code) > 0 ~ unname(ATC_SYSTEM_LOOKUP[substr(atc_code, 1, 1)]),
      TRUE ~ "Unclassified/No ATC"
    )
  )

# CHANGED: corrected search pattern. The former pattern
# "\\b(pain|analge|anaesth|anesth)\\b" required a word boundary after the
# stem and therefore matched only the whole word "pain", never "analgesic",
# "anesthetic" or "anaesthesia". The pattern is now defined once and used
# everywhere.
pain_pattern <- "\\b(pain|analg|anaesth|anesth)"

# Add whether pain or anesthesia is mentioned for the drug
# Pain-related if the description OR an ATC subgroup label matches; decided per drug
omim_drug_integrated <- omim_drug_integrated %>%
  dplyr::mutate(
    pain_hit = (!is.na(drug_description) &
                  grepl(pain_pattern, drug_description, ignore.case = TRUE)) |
      grepl(pain_pattern, paste(atc_level_1, atc_level_2, atc_level_3), ignore.case = TRUE)
  ) %>%
  dplyr::group_by(drug_id) %>%
  dplyr::mutate(
    pain.use = dplyr::if_else(!is.na(drug_id) & any(pain_hit),
                              "Pain or analgesia mentioned",
                              "Pain or analgesia not mentioned")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-pain_hit)


######################### Final Dataset and Summary Counts ######################################################

# CHANGED: all counts are now taken from rows that carry a retained drug.
# omim_drug_integrated keeps diseases without gene or drug as rows with NA
# (left joins) and has one row per disease-gene-drug-ATC code, so its row
# count is NOT the number of associations.
omim_drug_final <- omim_drug_integrated %>%
  dplyr::filter(!is.na(drug_id), !is.na(ENTREZID))

n_assoc <- omim_drug_final %>%
  dplyr::distinct(OMIM_clean, ENTREZID, drug_id) %>%
  nrow()

message(sprintf("  - Rows in integrated table: %d (of which with a retained drug: %d)",
                nrow(omim_drug_integrated), nrow(omim_drug_final)))
message(sprintf("  - Drugs: %d (must equal filtered drugs above)",
                dplyr::n_distinct(omim_drug_final$drug_id)))
message(sprintf("  - Genes targeted by these drugs: %d/%d",
                dplyr::n_distinct(omim_drug_final$ENTREZID), length(omim_entrez_ids)))
message(sprintf("  - OMIM diseases linked to these drugs: %d/%d",
                dplyr::n_distinct(omim_drug_final$OMIM_clean),
                dplyr::n_distinct(omim_with_genes$OMIM_clean)))
message(sprintf("  - Distinct disease-gene-drug associations: %d", n_assoc))

######################### Exploratory Analysis ######################################################

# CHANGED: pain screen counted per drug, not per table row
drug_level <- omim_drug_final %>%
  dplyr::distinct(drug_id, drug_name, drug_description, pain.use)

message(sprintf("  - Drugs with pain/analgesia mentioned in description: %d/%d (%.1f%%)",
                sum(drug_level$pain.use == "Pain or analgesia mentioned"),
                nrow(drug_level),
                100 * mean(drug_level$pain.use == "Pain or analgesia mentioned")))

# Search in ATC classifications (unique subgroup labels of retained drugs)
atc_labels <- unique(na.omit(c(
  omim_drug_final$atc_level_1,
  omim_drug_final$atc_level_2,
  omim_drug_final$atc_level_3
)))
pain_in_atc <- grep(pain_pattern, atc_labels, ignore.case = TRUE, value = TRUE)
message(sprintf("  - ATC subgroup labels with pain/analgesic terms: %d", length(pain_in_atc)))
print(pain_in_atc)

# Summary of ATC codes
# CHANGED: retained drugs only; level names corrected. In DrugBank/dbparser,
# level_1 is the most specific level (chemical subgroup), level_3 the least
# specific of the three (therapeutic subgroup).
message("\nUnique ATC classifications among retained drugs:")
message(sprintf("  - ATC codes: %d",
                dplyr::n_distinct(omim_drug_final$atc_code, na.rm = TRUE)))
message(sprintf("  - Chemical subgroups (level_1): %d",
                dplyr::n_distinct(omim_drug_final$atc_level_1, na.rm = TRUE)))
message(sprintf("  - Pharmacological subgroups (level_2): %d",
                dplyr::n_distinct(omim_drug_final$atc_level_2, na.rm = TRUE)))
message(sprintf("  - Therapeutic subgroups (level_3): %d",
                dplyr::n_distinct(omim_drug_final$atc_level_3, na.rm = TRUE)))

# Summary of drug therapeutic areas
# CHANGED: drugs per anatomical main group instead of table rows. A drug with
# codes in two main groups is counted in both, so the sum can exceed the
# number of drugs. All groups are printed.
atc_system_summary <- omim_drug_final %>%
  dplyr::distinct(drug_id, atc_system) %>%
  dplyr::count(atc_system, sort = TRUE)

message("\nDrugs per ATC anatomical main group:")
print(as.data.frame(atc_system_summary))

n_drugs_multi_system <- omim_drug_final %>%
  dplyr::distinct(drug_id, atc_system) %>%
  dplyr::count(drug_id) %>%
  dplyr::filter(n > 1) %>%
  nrow()
message(sprintf("  - Drugs assigned to more than one main group: %d", n_drugs_multi_system))

message("\nTop lines of complete dataframe:")
print(names(omim_drug_integrated))
print(head(omim_drug_integrated[, setdiff(names(omim_drug_integrated), "drug_description")], 10))

omim_drug_final %>%
  dplyr::filter(pain.use == "Pain or analgesia not mentioned",
                grepl(pain_pattern, paste(atc_level_1, atc_level_2, atc_level_3),
                      ignore.case = TRUE)) %>%
  dplyr::distinct(drug_name, atc_level_2) 



######################### Export Tables for Report ######################################################

message("Exporting tables for scientific report...")

write.csv(
  omim_diseases,
  file = file.path(PATH_R, "Table_OMIM_Diseases.csv"),
  row.names = FALSE
)
message(sprintf("  - Saved: Table_OMIM_Diseases.csv (%d rows)", nrow(omim_diseases)))

write.csv(
  omim_with_genes,
  file = file.path(PATH_R, "Table_OMIM_With_Genes.csv"),
  row.names = FALSE
)
message(sprintf("  - Saved: Table_OMIM_With_Genes.csv (%d rows)", nrow(omim_with_genes)))

write.csv(
  drug_targets_filtered,
  file = file.path(PATH_R, "Table_Drug_Targets_Filtered.csv"),
  row.names = FALSE
)
message(sprintf("  - Saved: Table_Drug_Targets_Filtered.csv (%d rows)", nrow(drug_targets_filtered)))

write.csv(
  drugs_filtered,
  file = file.path(PATH_R, "Table_Drugs_Filtered.csv"),
  row.names = FALSE
)
message(sprintf("  - Saved: Table_Drugs_Filtered.csv (%d rows)", nrow(drugs_filtered)))

write.csv(
  omim_drug_integrated,
  file = file.path(PATH_R, "Table_OMIM_Drug_Integrated.csv"),
  row.names = FALSE
)
message(sprintf("  - Saved: Table_OMIM_Drug_Integrated.csv (%d rows)", nrow(omim_drug_integrated)))

# CHANGED: one row per drug taken from the drug-level table
omim_drug_integrated_short_unique <- drug_level %>%
  dplyr::select(drug_name, pain.use) %>%
  dplyr::arrange(pain.use, drug_name)

print(
  omim_drug_integrated_short_unique %>%
    dplyr::group_by(pain.use) %>%
    dplyr::tally()
)

write.csv(
  omim_drug_integrated_short_unique,
  file = file.path(PATH_R, "Table_OMIM_Drug_Integrated_short_unique.csv"),
  row.names = FALSE
)
message(sprintf("  - Saved: Table_OMIM_Drug_Integrated_short_unique.csv (%d rows)", nrow(omim_drug_integrated_short_unique)))

t2 <- omim_drug_integrated_short_unique
m  <- t2$drug_name[t2$pain.use == "Pain or analgesia mentioned"]
n  <- t2$drug_name[t2$pain.use == "Pain or analgesia not mentioned"]
m  <- m[order(tolower(m))]
n  <- n[order(tolower(n))]

n_cols <- 5
n_rows <- max(length(m), ceiling(length(n) / n_cols))   # 34
pad    <- function(x, k) c(x, rep("", k - length(x)))

table2 <- cbind(pad(m, n_rows), matrix(pad(n, n_rows * n_cols), nrow = n_rows))
colnames(table2) <- c("Pain or analgesia mentioned",
                      "Pain or analgesia not mentioned", rep("", n_cols - 1))

write.table(table2, file = file.path(PATH_R, "Table2_wide.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
              
message(sprintf("  - Saved: Table2_wide.tsv (for Word inclusion) (%d rows)", nrow(table2)))


message("Table export complete.")

######################### Plot Data ######################################################

# CHANGED: one common data set for all three plots, reduced to distinct
# disease-gene-drug-main group combinations. Previously a drug with several
# ATC codes in the same main group was drawn with a wider flow.
plot_data <- omim_drug_final %>%
  dplyr::distinct(MIM.Number, SYMBOL, drug_name, atc_system, pain.use)

message(sprintf("  - Plot data: %d rows; OMIM entries: %d; genes: %d; drugs: %d",
                nrow(plot_data),
                dplyr::n_distinct(plot_data$MIM.Number),
                dplyr::n_distinct(plot_data$SYMBOL),
                dplyr::n_distinct(plot_data$drug_name)))

######################### Create Sankey Plot ######################################################

library(dplyr)
library(ggplot2)
library(ggsankey)

# 1. Clean and prepare
sankey_clean <- plot_data

sankey_stages <- sankey_clean %>%
  dplyr::rename(
    OMIM = MIM.Number,
    Gene = SYMBOL,
    Drug = drug_name,
    ATC  = atc_system,
    Pain = pain.use
  )

# 2. Long format (this ONLY returns x, node, next_x, next_node)
df_long <- sankey_stages %>%
  make_long(OMIM, Gene, Drug, ATC, Pain)

str(df_long)
# x, node, next_x, next_node

# 3. Plot: fill by node (works with the columns you actually have)
p_omim_drug_sankey <- ggplot(
  df_long,
  aes(
    x         = x,
    next_x    = next_x,
    node      = node,
    next_node = next_node,
    fill      = factor(node)
  )
) +
  geom_sankey(flow.alpha = 0.5, node.color = "grey30", show.legend = FALSE) +
  geom_sankey_label(
    aes(label = node),
    size  = 2.5,
    color = "black",
    fill  = "white",
    label.padding = unit(0.1, "lines"),  # Minimal padding
    label.size = 0  # Remove border around labels
  ) +
  # scale_fill_viridis_d(option = "A", alpha = 0.95, name = "Node") +
  theme_sankey(base_size = 14) +
  labs(
    title = "OMIM Disease-Gene-Drug-Therapeutic Area Flow",
    x = NULL, y = NULL
  )

print(p_omim_drug_sankey)

# 4. Save
ggsave(
  filename = file.path(PATH_R, "OMIM_Drug_Sankey_ggplot.pdf"),
  plot     = p_omim_drug_sankey,
  width    = 14,
  height   = 10,
  units    = "in"
)

ggsave(
  filename = file.path(PATH_R, "OMIM_Drug_Sankey_ggplot.png"),
  plot     = p_omim_drug_sankey,
  width    = 14,
  height   = 10,
  units    = "in",
  dpi      = 300
)

######################### Alternative: Plotly Sankey ######################################################
# Load required libraries
library(plotly)
library(htmlwidgets)
library(viridisLite)
library(dplyr)

message("\nCreating Plotly Sankey diagram...")

# Prepare data for Plotly
sankey_clean <- plot_data

# OMIM ORDER - establish once upfront
unique_omim <- sort(unique(sankey_clean$MIM.Number))
n_omim <- length(unique_omim)

if (n_omim <= 12) {
  omim_colors <- viridis(max(3, n_omim), option = "C")
} else {
  omim_colors <- viridis(n_omim, option = "C")
}
names(omim_colors) <- unique_omim

# # Generate colors using a color palette
# library(RColorBrewer)
# if (n_omim <= 12) {
#   omim_colors <- brewer.pal(max(3, n_omim), "Set3")
# } else if (n_omim <= 74) {
#   omim_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_omim)
# } else {
#   omim_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_omim)
# }
# 
# library(ggthemes)
# if (n_omim <= 8) {
#   omim_colors <- ggthemes::colorblind_pal()(max(3, n_omim))
# } else if (n_omim <= 74) {
#   omim_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_omim)
# } else {
#   omim_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_omim)
# }

names(omim_colors) <- unique_omim


# Create OMIM position column FIRST
sankey_clean$omim_pos <- match(sankey_clean$MIM.Number, unique_omim)

# Generate links maintaining omim_pos order - FIXED first() error
links1 <- sankey_clean %>%
  dplyr::group_by(omim_pos, source = MIM.Number, target = SYMBOL) %>%
  dplyr::summarise(value = dplyr::n(), omim_origin = MIM.Number[1], .groups = "drop") %>%
  dplyr::arrange(omim_pos)

links2 <- sankey_clean %>%
  dplyr::group_by(omim_pos, source = SYMBOL, target = drug_name) %>%
  dplyr::summarise(value = dplyr::n(), omim_origin = MIM.Number[1], .groups = "drop") %>%
  dplyr::arrange(omim_pos)

links3 <- sankey_clean %>%
  dplyr::group_by(omim_pos, source = drug_name, target = atc_system) %>%
  dplyr::summarise(value = dplyr::n(), omim_origin = MIM.Number[1], .groups = "drop") %>%
  dplyr::arrange(omim_pos)

# links4 REVERSED by omim_pos per source (compensates Plotly flip)
links4 <- sankey_clean %>%
  dplyr::group_by(source = atc_system, omim_pos, target = pain.use) %>%
  dplyr::summarise(value = dplyr::n(), omim_origin = MIM.Number[1], .groups = "drop") %>%
  dplyr::arrange(source, desc(omim_pos), target) %>%
  dplyr::ungroup()

# Combine maintaining relative order
all_links <- dplyr::bind_rows(links1, links2, links3, links4)

# Stable node ordering
nodes <- data.frame(
  name = unique(c(as.character(all_links$source), as.character(all_links$target))), 
  stringsAsFactors = FALSE
)

all_links$source_id <- match(all_links$source, nodes$name) - 1
all_links$target_id <- match(all_links$target, nodes$name) - 1

message(sprintf("  - Nodes: %d", nrow(nodes)))
message(sprintf("  - Links: %d", nrow(all_links)))

# Colors
link_colors_rgba <- sapply(1:nrow(all_links), function(i) {
  omim_id <- all_links$omim_origin[i]
  if (!is.na(omim_id) && omim_id %in% names(omim_colors)) {
    rgb_vals <- col2rgb(omim_colors[omim_id])
    sprintf("rgba(%d,%d,%d,0.5)", rgb_vals[1], rgb_vals[2], rgb_vals[3])
  } else {
    "rgba(150,150,150,0.4)"
  }
})

node_origin <- sapply(nodes$name, function(node_name) {
  rows <- which(all_links$source == node_name | all_links$target == node_name)
  if (length(rows) == 0) return(NA_character_)
  tab <- table(all_links$omim_origin[rows])
  names(tab)[which.max(tab)]
})

node_colors <- sapply(node_origin, function(omim_id) {
  if (!is.na(omim_id) && omim_id %in% names(omim_colors)) {
    omim_colors[omim_id]
  } else {
    "#D3D3D3"
  }
})

# Plot
fig_sankey <- plot_ly(
  type = "sankey", orientation = "h", arrangement = "snap",
  node = list(
    label = nodes$name, pad = 15, thickness = 20,
    line = list(color = "black", width = 0.5),
    color = node_colors
  ),
  link = list(
    source = all_links$source_id, target = all_links$target_id,
    value = all_links$value, color = link_colors_rgba,
    line = list(color = "rgba(0,0,0,0.1)", width = 0.5)
  )
) %>%
  layout(
    title = list(text = "OMIM Disease-Gene-Drug-Therapeutic Area Flow", font = list(size = 18)),
    font = list(size = 8, color = "black"),
    xaxis = list(showgrid = FALSE, zeroline = FALSE),
    yaxis = list(showgrid = FALSE, zeroline = FALSE),
    height = 1500, width = 1500,
    paper_bgcolor = "white", plot_bgcolor = "white",
    margin = list(b = 100, t = 80, l = 50, r = 50),
    annotations = list(
      list(x = 0.05, y = -0.05, text = "OMIM", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold")),
      list(x = 0.28, y = -0.05, text = "Gene", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold")),
      list(x = 0.51, y = -0.05, text = "Drug", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold")),
      list(x = 0.74, y = -0.05, text = "ATC System", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold")),
      list(x = 0.95, y = -0.05, text = "Pain/Analgesia", showarrow = FALSE,
           xref = "paper", yref = "paper", xanchor = "center",
           font = list(size = 16, color = "black", family = "Arial, bold"))
    )
  )

print(fig_sankey)


# Save as interactive HTML
html_file <- file.path(PATH_R, "OMIM_Drug_Sankey_plotly.html")
saveWidget(fig_sankey, html_file, selfcontained = TRUE)

library(plotly)
library(reticulate)

# Set environment FIRST, before library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/joern/.Datenplatte/Subversion/.virtualenvs/r-reticulate/bin/python")

# Now load reticulate - it will use the env var
library(reticulate)

# Check - MUST show r-reticulate path
py_config()
py_module_available("plotly")
py_module_available("kaleido")

virtualenv_root()  # shows base dir
virtualenv_list()  # lists all virtualenvs

Sys.setenv(RETICULATE_PYTHON = "/home/joern/.Datenplatte/Subversion/.virtualenvs/r-reticulate/bin/python")

library(plotly)
library(reticulate)

svg_file <- file.path(PATH_R, "OMIM_Drug_Sankey_plotly.svg")

reticulate::py_install("kaleido==0.2.1", pip = TRUE)  # older stable version
plotly::save_image(fig_sankey, svg_file, width = 1500, height = 1500)

message(sprintf("  - Saved interactive HTML: %s", html_file))
message("\nTo export PNG/SVG, oopen it in a browser and export from there")

######################### Alternative: Chord Diagram ######################################################

# Load required library
library(circlize)

message("\nCreating Chord diagram...")

# Prepare data for chord diagram
# Chord diagrams work best with pairwise connections, so we'll create adjacency matrix
# We'll show connections between consecutive stages

# Create adjacency matrices for each transition
chord_data <- plot_data

chord_data$drug_name <- stringr::str_wrap(chord_data$drug_name, width = 25)

# Create a combined connection matrix
# Transition 1: OMIM -> Gene
mat1 <- chord_data %>%
  dplyr::group_by(MIM.Number, SYMBOL) %>%
  dplyr::summarise(value = dplyr::n(), .groups = "drop")

# Transition 2: Gene -> Drug
mat2 <- chord_data %>%
  dplyr::group_by(SYMBOL, drug_name) %>%
  dplyr::summarise(value = dplyr::n(), .groups = "drop")

# Transition 3: Drug -> ATC
mat3 <- chord_data %>%
  dplyr::group_by(drug_name, atc_system) %>%
  dplyr::summarise(value = dplyr::n(), .groups = "drop")

# Transition 4: ATC -> Pain
mat4 <- chord_data %>%
  dplyr::group_by(atc_system, pain.use) %>%
  dplyr::summarise(value = dplyr::n(), .groups = "drop")

# Combine all transitions into one data frame
colnames(mat1) <- c("from", "to", "value")
colnames(mat2) <- c("from", "to", "value")
colnames(mat3) <- c("from", "to", "value")
colnames(mat4) <- c("from", "to", "value")

chord_links <- dplyr::bind_rows(mat1, mat2, mat3, mat4)

message(sprintf("  - Total connections: %d", nrow(chord_links)))
message(sprintf("  - Unique nodes: %d", length(unique(c(chord_links$from, chord_links$to)))))

# Create color scheme for different node types and track categories
all_nodes <- unique(c(chord_links$from, chord_links$to))
node_colors <- rep("#999999", length(all_nodes))
names(node_colors) <- all_nodes

# Create category assignment for each node
node_categories <- rep("Other", length(all_nodes))
names(node_categories) <- all_nodes

# Color by category and assign category labels
omim_nodes <- all_nodes %in% unique(chord_data$MIM.Number)
gene_nodes <- all_nodes %in% unique(chord_data$SYMBOL)
drug_nodes <- all_nodes %in% unique(chord_data$drug_name)
atc_nodes <- all_nodes %in% unique(chord_data$atc_system)
pain_nodes <- all_nodes %in% unique(chord_data$pain.use)

node_colors[omim_nodes] <- "#E41A1C"  # Red for OMIM
node_colors[gene_nodes] <- "#377EB8"   # Blue for genes
node_colors[drug_nodes] <- "#4DAF4A"   # Green for drugs
node_colors[atc_nodes] <- "#FF7F00"    # Orange for ATC
node_colors[pain_nodes] <- "#984EA3"   # Purple for pain

node_categories[omim_nodes] <- "OMIM"
node_categories[gene_nodes] <- "Gene"
node_categories[drug_nodes] <- "Drug"
node_categories[atc_nodes] <- "ATC"
node_categories[pain_nodes] <- "Pain"

# Create a color function for links based on source node
link_colors <- node_colors[chord_links$from]

# Save as PDF
pdf_file <- file.path(PATH_R, "OMIM_Drug_Chord.pdf")
pdf(pdf_file, width = 8, height = 8)

# Create chord diagram
circos.clear()
chordDiagram(
  chord_links,
  grid.col = node_colors,
  col = link_colors,
  transparency = 0,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1),
  link.lwd = 2
)

# Add labels with category indicators
circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    category <- node_categories[sector_name]
    
    # Add category label in parentheses for first few items of each category
    label_text <- paste0(sector_name, " [", category, "]")
    
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      label_text,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.3
    )
  },
  bg.border = NA
)

# Add title
title("OMIM Disease-Gene-Drug-Therapeutic Area Connections", cex.main = 1)

# Add comprehensive legend
legend(
  "topleft",
  legend = c(
    "OMIM Disease (Red)",
    "Gene Symbol (Blue)",
    "Drug Name (Green)",
    "Therapeutic Area (Orange)",
    "Pain/Analgesia (Purple)"
  ),
  fill = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3"),
  cex = 0.5,
  bty = "n",
  title = "Node Categories"
)

dev.off()
circos.clear()
message(sprintf("  - Saved: %s", pdf_file))

# Save as SVG
svg_file <- file.path(PATH_R, "OMIM_Drug_Chord.svg")
svg(svg_file, width = 8, height = 8)

circos.clear()
chordDiagram(
  chord_links,
  grid.col = node_colors,
  col = link_colors,
  transparency = 0.2,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1),
  link.lwd = 2
)

circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    category <- node_categories[sector_name]
    
    # Add category label in parentheses for first few items of each category
    label_text <- paste0(sector_name, " [", category, "]")
    
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      label_text,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.3
    )
  },
  bg.border = NA
)

title("OMIM Disease-Gene-Drug-Therapeutic Area Connections", cex.main = 1)

legend(
  "topleft",
  legend = c(
    "OMIM Disease (Red)",
    "Gene Symbol (Blue)",
    "Drug Name (Green)",
    "Therapeutic Area (Orange)",
    "Pain/Analgesia (Purple)"
  ),
  fill = c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3"),
  cex = 0.5,
  bty = "n",
  title = "Node Categories"
)

dev.off()
circos.clear()
message(sprintf("  - Saved: %s", svg_file))

message("\nChord diagram shows all connections in circular layout.")
message("Note: Sequential flow is less clear than in Sankey/parallel sets.")

# omim_drug_final %>%
#   dplyr::filter(pain.use == "Pain or analgesia not mentioned",
#                 grepl(pain_pattern, paste(atc_level_1, atc_level_2, atc_level_3),
#                       ignore.case = TRUE)) %>%
#   dplyr::distinct(drug_name, atc_level_2)

#################################### Literature check: evidence per drug ########################################################
# Plots the category 1 records of the PubMed screening (pubmed_screening_complete.csv)
# per drug, stacked by type of evidence. The CSV is not produced by this script; it
# holds the record-level results of the literature check and is expected in PATH_R.
# PLOT_DRUG_SET: "human" = only drugs with category 1 evidence from humans
#                          (controlled clinical, uncontrolled clinical, or volunteer studies);
#                "all"   = all drugs with category 1 evidence, in two panels by tier.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridisLite)
})

FILE_SCREENING <- file.path(PATH_R, "pubmed_screening_complete.csv")
PLOT_DRUG_SET  <- "human"   # "human" or "all"
X_CAP          <- 50        # axis cap; longer bars are truncated and labelled with their total

check_file_exists(FILE_SCREENING, "PubMed screening table")
screening <- read.csv(FILE_SCREENING, stringsAsFactors = FALSE)

evidence_levels <- c(CC = "Controlled clinical",
                     CU = "Uncontrolled clinical",
                     HV = "Human volunteers",
                     PC = "Preclinical",
                     RV = "Review")
human_codes <- c("CC", "CU", "HV")

cat1 <- screening %>%
  dplyr::filter(category == 1, evidence_code %in% names(evidence_levels))

drug_tier <- cat1 %>%
  dplyr::group_by(drug) %>%
  dplyr::summarise(n_human = sum(evidence_code %in% human_codes),
                   n_total = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(tier = dplyr::if_else(n_human > 0, "Human data", "Preclinical or review only"))

if (PLOT_DRUG_SET == "human") {
  drug_tier <- dplyr::filter(drug_tier, tier == "Human data")
}

message(sprintf("Literature check plot: %d drugs (%s)", nrow(drug_tier), PLOT_DRUG_SET))

plot_lit <- cat1 %>%
  dplyr::inner_join(drug_tier, by = "drug") %>%
  dplyr::count(drug, tier, n_human, n_total, evidence_code, name = "n") %>%
  dplyr::mutate(
    evidence = factor(evidence_levels[evidence_code], levels = rev(evidence_levels)),
    drug = factor(drug, levels = drug_tier %>%
                    dplyr::arrange(n_human, n_total) %>%
                    dplyr::pull(drug))
  )

labels_trunc <- drug_tier %>%
  dplyr::filter(n_total > X_CAP) %>%
  dplyr::mutate(label = paste0("n = ", n_total, " \u2192"))

p_lit <- ggplot(plot_lit, aes(x = n, y = drug, fill = evidence)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  geom_text(data = labels_trunc, aes(x = X_CAP, y = drug, label = label),
            inherit.aes = FALSE, hjust = 1.05, size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = setNames(viridisLite::viridis(length(evidence_levels), option = "D", end = 0.95),
                                      evidence_levels),
                    breaks = evidence_levels, name = "Type of evidence") +
  coord_cartesian(xlim = c(0, X_CAP), expand = FALSE) +
  labs(x = "Records reporting an analgesic or antinociceptive effect (category 1)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

if (PLOT_DRUG_SET == "all") {
  p_lit <- p_lit + facet_grid(tier ~ ., scales = "free_y", space = "free_y")
}

print(p_lit)

plot_height <- if (PLOT_DRUG_SET == "all") 12 else 5.5
ggsave(file.path(PATH_R, paste0("Literature_Evidence_", PLOT_DRUG_SET, ".svg")),
       p_lit, width = 8, height = plot_height)
ggsave(file.path(PATH_R, paste0("Literature_Evidence_", PLOT_DRUG_SET, ".png")),
       p_lit, width = 8, height = plot_height, dpi = 300)
message(sprintf("  - Saved: Literature_Evidence_%s.pdf/.png", PLOT_DRUG_SET))

