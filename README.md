# Existing drugs targeting genes of hereditary pain syndromes as candidates for analgesic repurposing

This repository accompanies the scientific report **"Existing drugs targeting genes of hereditary pain syndromes as candidates for analgesic repurposing: a computational knowledge discovery analysis."** It provides the R workflow and its tabular outputs; the report remains the primary source for the study rationale, methods, interpretation, and limitations.

**Suggested citation:** Jörn Lötsch, André Himmelspach, and Dario Kringel. *Existing drugs targeting genes of hereditary pain syndromes as candidates for analgesic repurposing: a computational knowledge discovery analysis* (in preparation).

![Chord diagram linking OMIM phenotypes, genes, drugs, ATC therapeutic areas, and pain annotation](R/OMIM_Drug_Chord.svg)

*Figure. Connections between OMIM phenotype entries, genes, drugs, ATC anatomical main groups, and the pain/analgesia text annotation. Node colors indicate category: red, OMIM phenotype; blue, gene; green, drug; orange, ATC group; purple, pain/analgesia annotation. Connections represent database-derived associations, not evidence that a drug treats pain.*

## Analysis

The analysis is implemented in [`R/omim_drugbank.R`](R/omim_drugbank.R). It links OMIM phenotype entries with a known molecular basis to human genes, then identifies approved or investigational DrugBank compounds that target those genes. The script integrates phenotype, gene, drug, and ATC classification data and labels whether pain-, analgesia-, or anesthesia-related terms appear in a drug description or ATC label. This text-based label is not evidence of an approved analgesic indication or clinical efficacy.

The analysis used an OMIM text-search export and DrugBank full database XML (version 5.1.14). These source files are not included. Obtain them from OMIM and DrugBank under their respective terms; DrugBank XML must not be redistributed. Edit the file paths and `BASE_PATH` in the configuration section of the R script for your local setup.

## Running the script

Install the required R packages: `dbparser` **version 1.2.0 is required** (other versions did not work for us), `dplyr`, `stringr`, `AnnotationDbi`, `org.Hs.eg.db`, `ggplot2`, `ggsankey`, `plotly`, `htmlwidgets`, `reticulate`, `circlize`, and `viridisLite`. `org.Hs.eg.db` is distributed through Bioconductor.

The DrugBank XML is copyrighted and is not included. Obtain it from DrugBank under your own license and update `FILE_DRUGBANK_XML` in the script. For a fresh R session, set `read_DrugBank <- TRUE` so the script reads and parses that XML, then run from the project directory:

```sh
Rscript R/omim_drugbank.R
```

The CSV tables are written to `R/`. The script also creates plots. Its Plotly SVG export uses `reticulate` and Kaleido and currently contains a machine-specific Python path; adjust that configuration for your environment if needed. The CSV tables are exported before the plotting code runs.

## CSV outputs

The following files are written by the script to `R/`:

| File | Contents |
| --- | --- |
| `Table_OMIM_Diseases.csv` | OMIM phenotype entries selected as having a known molecular basis. |
| `Table_OMIM_With_Genes.csv` | Selected OMIM entries with their Entrez Gene IDs and gene symbols, where mapped. |
| `Table_Drug_Targets_Filtered.csv` | DrugBank target records whose genes overlap the OMIM-derived gene set. |
| `Table_Drugs_Filtered.csv` | Retained approved or investigational human drugs targeting those genes, after excluding specified DrugBank groups. |
| `Table_OMIM_Drug_Integrated.csv` | Integrated phenotype-gene-drug records with drug description, ATC annotations, and the text-based pain/analgesia label. |
| `Table_OMIM_Drug_Integrated_short_unique.csv` | One row per drug with its name and pain/analgesia annotation. |


## Reader and reviewer note

This repository is a technical companion to the manuscript. It makes the R workflow and generated outputs inspectable, but does not replace the manuscript, which remains authoritative for the study's methods, scientific interpretation, and discussion. The analysis requires separately obtained OMIM and licensed DrugBank source data. Database-derived target associations and text-based pain labels are hypothesis-generating and do not establish analgesic efficacy.

## Citation

Use this repository as the computational companion to the paper and the implementation used to generate the reported analyses. Please cite the associated manuscript as:

Lötsch J, Himmelspach A, Kringel D. *Existing drugs targeting genes of hereditary pain syndromes as candidates for analgesic repurposing: a computational knowledge discovery analysis.* (in preparation)
