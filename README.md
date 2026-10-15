# Existing drugs targeting genes of hereditary pain syndromes as candidates for analgesic repurposing

This repository is the computational companion to the manuscript **"Existing drugs targeting genes of hereditary pain syndromes as candidates for analgesic repurposing: a computational knowledge discovery analysis."** It provides the R workflow and its tabular outputs; the manuscript remains the source for the study rationale, methods, interpretation, and limitations.

## Analysis

The analysis is implemented in [R/omim_drugbank.R](R/omim_drugbank.R). It links OMIM phenotype entries with a known molecular basis to human genes, then identifies approved or investigational DrugBank compounds that target those genes. The script integrates phenotype, gene, drug, and ATC classification data and labels whether pain-, analgesia-, or anesthesia-related terms appear in a drug description or ATC label. This text-based label is not evidence of an approved analgesic indication or clinical efficacy.

The analysis used an OMIM text-search export and DrugBank full database XML (version 5.1.14). These source files are not included. Obtain them from OMIM and DrugBank under their respective terms; DrugBank XML must not be redistributed. Edit the file paths and `BASE_PATH` in the configuration section of the R script for your local setup.

![Chord diagram linking OMIM phenotypes, genes, drugs, ATC therapeutic areas, and pain annotation](R/OMIM_Drug_Chord.svg)

*Figure. Connections between OMIM phenotype entries, genes, drugs, ATC anatomical main groups, and the pain/analgesia text annotation. Node colors indicate category: red, OMIM phenotype; blue, gene; green, drug; orange, ATC group; purple, pain/analgesia annotation. Connections represent database-derived associations, not evidence that a drug treats pain.*

## Running the script

Install the required R packages: `dbparser` **version 1.2.0 is required** (other versions did not work for us), `dplyr`, `stringr`, `AnnotationDbi`, `org.Hs.eg.db`, `ggplot2`, `ggsankey`, `plotly`, `htmlwidgets`, `reticulate`, `circlize`, and `viridisLite`. `org.Hs.eg.db` is distributed through Bioconductor.

The DrugBank XML is copyrighted and is not included. Obtain it from DrugBank under your own license and update `FILE_DRUGBANK_XML` in the script. On the first run, set `read_DrugBank <- TRUE` so the script parses that XML (this is slow) and stores the parsed objects in the workspace image `FILE_WORKSPACE`. On later runs, `read_DrugBank <- FALSE` loads the parsed objects from that image instead. The image is derived from the DrugBank XML and is therefore not included either. Run from the project directory:

```sh
Rscript R/omim_drugbank.R
```

The CSV tables are written to `R/`. The script also creates plots. Its Plotly SVG export uses `reticulate` and Kaleido and currently contains a machine-specific Python path; adjust that configuration for your environment if needed. The CSV tables are exported before the plotting code runs.

The script also reads `pubmed_screening_complete.csv` (see below), which is included in the repository, to create a stacked bar plot of category 1 records by drug and evidence type. By default, it includes only drugs with human evidence (`PLOT_DRUG_SET <- "human"`); setting this option to `"all"` includes all category 1 drugs and separates those with human data from those with preclinical or review evidence only. The plot is saved locally to `R/Literature_Evidence_human.pdf` and `.png` (or `Literature_Evidence_all.pdf` and `.png`); these plot files are not included in the repository.

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

## Further CSV results

The following file is not written by the script. It contains the results of the literature check, obtained by screening the PubMed records with the assistance of a large language model, with verification by the authors:

| File | Contents |
| --- | --- |
| `pubmed_screening_complete.csv` | Screening assignments for the PubMed records matching the candidate drugs. The file documents the drug assignment, screening category, evidence or exclusion code, and supporting bibliographic details. |

### `pubmed_screening_complete.csv` legend

This file reports the literature check of the 169 drugs whose DrugBank descriptions and ATC labels did not mention pain or analgesia. Each row is one PubMed record–candidate drug combination. A record naming multiple candidate drugs can therefore appear more than once: 3,662 PubMed records are represented by 3,683 rows. Counts reported in the paper are record counts, with each record counted once under its most favorable category.

| Column | Meaning |
| --- | --- |
| `PMID` | PubMed identifier. |
| `year` | Publication year. |
| `drug` | Candidate drug assigned to the record. Multiple names separated by `;` indicate a record excluded under the name-overlap rules, with all matching names listed. `unassigned` indicates one of 20 records for which no candidate drug could be identified; these records are all category 3. |
| `category` | Screening category, defined below. |
| `evidence_code` | Evidence type for categories 1 and 2, or reason for category 3. |
| `code_meaning` | Plain-text explanation of `evidence_code`. |
| `note` | Brief summary of the finding relevant to the assignment; often empty for category 3. |
| `title` | Article title. |
| `journal` | Journal abbreviation. |
| `screening_round` | `round1`: first screening batch (drug names N to Z); `round2`: second batch (remaining drugs); `auto`: excluded before screening by the name-overlap rules. |

**Categories**

- **1:** The candidate drug itself reduced pain or nociception.
- **2:** The drug was tested, but the finding was inconclusive or negative; it was given only in combination, as a derivative, or in a special formulation; or its effect was shown only at the electrophysiological or cellular level.
- **3:** Not counted as evidence of an analgesic effect of the candidate drug.

**Evidence codes (categories 1 and 2)**

- **CC:** Controlled clinical study in patients.
- **CU:** Uncontrolled clinical report (open study, case series, case report, or retrospective analysis).
- **HV:** Experimental study in healthy human volunteers.
- **PC:** Preclinical study (animal or in vitro).
- **RV:** Review, systematic review, or meta-analysis.

**Reason codes (category 3)**

- **T:** Drug used as a pharmacological tool to probe a mechanism.
- **P:** Pain was an adverse effect of the drug, or the drug was pronociceptive.
- **C:** Drug treated an underlying disease, such as angina, from which pain relief followed.
- **M:** Substance was measured, for example as a biomarker, but not administered.
- **I:** Drug was mentioned incidentally.
- **N:** No pain outcome was reported.
- **O:** Other compound, name overlap, or unspecific match. This includes 425 records excluded before screening under predefined rules: D-aspartic acid matched only via “N-methyl-D-aspartate”, formic acid only via “formalin” or “formaldehyde”, ethanol only as an extraction solvent, and copper only via “Cu/Zn superoxide dismutase”.

Screening was performed by a single reader assisted by a large language model, on title and abstract (first batch) or on title and the closing part of the abstract, which usually contains results and conclusions (second batch). The authors verified the category 1 and 2 assignments.

## PubMed search strategy

The following query was used for the PubMed search on October 10, 2026. It combines the search string used for the OMIM query with the names of the 169 drugs without mention of pain or analgesia. It returned 3,662 records. No PubMed filters were applied.

<details>
<summary>Search query</summary>

```text
 ("pain perception" OR "altered pain perception" OR "abnormal pain perception" OR "pain sensitivity" OR "abnormal pain sensitivity" OR "altered pain sensitivity" OR "pain insensitivity" OR "insensitivity to pain" OR "sensitivit* to pain" OR "insensitivit* to pain" OR "reduced sensitivity to pain" OR "increased sensitivity to pain" OR "reduced response to pain" OR "increased response to pain" OR "decreased sensitivity to pain" OR "enhanced sensitivity to pain" OR "loss of pain perception" OR "impaired pain perception" OR "diminished pain perception" OR "heightened pain perception" OR "pain threshold" OR "pain thresholds" OR "abnormal pain threshold" OR "altered pain threshold" OR "high pain threshold" OR "low pain threshold" OR hyperalgesia OR hypoalgesia OR allodynia OR nociception OR nociceptive OR "nociceptive pain" OR "nociceptive threshold" OR "congenital insensitivity to pain" OR "congenital indifference to pain" OR "loss of pain" OR "absent pain" OR "enhanced pain" OR "hereditary sensory neuropathy" OR "hereditary sensory and autonomic neuropathy" OR "primary erythromelalgia" OR "paroxysmal extreme pain disorder" OR "small fiber neuropathy") AND ("ABT-510" OR "AV-412" OR "Abivertinib" OR "Acamprosate" OR "Acetylcysteine" OR "Acoramidis" OR "Afamelanotide" OR "Afatinib" OR "Aflibercept" OR "Almonertinib" OR "Altiratinib" OR "Alvocidib" OR "Amivantamab" OR "Amlodipine" OR "Antroquinonol" OR "Aripiprazole" OR "Atagabalin" OR "Atomoxetine" OR "Azacitidine" OR "BMS-599626" OR "BMS-690514" OR "Bevacizumab" OR "Bevasiranib" OR "Borneol" OR "Bremelanotide" OR "Brigatinib" OR "Brolucizumab" OR "CUDC-101" OR "Cabozantinib" OR "Canertinib" OR "Cannabidivarin" OR "Carvacrol" OR "Carvedilol" OR "Cenegermin" OR "Cetuximab" OR "Cinnamaldehyde" OR "Copper" OR "Cycloserine" OR "D-Aspartic Acid" OR "Dacomitinib" OR "Dalteparin" OR "Decitabine" OR "Depatuxizumab mafodotin" OR "Dersimelagon" OR "Diltiazem" OR "Donepezil" OR "Eliprodil" OR "Emvododstat" OR "Entrectinib" OR "Epigallocatechin gallate" OR "Eplontersen" OR "Eptifibatide" OR "Erlotinib" OR "Eslicarbazepine" OR "Ethanol" OR "Falnidamol" OR "Faricimab" OR "Felodipine" OR "Fluciclovine (18F)" OR "Flucytosine" OR "Fluorocyclopentenylcytosine" OR "Foreskin keratinocyte (neonatal)" OR "Formic acid" OR "Fostamatinib" OR "Funapide" OR "GW 468816" OR "Gabapentin enacarbil" OR "Galunisertib" OR "Gefitinib" OR "Gingerol" OR "Gliclazide" OR "Glyburide" OR "Guadecitabine" OR "Ibutilide" OR "Icotinib" OR "Imagabalin" OR "Imidazole" OR "Inotersen" OR "Isradipine" OR "Ivonescimab" OR "Kynurenic Acid" OR "LY-3200882" OR "Lacosamide" OR "Lapatinib" OR "Larotrectinib" OR "Lazertinib" OR "Lorlatinib" OR "Magnesium carbonate" OR "Martinostat" OR "Matuzumab" OR "Memantine" OR "Merimepodib" OR "Methyl isocyanate" OR "Methylparaben" OR "Minocycline" OR "Mubritinib" OR "NM-3" OR "Naquotinib" OR "Nazartinib" OR "Necitumumab" OR "Nepidermin" OR "Neratinib" OR "Nerispirdine" OR "Nicardipine" OR "Nilvadipine" OR "Nisoldipine" OR "Olafertinib" OR "Oleandrin" OR "Olmutinib" OR "Osimertinib" OR "Oxcarbazepine" OR "PF-05089771" OR "PF-05241328" OR "Palifosfamide" OR "Panitumumab" OR "Patisiran" OR "Pegcantratinib" OR "Pelitinib" OR "Phenobarbital" OR "Phorbol 12-myristate 13-acetate diester" OR "Pidolic acid" OR "Polaprezinc" OR "Poziotinib" OR "Pralsetinib" OR "Procainamide" OR "Pyrotinib" OR "Ranibizumab" OR "Regorafenib" OR "Repotrectinib" OR "Revusiran" OR "Rilpivirine" OR "Riluzole" OR "Rindopepimut" OR "Risuteganib" OR "Rociletinib" OR "Romidepsin" OR "Rufinamide" OR "SNS-032" OR "Seletracetam" OR "Selitrectinib" OR "Sevabertinib" OR "Sorafenib" OR "Squalamine" OR "Sunvozertinib" OR "TG-100801" OR "Tafamidis" OR "Taletrectinib" OR "Tarloxotinib" OR "Tavilermide" OR "Tetraethylammonium" OR "Tiapride" OR "Trilaciclib" OR "Tuxobertinib" OR "Vactosertib" OR "Vandetanib" OR "Vanucizumab" OR "Varlitinib" OR "Vatalanib" OR "Vixotrigine" OR "Vutrisiran" OR "Xiliertinib" OR "Z-160" OR "Zalutumumab" OR "Zanubrutinib" OR "Zinc" OR "Zinc acetate" OR "Zinc chloride" OR "Zinc sulfate, unspecified form" OR "Zonisamide")
```

</details>

## Reader and reviewer note

This repository is a technical companion to the manuscript. It makes the R workflow and generated outputs inspectable, but does not replace the manuscript, which remains authoritative for the study's methods, scientific interpretation, and discussion. The analysis requires separately obtained OMIM and licensed DrugBank source data. Database-derived target associations and text-based pain labels are hypothesis-generating and do not establish analgesic efficacy.

## Citation

Use this repository as the computational companion to the paper and the implementation used to generate the reported analyses. Please cite the associated manuscript as:

Lötsch J, Himmelspach A, Kringel D. *Existing drugs targeting genes of hereditary pain syndromes as candidates for analgesic repurposing: a computational knowledge discovery analysis.* (in preparation)
