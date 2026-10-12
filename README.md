# Existing drugs targeting genes of hereditary pain syndromes as candidates for analgesic repurposing

This repository is the computational companion to the manuscript **"Existing drugs targeting genes of hereditary pain syndromes as candidates for analgesic repurposing: a computational knowledge discovery analysis."** It provides the R workflow and its tabular outputs; the manuscript remains the source for the study rationale, methods, interpretation, and limitations.

## Analysis

The analysis is implemented in [R/omim_drugbank.R](R/omim_drugbank.R). It links OMIM phenotype entries with a known molecular basis to human genes, then identifies approved or investigational DrugBank compounds that target those genes. The script integrates phenotype, gene, drug, and ATC classification data and labels whether pain-, analgesia-, or anesthesia-related terms appear in a drug description or ATC label. This text-based label is not evidence of an approved analgesic indication or clinical efficacy.

The analysis used an OMIM text-search export and DrugBank full database XML (version 5.1.14). These source files are not included. Obtain them from OMIM and DrugBank under their respective terms; DrugBank XML must not be redistributed. Edit the file paths and `BASE_PATH` in the configuration section of the R script for your local setup.

![Chord diagram linking OMIM phenotypes, genes, drugs, ATC therapeutic areas, and pain annotation](R/OMIM_Drug_Chord.svg)

*Figure. Connections between OMIM phenotype entries, genes, drugs, ATC anatomical main groups, and the pain/analgesia text annotation. Node colors indicate category: red, OMIM phenotype; blue, gene; green, drug; orange, ATC group; purple, pain/analgesia annotation. Connections represent database-derived associations, not evidence that a drug treats pain.*

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

## Further CSV results

The following files are not written by the script but created manually following PubMed search:

| File | Contents |
| --- | --- |
| `pubmed_screening_complete.csv` | Assignment of each record with drug, category, and type of evidence. |

## PubMed search strategy

The following query was used for the PubMed search underlying the manual screening on October 12, 2026. For full reproducibility, also report any PubMed filters used.

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
