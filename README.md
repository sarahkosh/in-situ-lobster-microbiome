# Lobster in situ microbiomes: 16S rRNA analysis

Analysis code for a 16S rRNA amplicon survey of the microbial communities associated with American lobster (*Homarus americanus*) embryos and with the adult female's pleopods, carapace, and hemolymph, sampled in situ at multiple locations and embryo developmental stages.

**Status:** In preparation (unpublished dissertation chapter)

## Questions addressed

- How do community diversity and composition differ among sample types (early and late embryos, pleopods, carapace, hemolymph)?
- Do communities differ by embryo developmental stage and by sampling location?
- Which genera are enriched on pleopods compared with embryos, and on carapace compared with hemolymph?
- Does epizootic shell disease (ESD) status relate to carapace community composition?

## Repository contents

| File | Purpose |
|---|---|
| `build_phyloseq.R` | Combines the DADA2 ASV and taxonomy tables with sample metadata, removes cyanobacteria and mitochondrial sequences, and saves a phyloseq object |
| `community_analysis.R` | Alpha and beta diversity, PERMANOVA, composition barplots, shared taxa, and differential abundance (ANCOM-BC2) |
| `core.R` | Calculate core community of each sample type |
| `core_picrust.R` | Predicted function of each core community |

## Inputs

Place these in `data/`:

- `seqtab.nochim.rds`: DADA2 ASV table (samples x ASVs, chimeras removed)
- `taxa.rds`: DADA2 taxonomy table
- `map.csv`: sample metadata with the columns below

| Column | Description |
|---|---|
| `Sample` | Sample ID (matches the ASV table) |
| `Sample_Type` | Early_Embryo, Late_Embryo, Pleo_Swab, Shell, Cara_Swab, Hemolymph |
| `Sample_Kind` | Embryo, Pleopod, Carapace, Hemolymph |
| `Type` | Plot labels for sample type (e.g. "Early Embryo") |
| `Embryogenesis` | Early or Late |
| `Location` | DownEast (Maine), NH, SME [define SME] |
| `ESD_Status` | Positive or Negative |

## Usage

```r
source("01_build_phyloseq.R")        # writes data/phyob.rds
source("02_community_analysis.R")    # writes tables to results/ and figures to results/figures/
```

## Methods summary

- **Filtering:** cyanobacteria and mitochondrial sequences removed; one sample with 125 reads removed.
- **Alpha diversity:** observed richness and Shannon diversity on non-rarefied counts; Kruskal-Wallis tests, with Benjamini-Hochberg-adjusted pairwise Wilcoxon tests across sample types.
- **Beta diversity:** Bray-Curtis distances. PCoA plots use data rarefied (with replacement) to the smallest library in each subset; rarefaction depths and random seeds are recorded in the script and in `results/rarefaction_depths.csv`. PERMANOVA (`vegan::adonis2`) and dispersion tests (`betadisper`) use non-rarefied counts.
- **Shared taxa:** genera and ASVs present in at least two samples per group.
- **Differential abundance:** ANCOM-BC2 at genus level (prevalence cutoff 0.1, alpha 0.05, sensitivity analysis on). Model formulas for each comparison are specified in the script.


