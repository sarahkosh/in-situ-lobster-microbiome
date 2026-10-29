## Build the phyloseq object from DADA2 output
##
## Input (in data/):
##   seqtab.nochim.rds  DADA2 ASV table (samples x ASVs, chimeras removed)
##   taxa.rds           DADA2 taxonomy table
##   map.csv            sample metadata; must include a "Sample" column that
##                      matches the sample names in the ASV table
## Output:
##   data/phyob.rds     filtered phyloseq object with ASV sequences attached

library(phyloseq)
library(dplyr)
library(Biostrings)

seqtab.nochim <- readRDS("data/seqtab.nochim.rds")
taxa          <- readRDS("data/taxa.rds")
map           <- read.csv("data/map.csv")
rownames(map) <- map$Sample

# Every sequenced sample needs metadata
missing_meta <- setdiff(rownames(seqtab.nochim), map$Sample)
if (length(missing_meta) > 0) {
  stop("Samples without metadata: ", paste(missing_meta, collapse = ", "))
}

phyob <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                  tax_table(taxa),
                  sample_data(map))

# Remove cyanobacteria and mitochondrial sequences (keeping taxa unassigned at those ranks)
phyob <- phyob %>%
  subset_taxa(Phylum != "Cyanobacteria" | is.na(Phylum)) %>%
  subset_taxa(Family != "Mitochondria" | is.na(Family))

# Remove a sample with only 125 reads
phyob <- subset_samples(phyob, Sample != "E_DE13_pswab")

# Attach the ASV sequences, then replace the sequence IDs with ASV1, ASV2, ...
# (done in this order so the sequences are stored in the saved object)
dna <- Biostrings::DNAStringSet(taxa_names(phyob))
names(dna) <- taxa_names(phyob)
phyob <- merge_phyloseq(phyob, dna)
taxa_names(phyob) <- paste0("ASV", seq(ntaxa(phyob)))

print(phyob)
message("Reads per sample: min = ", min(sample_sums(phyob)),
        ", median = ", median(sample_sums(phyob)),
        ", max = ", max(sample_sums(phyob)))

saveRDS(phyob, "data/phyob.rds")
