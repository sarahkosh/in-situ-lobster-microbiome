### Core Picrust ###
# Run picrust on all  core communities to see how functions differ#

### Make fasta file and ASV count table for Picrust analysis ###
library(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)

#load in core ps objects
setwd("/Location/to/ps_objects")
emb.core <- readRDS("emb_core90.rds")
hem.core <- readRDS("hem_core70.rds")
pleo.core <- readRDS("pleo_core90.rds")
cara.core <- readRDS("cara_core90.rds")
#merge ps objects
core.all <- merge_phyloseq(emb.core, hem.core)
core.all <- merge_phyloseq(core.all, pleo.core)
core.all <- merge_phyloseq(core.all, cara.core)

# Extract OTU/ASV table from desired phyloseq object(s)
asv.tab.emb <- as(otu_table(emb.core), "matrix")
asv.tab.hem <- as(otu_table(hem.core), "matrix")
asv.tab.pleo <- as(otu_table(pleo.core), "matrix")
asv.tab.cara <- as(otu_table(cara.core), "matrix")

asv.tab.core <- as(otu_table(core.all), "matrix")

# Ensure taxa are rows
if (!taxa_are_rows(emb.core)) {asv.tab.emb <- t(asv.tab.emb)}

if (!taxa_are_rows(hem.core)) {asv.tab.hem <- t(asv.tab.hem)}

if (!taxa_are_rows(pleo.core)) {asv.tab.pleo <- t(asv.tab.pleo)}

if (!taxa_are_rows(cara.core)) {asv.tab.cara <- t(asv.tab.cara)}

if (!taxa_are_rows(core.all)) {asv.tab.core <- t(asv.tab.core)}

# Write table (TAB-delimited, no quotes)
write.table(asv.tab.emb,
  file = "/core.picrust/asv_abund_emb_core.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA)

write.table(asv.tab.hem,
            file = "/core.picrust/asv_abund_hem_core.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

write.table(asv.tab.pleo,
            file = "/core.picrust/asv_abund_pleo_core.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

write.table(asv.tab.cara,
            file = "core.picrust/asv_abund_cara_core.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

write.table(asv.tab.core,
            file = "core.picrust/asv_abund_core.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = NA)


##picrust requires exact ASVs

# Extract sequences
library(Biostrings)

# Extract sequences
dna.emb <- refseq(emb.core)
dna.hem <- refseq(hem.core)
dna.pleo <- refseq(pleo.core)
dna.cara <- refseq(cara.core)
dna.core <- refseq(core.all)

# Ensure names match ASV IDs in the table
names(dna.emb) <- taxa_names(emb.core)
names(dna.hem) <- taxa_names(hem.core)
names(dna.pleo) <- taxa_names(pleo.core)
names(dna.cara) <- taxa_names(cara.core)

names(dna.core) <- taxa_names(core.all)

# Write FASTA
writeXStringSet(dna.emb, "/core.picrust/emb_sequences.fasta")
writeXStringSet(dna.hem, "/core.picrust/hem_sequences.fasta")
writeXStringSet(dna.pleo, "/core.picrust/pleo_sequences.fasta")
writeXStringSet(dna.cara, "/core.picrust/cara_sequences.fasta")
writeXStringSet(dna.core, "core.picrust/core_sequences.fasta")

# Ensure names match ASV IDs in the table
identical(rownames(asv.tab.emb), names(dna.emb))
identical(rownames(asv.tab.hem), names(dna.hem))
identical(rownames(asv.tab.pleo), names(dna.pleo))
identical(rownames(asv.tab.cara), names(dna.cara))
identical(rownames(asv.tab.core), names(dna.core))

# Check Dimensions
dim(asv.tab.cara)   # ASVs × samples
length(dna.cara)    # number of ASVs


### Picrust2 pipeline: This takes about 15 minutes to run

## in terminal ##
# conda activate picrust2
# picrust2_pipeline.py -s emb_sequences.fasta -i asv_abund_emb_core.tsv -o picrust2_out_emb -p 1
### Continue with ANCOMBC2 and ggpicrust2 for visualizations and analysis ####

