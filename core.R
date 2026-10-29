### in situ core microbiome ###
#from community_analysis.R

phyob <- readRDS("phyob.rds")

library(phyloseq)
library(microbiome)
library(dplyr)
library(stringr)
library(VennDiagram)
library(grid)

##agglomerate phyloseq object to genera

phyob_genus <- tax_glom(phyob, taxrank = "Genus", NArm = TRUE)

#phyloseq objects of each sample type
phyob_embryo <- subset_samples(phyob_genus, Sample_Kind == "Embryo")
phyob_embryo <- prune_taxa(taxa_sums(phyob_embryo) > 0, phyob_embryo) #prune genera with 0 
prev.90.embryo <- core_members(phyob_embryo, detection = 0, prevalence = 0.89)
emb_core_ps <- prune_taxa(prev.90.embryo, phyob_embryo)
emb_taxa_core <- as.data.frame(tax_table(emb_core_ps))


phyob_pleopod <- subset_samples(phyob_genus, Sample_Kind == "Pleopod")
phyob_pleopod <- prune_taxa(taxa_sums(phyob_pleopod) > 0, phyob_pleopod)
prev.90.pleopod <- core_members(phyob_pleopod, detection = 0, prevalence = 0.89)

phyob_carapace <- subset_samples(phyob_genus, Sample_Kind == "Carapace")
phyob_carapace <- prune_taxa(taxa_sums(phyob_carapace) > 0, phyob_carapace)
prev.98.carapace <- core_members(phyob_carapace, detection = 0, prevalence = .97) # no genera present on 99 or 100% of samples. 
prev.90.carapace <- core_members(phyob_carapace, detection = 0, prevalence = 0.89)

phyob_hemolymph <- subset_samples(phyob_genus, Sample_Kind == "Hemolymph")
phyob_hemolymph <- prune_taxa(taxa_sums(phyob_hemolymph) > 0, phyob_hemolymph)
prev.89.hemolymph <- core_members(phyob_hemolymph, detection = 0, prevalence = .49)
prev.70.hemolymph <- core_members(phyob_hemolymph, detection = 0, prevalence = .69)

## Check intersections for Venn diagram 
n1 <- length(prev.90.embryo)
n2 <- length(prev.90.pleopod)
n3 <- length(prev.90.carapace)
n4 <- length(prev.70.hemolymph)

n12 <- length(intersect(prev.90.embryo, prev.90.pleopod))
n13 <- length(intersect(prev.90.embryo, prev.90.carapace))
n14 <- length(intersect(prev.90.embryo, prev.70.hemolymph))
n23 <- length(intersect(prev.90.pleopod, prev.90.carapace))
n24 <- length(intersect(prev.90.pleopod, prev.70.hemolymph))
n34 <- length(intersect(prev.90.carapace, prev.70.hemolymph))
list123 <- list(prev.90.embryo, prev.90.pleopod, prev.90.carapace) #combine list of 3
asv123 <- Reduce(intersect, list123) #list of shared taxa
n123 <- length(Reduce(intersect, list123)) #count of shared taxa for venn 

list124 <- list(prev.90.embryo, prev.90.pleopod, prev.70.hemolymph)
n124 <- length(Reduce(intersect, list124))

list134 <- list(prev.90.embryo, prev.90.carapace, prev.70.hemolymph)
n134 <- length(Reduce(intersect, list134))

list234 <- list(prev.90.pleopod, prev.90.carapace, prev.70.hemolymph)
n234 <- length(Reduce(intersect, list234))

list1234 <- list(prev.90.embryo, prev.90.pleopod, prev.90.carapace, prev.70.hemolymph)
n1234 <- length(Reduce(intersect, list1234))



venn.1 <- draw.quad.venn(area1 = n1,
                         area2 = n2,
                         area3 = n3,
                         area4 = n4,
                         n12 = n12,
                         n13 = n13,
                         n14 = n14,
                         n23 = n23,
                         n24 = n24,
                         n34 = n34,
                         n123 = n123,
                         n124 = n124,
                         n134 = n134,
                         n234 = n234,
                         n1234 = n1234,
                         category = c("Embryo", "Pleopod", "Carapace", "Hemolymph"),
                         fill = c("orange", "red", "green", "blue"))









## shared genera - no prevalence cut off

taxa.embryo <- taxa_names(phyob_embryo)
taxa.pleopod <- taxa_names(phyob_pleopod)
taxa.carapace <- taxa_names(phyob_carapace)
taxa.hemolymph <- taxa_names(phyob_hemolymph)

## Check intersections for Venn diagram 
m1 <- length(taxa.embryo)
m2 <- length(taxa.pleopod)
m3 <- length(taxa.carapace)
m4 <- length(taxa.hemolymph)

m12 <- length(intersect(taxa.embryo, taxa.pleopod))
m13 <- length(intersect(taxa.embryo, taxa.carapace))
m14 <- length(intersect(taxa.embryo, taxa.hemolymph))
m23 <- length(intersect(taxa.pleopod, taxa.carapace))
m24 <- length(intersect(taxa.pleopod, taxa.hemolymph))
m34 <- length(intersect(taxa.carapace, taxa.hemolymph))
mlist123 <- list(taxa.embryo, taxa.pleopod, taxa.carapace) #combine list of 3
masv123 <- Reduce(intersect, mlist123) #list of shared taxa
m123 <- length(Reduce(intersect, mlist123)) #count of shared taxa for venn 

mlist124 <- list(taxa.embryo, taxa.pleopod, taxa.hemolymph)
m124 <- length(Reduce(intersect, mlist124))

mlist134 <- list(taxa.embryo, taxa.carapace, taxa.hemolymph)
m134 <- length(Reduce(intersect, mlist134))

mlist234 <- list(taxa.pleopod, taxa.carapace, taxa.hemolymph)
m234 <- length(Reduce(intersect, mlist234))

mlist1234 <- list(taxa.embryo, taxa.pleopod, taxa.carapace, taxa.hemolymph)
m1234 <- length(Reduce(intersect, mlist1234))



venn.2 <- draw.quad.venn(area1 = m1,
                         area2 = m2,
                         area3 = m3,
                         area4 = m4,
                         n12 = m12,
                         n13 = m13,
                         n14 = m14,
                         n23 = m23,
                         n24 = m24,
                         n34 = m34,
                         n123 = m123,
                         n124 = m124,
                         n134 = m134,
                         n234 = m234,
                         n1234 = m1234,
                         category = c("Embryo", "Pleopod", "Carapace", "Hemolymph"),
                         fill = c("orange", "red", "green", "blue"))



## identify 207 shared taxa from m1234
shared.genera.all <-Reduce(intersect, mlist1234)
phyob_shared_all <- prune_taxa(shared.genera.all, phyob_genus) #make phyloseq object of only shared genera

#save taxa table as csv
taxa_df_all <- as.data.frame(tax_table(phyob_shared_all))
write.csv(taxa_df_all, "shared_all_taxa_table.csv", row.names = FALSE)




