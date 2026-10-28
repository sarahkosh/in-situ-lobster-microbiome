setwd("/Users/jskoshak/Desktop/in_situ")

taxa <- readRDS("taxa.rds") # load in your taxa object 
seqtab.nochim <- readRDS("seqtab.nochim.rds") # load in your asv object


phyob <- readRDS("phyob.rds")



physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   tax_table(taxa))
physeq

map <- read.csv("map.csv")
map <- sample_data(map)
rownames(map) <- map$Sample 
print(map)

phyloseq_merge <- merge_phyloseq(physeq, map)
print(phyloseq_merge) 

####Subsetting####

phyob <- phyloseq_merge %>%  # remove cyanobacteria and mitochondrial DNA
  subset_taxa(Phylum !="Cyanobacteria"| is.na(Phylum)) %>%
  subset_taxa(Family != "Mitochondria" | is.na(Family))

phyob <- subset_samples(phyob, Sample != "E_DE13_pswab") #this sample only had 125 reads so I removed it
phyob_early <-subset_samples(phyob, Embryogenesis =="Early")
phyob_late <- subset_samples(phyob, Embryogenesis =="Late")

#store sequences in refseq of phyloseq object but change to ASV1, ASV2 etc.
dna <- Biostrings::DNAStringSet(taxa_names(phyob))
names(dna) <- taxa_names(phyob)
phyob1 <- merge_phyloseq(phyob, dna)
taxa_names(phyob) <- paste0("ASV", seq(ntaxa(phyob)))


### update sample_data with new map file##
map <- read_csv("map.csv")
map_df <- as.data.frame(map)
rownames(map_df) <- map_df$Sample 
sample_data(phyob) <- sample_data(map_df)
phyob
saveRDS(phyob, file = "phyob.rds")

library(phyloseq) 
library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(microbiome)

#### Alpha Diversity #####
sample_order <- c(
  "Early_Embryo",
  "Hemolymph",
  "Late_Embryo",
  "Pleo_Swab",
  "Shell",
  "Cara_Swab")

phyob@sam_data$Sample_Type <- factor(
  phyob@sam_data$Sample_Type,
  levels = sample_order)


richness_plot <- plot_richness(phyob, x="Sample_Type", color="Sample_Type", measures=c("Observed", "Shannon"),  title = "Alpha Diversity (non-rarefied)", nrow = 3, scales = "free_y") +
  theme(axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 17),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 17, angle = 35, hjust = 0.9, vjust = 0.88),
        strip.text = element_text(size = 17),
        axis.ticks.length = unit(.3 ,"cm"),
        plot.title.position = "panel",
        plot.title = element_text(size = 19)) +
  geom_point(aes(color= Sample_Type), alpha=1, size=7, shape = 20) +
  scale_shape_manual(values = 15, name = "") +
  scale_color_manual(values = c("Early_Embryo" = "#779fe5", "Late_Embryo" = "#3f78a6", "Hemolymph" = "#519c3d", "Pleo_Swab" = "#7272cc", "Shell" = "#c17fee", "Cara_Swab" = "#345fa9", name = "")) +
  theme(legend.position = "none",
        panel.grid = element_line(color = "grey87", linetype = 4),
        panel.background = element_rect(fill = "grey95"))+
  scale_x_discrete(labels = c(Early_Embryo = 'Early Embryo',Late_Embryo = 'Late Embryo', Hemolymph = 'Hemolymph', Cara_Swab = 'Carapace Swab', Shell = 'Shell', Pleo_Swab = 'Pleopod Swab'))+
  stat_summary(fun.data=mean_se, geom ="errorbar", colour="gray10", size=0.9, width = 0.25) 

##box and whisker plot instead 
library(ggpubr)

ggboxplot(rich.stats, x = "Sample_Type", y = "Shannon", 
          color = "Sample_Type", palette = c("#779fe5", "#3f78a6", "#519c3d", "#7272cc", "#c17fee", "#345fa9"),
          order = c("Early_Embryo", "Hemolymph", "Late_Embryo", "Pleo_Swab", "Shell", "Cara_Swab"),
          ylab = "Alpha Diversity", xlab = "Sample Type")

### Alpha Diversity Statistics ###
library(FSA)
rich.stats <- estimate_richness(phyob, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(phyob)), by = "Sample")

kruskal.test(Shannon ~ Sample_Type, data = rich.stats) #Alpha diversity significantly different between sample types, but don't know which groups are different

pairwise.wilcox.test(rich.stats$Shannon, rich.stats$Sample_Type, p.adjust.method = "BH") #pairwise comparison - All are significantly different from each other except: Pleo_swab&Shell, Hemolymph&Early, Late&Hemolymph, Shell&Cara_Swab, Shell&Pleo_Swab


## average shannon of each sample type
e_emb_rich <- rich.stats %>%
  select(Shannon, Observed, Sample_Type) %>%
  filter(Sample_Type == "Early_Embryo") 
e_emb_shan <- mean(e_emb_rich$Shannon) #2.0
e_emb_obs <- mean(e_emb_rich$Observed) #191.4

l_emb_rich <- rich.stats %>%
  select(Shannon, Observed, Sample_Type) %>%
  filter(Sample_Type == "Late_Embryo") 
l_emb_shan <- mean(l_emb_rich$Shannon) #3.1
l_emb_obs <- mean(l_emb_rich$Observed) #367.6

pleo_rich <- rich.stats %>%
  select(Observed, Shannon, Sample_Kind) %>%
  filter(Sample_Kind == "Pleopod") 
pleo_shan <- mean(pleo_rich$Shannon) #3.8
pleo_obs <- mean(pleo_rich$Observed) #615.3

cara_rich <- rich.stats %>%
  select(Observed, Shannon, Sample_Kind) %>%
  filter(Sample_Kind == "Carapace") 
cara_shan <- mean(cara_rich$Shannon) #4.3
cara_obs <- mean(cara_rich$Observed) #695.1

hem_rich <- rich.stats %>%
  select(Observed, Shannon, Sample_Kind) %>%
  filter(Sample_Kind == "Hemolymph") 
hem_shan <- mean(hem_rich$Shannon) #2.8
hem_obs <- mean(hem_rich$Observed) #219.8




### PCoA ###
min_lib <- min(sample_sums(phyob)) #find minimum number of reads in sample with lowest number of reads

min_lib #3738
#subsample
set.seed(3)
subhab_all <- rarefy_even_depth(phyob, sample.size = min_lib, verbose=FALSE, replace=TRUE) #randomly subsample 
subhab_all #new ps object with rarefied and subsampled ASVs 

hab_pcoa_all <- ordinate(physeq = subhab_all, method ='PCoA', distance = 'bray')

all_PCOA <- plot_ordination(physeq = subhab_all, ordination = hab_pcoa_all, color='Type', shape = 'Type') +
  geom_point(aes(color = Type), alpha=1, size=4) +
  stat_ellipse(aes(group = Type), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals1, name = "") +
  scale_shape_manual(values = vals2, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

vals1 <- c("Early Embryo" = "#E69F00", "Late Embryo"  = "#56B4E9","Hemolymph"   = "#009E73", "Carapace"    = "#0072B2","Pleopod"     = "#CC79A7")
vals2 <- c('Early Embryo' = 16, 'Late Embryo' = 16, 'Carapace' = 17, 'Pleopod' = 15, 'Hemolymph' = 18)

## PERMANOVA & Betadisper : All samples ~ Sample_Type##
library(vegan)
all_bray <- phyloseq::distance(phyob, method = "bray") 
df_all <- data.frame(sample_data(phyob))
adon_all <- adonis2(all_bray ~ Sample_Type, data = df_all) # p = 0.001

disp_all <- betadisper(all_bray, df_all$Sample_Type)
anova_all_disp <- anova(disp_all) #p= 1.002 e-09. dispersions are different
all_tukey <- TukeyHSD(disp_all)
plot(disp_all)


#### Taxonomy Barplot - All Samples ####
#### Phylum barplot #####
phyob_phylum <- phyob %>%
  tax_glom(taxrank='Phylum') %>% #select/combine at specified taxa level
  transform_sample_counts(function(x) {(x/sum(x))*100} ) %>% #convert to relative abundance
  psmelt() %>% #switch to ggplot usable format
  filter(Abundance > 1) %>% #filter by taxa with relative abundance >1%
  arrange(Phylum) #sort alphabetically 

phylum_bp <- ggplot(phyob_phylum, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  facet_grid(~factor(Type, levels=c("Early Embryo", "Late Embryo", "Pleopod", "Carapace", "Hemolymph")), scales = "free", space = "free")+ #makes panels on your plot
  guides(x = guide_axis(angle = 90)) + #sample names turn 90
  geom_bar(stat='identity', color = "black") + # black around bars
  scale_fill_manual(values = set_colors40) + #colors set below
  guides(fill = guide_legend(reverse = FALSE, keywidth = 2, keyheight = 1, ncol = 1)) + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 11))+
  ylab('Relative Phylum Abundance (>1%) \n') +
  xlab('Sample ID')

phylum_bp

### To make it interactive ###
library(plotly)
ggplotly(phylum_bp)




########################## EMBRYOS ###############################
##################################################################


### Early Embryos ###
ps_early_emb <- subset_samples(phyob, Sample_Type == "Early_Embryo")
ps_early_emb <- prune_taxa(taxa_sums(ps_early_emb)>0, ps_early_emb)

##PCOA##


ml_early <- min(sample_sums(ps_early_emb))
set.seed(4)
sub_earl <- rarefy_even_depth(ps_early_emb, sample.size = ml_early, verbose=FALSE, replace=TRUE)
sub_earl 

pcoa_earl <- ordinate(physeq = sub_earl, method ='PCoA', distance = 'bray')

earl_PCOA <- plot_ordination(physeq = sub_earl, ordination = pcoa_earl, color='Location', shape = 'Location') +
  geom_point(aes(color = Location), alpha=1, size=4) +
  stat_ellipse(aes(group = Location), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals3, name = "") +
  scale_shape_manual(values = vals4, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

vals3 <- c("DownEast" = "#779fe5", "NH" = "#519c3d")
vals4 <- c('DownEast' = 16, 'NH' = 17)

## PERMANOVA & Betadisper##
early_bray <- phyloseq::distance(ps_early_emb, method = "bray") 
df_earl_emb <- data.frame(sample_data(ps_early_emb))
adon_earl_emb <- adonis2(early_bray ~ Location, data = df_earl_emb) # p=0.264. no difference from location in early embryos

disp_earl_emb_loc <- betadisper(early_bray, df_earl_emb$Location)
anova_early_emb_loc <- anova(disp_earl_emb_loc) #p=0.1946. dispersions are not different

##Alpha div by location early emb ##
rich.earl <- estimate_richness(ps_early_emb, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(ps_early_emb)), by = "Sample")

kruskal.test(Shannon ~ Location, data = rich.earl) #chi-squared = 0.077143, df = 1, pval = 0.3798




### Late Embryos ###
ps_late_emb <- subset_samples(phyob, Sample_Type == "Late_Embryo")
ps_late_emb <- prune_taxa(taxa_sums(ps_late_emb)>0, ps_late_emb)

ml_late <- min(sample_sums(ps_late_emb))
set.seed(5)
sub_late <- rarefy_even_depth(ps_late_emb, sample.size = ml_late, verbose=FALSE, replace=TRUE)
sub_late 

pcoa_late <- ordinate(physeq = sub_late, method ='PCoA', distance = 'bray')

late_PCOA <- plot_ordination(physeq = sub_late, ordination = pcoa_late, color='Location', shape = 'Location') +
  geom_point(aes(color = Location), alpha=1, size=4) +
  stat_ellipse(aes(group = Location), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals3, name = "") +
  scale_shape_manual(values = vals4, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 


## PERMANOVA & Betadisper##
ps_emb_nosme <- subset_samples(ps_late_emb, Location != "SME")
ps_emb_nosme <- prune_taxa(taxa_sums(ps_emb_nosme)>0, ps_emb_nosme)
late_bray <- phyloseq::distance(ps_emb_nosme, method = "bray") 
df_late_emb <- data.frame(sample_data(ps_emb_nosme))
adon_late_emb <- adonis2(late_bray ~ Location, data = df_late_emb) # p=0.205

disp_late_emb_loc <- betadisper(late_bray, df_late_emb$Location)
anova_late_emb_loc <- anova(disp_late_emb_loc) #p=0.1532. dispersions are not different

##Alpha div by location late emb ##
rich.late <- estimate_richness(ps_emb_nosme, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(ps_emb_nosme)), by = "Sample")

kruskal.test(Shannon ~ Location, data = rich.late) #chi-squared = 0.077143, df = 1, pval = 0.3798



## All embryos

ps_emb <- subset_samples(phyob, Sample_Kind == "Embryo")
ps_emb <- prune_taxa(taxa_sums(ps_emb)>0, ps_emb)

ml_emb <- min(sample_sums(ps_emb))
set.seed(6)
sub_emb <- rarefy_even_depth(ps_emb, sample.size = ml_emb, verbose=FALSE, replace=TRUE)
sub_emb

pcoa_emb <- ordinate(physeq = sub_emb, method ='PCoA', distance = 'bray')

emb_PCOA <- plot_ordination(physeq = sub_emb, ordination = pcoa_emb, color='Sample_Type', shape = 'Sample_Type') +
  geom_point(aes(color = Sample_Type), alpha=1, size=4) +
  stat_ellipse(aes(group = Sample_Type), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals5, name = "") +
  scale_shape_manual(values = vals6, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.25, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

vals5 <- c("Late_Embryo" = "#779fe5", "Early_Embryo" = "#519c3d")
vals6 <- c('Late_Embryo' = 16, 'Early_Embryo' = 17)

## PERMANOVA & Betadisper##
emb_bray <- phyloseq::distance(ps_emb, method = "bray") 
df_emb <- data.frame(sample_data(ps_emb))
adon_emb <- adonis2(emb_bray ~ Embryogenesis, data = df_emb) 

disp_emb_loc <- betadisper(emb_bray, df_emb$Embryogenesis)
anova_emb_loc <- anova(disp_emb_loc) 
plot(disp_emb_loc)

## Alpha Diversity between early and late
rich.emb <- estimate_richness(ps_emb, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(ps_emb)), by = "Sample")

kruskal.test(Shannon ~ Sample_Type, data = rich.emb) #pval 4.662e-0.08


### Shared genera of early and late ###
# Remove taxa without genus assignment
emb_ps_genus <- subset_taxa(ps_emb, !is.na(Genus) & Genus != "") #remove asvs where genera is na, 2287 asvs left
emb_ps_genus <- tax_glom(emb_ps_genus, taxrank = "Genus", NArm = TRUE) #agglomerate to genus, 588 genera

emb_earl <- subset_samples(emb_ps_genus, Embryogenesis == "Early")
emb_late <- subset_samples(emb_ps_genus, Embryogenesis == "Late")

emb_earl <- prune_taxa(taxa_sums(emb_earl) > 0, emb_earl) #452 genera
emb_late <- prune_taxa(taxa_sums(emb_late) > 0, emb_late) #364 genera

genus_emb_earl <- taxa_names(emb_earl)
genus_emb_late <- taxa_names(emb_late)

asv_emb_earl <- otu_table(emb_earl)
asv_emb_late <- otu_table(emb_late)

genus_emb_earl <- taxa_names(emb_earl)[rowSums(asv_emb_earl > 0) >= 2] #present in at least 2 samples, 452
genus_emb_late <- taxa_names(emb_late)[rowSums(asv_emb_late > 0) >= 2] #364


#### Make Venn Diagram
library(VennDiagram)
library(grid)

draw.pairwise.venn(
  area1 = length(genus_emb_earl),
  area2 = length(genus_emb_late),
  cross.area = length(intersect(genus_emb_earl, genus_emb_late)),
  category = c("Early Embryos", "Late Embryos"),
  fill = c("#519c3d", "tomato"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.2
)

## make a list of shared genera
genus_earl_list <- as.character(tax_table(emb_earl)[, "Genus"])
genus_late_list <- as.character(tax_table(emb_late)[, "Genus"])
intersect(genus_earl_list, genus_late_list)



### Shared ASVs of Early and Late Embryos
asv_early_embryo <- taxa_names(ps_early_emb)
asv_late_embryo <- taxa_names(ps_late_emb)

otu_A <- otu_table(ps_early_emb)
otu_B <- otu_table(ps_late_emb)

asv_A <- taxa_names(asv_early_embryo)[rowSums(otu_A > 0) >= 2]
asv_B <- taxa_names(asv_late_embryo)[rowSums(otu_B > 0) >= 2]

venn.plot <- draw.pairwise.venn(
  area1 = length(asv_early_embryo),
  area2 = length(asv_late_embryo),
  cross.area = length(intersect(asv_early_embryo, asv_late_embryo)),
  category = c("Early Embryo", "Late Embryo"),
  fill = c("#519c3d", "tomato"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.2
)

grid.draw(venn.plot)









###### ANCOM BC2 on Early and Late in situ embryos ###########
BiocManager::install("ANCOMBC")
library(ANCOMBC)
library(tidyverse)
library(DT)
library(dplyr)
library(nloptr)

set.seed(123)
ps_emb_nosme <- subset_samples(ps_emb, Location != "SME")

###### LCF Embryogenesis + Location ########
ancom_emb2 <- ancombc2(ps_emb_nosme, 
                       taxa_are_rows = FALSE, 
                       tax_level = "Genus", 
                       fix_formula = "Embryogenesis",
                       pseudo = FALSE,
                       pseudo_sens = TRUE,
                       group = "Embryogenesis", 
                       global = FALSE, #only for >3 groups
                       pairwise = FALSE, #only for >3 groups
                       dunnet = FALSE,
                       alpha= 0.05, 
                       prv_cut = 0.1,
                       verbose = TRUE,
                       struc_zero = TRUE,
                       trend = TRUE)

## inspect results
res_prim2 <- ancom_emb2$res

ancom_emb2$res %>%
  dplyr::filter(diff_EmbryogenesisLate == "TRUE") %>%
  dplyr::select(taxon, starts_with("lfc_E"), starts_with("lfc_L"), starts_with("q_E"), starts_with("q_L")) %>%
  mutate_if(is.numeric, function(x) round(x, 2)) %>%
  datatable(caption = "Embryogenesis")




### plot lefc : Embryos ~ Embryogenesis
df_emb2 <- res_prim2 %>%
  dplyr::select(taxon, ends_with("Late")) 

df_fig_emb2 <- df_emb2 %>%
  dplyr::filter(diff_EmbryogenesisLate == "TRUE") %>% 
  dplyr::arrange(desc(lfc_EmbryogenesisLate)) %>%
  dplyr::mutate(direct = ifelse(lfc_EmbryogenesisLate > 0, "Enriched on Late", "Enriched on Early"))

df_fig_emb2$taxon <- factor(df_fig_emb2$taxon, levels = df_fig_emb2$taxon)

df_fig_emb2$direct <- factor(df_fig_emb2$direct, 
                             levels = c("Enriched on Late", "Enriched on Early"))

fig_emb2 <- df_fig_emb2 %>%
  ggplot(aes(x = taxon, y = lfc_EmbryogenesisLate, fill = direct)) + 
  geom_bar(
    stat = "identity",
    width = 0.7,
    color = "black",
    position = position_dodge(width = 0.4)) +
  geom_errorbar(
    aes(
      ymin = lfc_EmbryogenesisLate - se_EmbryogenesisLate,
      ymax = lfc_EmbryogenesisLate + se_EmbryogenesisLate),
    width = 0.2,
    position = position_dodge(0.05),
    color = "black") + 
  labs(
    x = NULL,
    y = "Log fold change",
    title = "Genera log fold change on late vs. early-stage embryos") +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Enriched on Early" = "#519c3d",
      "Enriched on Late"  = "#779fe5")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      color = df_fig_emb2$color))

fig_emb2















### plot lefc : Embryogenesis + Location ~ Location
df_emb3 <- res_prim2 %>%
  dplyr::select(taxon, ends_with("NH")) 
df_fig_emb3 <- df_emb3 %>%
  dplyr::filter(diff_LocationNH == "TRUE") %>% 
  dplyr::arrange(desc(lfc_LocationNH)) %>%
  dplyr::mutate(direct = ifelse(lfc_LocationNH > 0, "Enriched in NH", "Enriched in ME"),
                color = ifelse(diff_robust_LocationNH, "aquamarine3", "black"))
df_fig_emb3$taxon <- factor(df_fig_emb3$taxon, levels = df_fig_emb3$taxon)
df_fig_emb3$direct <- factor(df_fig_emb3$direct, 
                             levels = c("Enriched in NH", "Enriched in ME"))

fig_emb3 <- df_fig_emb3 %>%
  ggplot(aes(x = taxon, y = lfc_LocationNH, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_LocationNH - se_LocationNH, ymax = lfc_LocationNH + se_LocationNH), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Genera log fold change in NH vs. ME in all embryos") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig_emb3$color))
fig_emb3







######## LFC Embryogenesis * Location #########
ancom_emb <- ancombc2(ps_emb_nosme, 
                      taxa_are_rows = FALSE, 
                      tax_level = "Genus", 
                      fix_formula = "Embryogenesis * Location",
                      pseudo = TRUE,
                      pseudo_sens = TRUE,
                      group = "Location", 
                      global = FALSE, #only for >3 groups
                      pairwise = FALSE, #only for >3 groups
                      dunnet = FALSE,
                      alpha= 0.05, 
                      prv_cut = 0.1,
                      verbose = TRUE,
                      struc_zero = TRUE,
                      trend = TRUE)

## inspect results
res_prim <- ancom_emb$res %>%
  mutate(lfc_NH_vs_Maine_Late = lfc_LocationNH + `lfc_EmbryogenesisLate:LocationNH`) #add column for Late vs. Early in NH

ancom_emb$res %>%
  dplyr::filter(`diff_EmbryogenesisLate:LocationNH` == "TRUE" | diff_EmbryogenesisLate == "TRUE" |diff_LocationNH == "TRUE") %>%
  dplyr::select(taxon, starts_with("lfc_E"), starts_with("lfc_L"), starts_with("q_E"), starts_with("q_L")) %>%
  mutate_if(is.numeric, function(x) round(x, 2)) %>%
  datatable(caption = "Embryogenesis * Location")
  mutate(lfc_NH_vs_Maine_Late = lfc_LocationNH + lfc_EmbryogenesisLate:LocationNH)




### plot lefc : Embryogenesis ~ Maine (ie. late vs. early embryos in Maine)
df_emb <- res_prim %>%
  dplyr::select(taxon, ends_with("Late")) 
df_fig_emb <- df_emb %>%
  dplyr::filter(diff_EmbryogenesisLate == "TRUE") %>% 
  dplyr::arrange(desc(lfc_EmbryogenesisLate)) %>%
  dplyr::mutate(direct = ifelse(lfc_EmbryogenesisLate > 0, "Enriched on Late", "Enriched on Early"),
                color = ifelse(diff_robust_EmbryogenesisLate, "aquamarine3", "black"))
df_fig_emb$taxon <- factor(df_fig_emb$taxon, levels = df_fig_emb$taxon)
df_fig_emb$direct <- factor(df_fig_emb$direct, 
                           levels = c("Enriched on Late", "Enriched on Early"))

fig_emb <- df_fig_emb %>%
  ggplot(aes(x = taxon, y = lfc_EmbryogenesisLate, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_EmbryogenesisLate - se_EmbryogenesisLate, ymax = lfc_EmbryogenesisLate + se_EmbryogenesisLate), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Genera log fold change on late vs. early-stage embryos in Maine") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig_emb$color))
fig_emb







### plot lefc : Early Embryos ~ Location
df2 <- res_prim %>%
  dplyr::select(taxon, ends_with("_LocationNH")) 
df_fig2 <- df2 %>%
  dplyr::filter(diff_LocationNH == "TRUE") %>% 
  dplyr::arrange(desc(lfc_LocationNH)) %>%
  dplyr::mutate(direct = ifelse(lfc_LocationNH > 0, "Enriched in NH", "Enriched in ME"),
                color = ifelse(diff_robust_LocationNH, "aquamarine3", "black"))
df_fig2$taxon <- factor(df_fig2$taxon, levels = df_fig2$taxon)
df_fig2$direct <- factor(df_fig2$direct, 
                            levels = c("Enriched in NH", "Enriched in ME"))

fig_emb2 <- df_fig2 %>%
  ggplot(aes(x = taxon, y = lfc_LocationNH, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_LocationNH - se_LocationNH, ymax = lfc_LocationNH + se_LocationNH), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Genera log fold change in Maine vs. NH in early-stage embryos") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig2$color))
fig_emb2


df3 <- res_prim %>%
  select(taxon, lfc_NH_vs_Maine_Late | `diff_EmbryogenesisLate:LocationNH` | `q_EmbryogenesisLate:LocationNH`) %>%
  filter(`diff_EmbryogenesisLate:LocationNH`== "TRUE") #only one taxon, Sphingorhabdus, had -0.585 LFC with a interaction adj p value of 0.002, indicating this taxon was enriched in late embryos in New Hampshire


### how LFC between early and late changes relative to NH
emb_loc_dep <- res_prim %>%
  dplyr::select(taxon, `lfc_EmbryogenesisLate:LocationNH`, `diff_EmbryogenesisLate:LocationNH`, `q_EmbryogenesisLate:LocationNH`) %>%
  dplyr::filter(`diff_EmbryogenesisLate:LocationNH` == "TRUE") ### One genera, Sphingorhabdus, had a significant negative LFC difference of -4.41617 indicating this taxa was enriched in late Embryos in NH (adj p = 0.00229), and contrasts with Maine enrichment observed in early embryos












################################################################################
###################### PLEOPODS ################################################
################################################################################
phyob_pleo <- subset_samples(phyob, Sample_Type == "Pleo_Swab")
phyob_pleo <- prune_taxa(taxa_sums(phyob_pleo) > 0, phyob_pleo)

pleo_early <- subset_samples(phyob_pleo, Embryogenesis == "Early")
pleo_early <- prune_taxa(taxa_sums(pleo_early) >0, pleo_early)

pleo_late <- subset_samples(phyob_pleo, Embryogenesis == "Late")
pleo_late <- prune_taxa(taxa_sums(pleo_late) >0, pleo_late)


##### Alpha Diversity ######
rich.pleo <- estimate_richness(phyob_pleo, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(phyob_pleo)), by = "Sample")

kruskal.test(Shannon ~ Embryogenesis, data = rich.pleo) #chi-squared = 17.413, df = 1, p-value = 3.007e-05

ggboxplot(rich.pleo, x = "Embryogenesis", y = "Shannon", 
          color = "Embryogenesis", palette = c("#519c3d", "#c17fee"),
          order = c("Early", "Late"),
          ylab = "Alpha Diversity", xlab = "Sample Type")





#Early pleopod by location 
rich.pleo.earl <- estimate_richness(pleo_early, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(pleo_early)), by = "Sample")

kruskal.test(Shannon ~ Location, data = rich.pleo.earl) #chi-squared = 3.0582, df = 1, p-value = 0.08033

#Late pleopod by location
rich.pleo.late <- estimate_richness(pleo_late, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(pleo_late)), by = "Sample")

kruskal.test(Shannon ~ Location, data = rich.pleo.late) #chi-squared = 0.11797, df = 2, p-value = 0.9427
pairwise.wilcox.test(rich.pleo.late$Shannon, rich.pleo.late$Location, p.adjust.method = "BH")



### PCoA ###
min_lib_pleo <- min(sample_sums(phyob_pleo))

min_lib_pleo #47635
#subsample
set.seed(446)
pleo_rare <- rarefy_even_depth(phyob_pleo, sample.size = min_lib_pleo, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number

pleo_ord <- ordinate(physeq = pleo_rare, method ='PCoA', distance = 'bray')

pleo_PCOA <- plot_ordination(physeq = pleo_rare, ordination = pleo_ord, color='Embryogenesis', shape = 'Location') +
  geom_point(aes(color = Embryogenesis), alpha=1, size=4) +
  stat_ellipse(aes(group = Embryogenesis), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals9, name = "") +
  scale_shape_manual(values = vals10, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.15, 0.8),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

pleo_PCOA

vals9 <- c("Early" = "#c17fee", "Late" = "#3f78a6")
vals10 <- c('DownEast' = 16, 'NH' = 17, 'SME' = 15)
vals11 <- c("DownEast" = "#c17fee", "NH" = "#3f78a6", "SME" = "#519c3d")
vals12 <- c('Early' = 16, 'Late' = 17)
##### PERMANOVA ######
set.seed(47)

# Calculate bray curtis distance matrix
bray_pleo <- phyloseq::distance(phyob_pleo, method = "bray")

df_pleo <- data.frame(sample_data(phyob_pleo))
adon_pleo <- adonis2(bray_pleo ~ Embryogenesis, data = df_pleo) #p=0.001
#check dispersions
disp_pleo <- betadisper(bray_pleo, df_pleo$Embryogenesis)
anova_pleo <- anova(disp_pleo) #p=0.1779. dispersions are not different



### PCoA Early pleopods ###
min_lib_pleo_early <- min(sample_sums(pleo_early))

min_lib_pleo_early #48738
#subsample
set.seed(449)
pleo_rare_early <- rarefy_even_depth(pleo_early, sample.size = min_lib_pleo_early, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number

pleo_ord_early <- ordinate(physeq = pleo_rare_early, method ='PCoA', distance = 'bray')

pleo_PCOA_early <- plot_ordination(physeq = pleo_rare_early, ordination = pleo_ord_early, color='Location', shape = 'Location') +
  geom_point(aes(color = Location), alpha=1, size=4) +
  stat_ellipse(aes(group = Location), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals14, name = "") +
  scale_shape_manual(values = vals13, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.15, 0.8),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

pleo_PCOA_early

vals9 <- c("Early" = "#c17fee", "Late" = "#3f78a6")
vals10 <- c('DownEast' = 16, 'NH' = 17, 'SME' = 15)
vals11 <- c("DownEast" = "#c17fee", "NH" = "#3f78a6", "SME" = "#519c3d")
vals12 <- c('Early' = 16, 'Late' = 17)
vals13 <- c('DownEast' = 16, 'NH' = 17)
vals14 <- c("DownEast" = "#c17fee", "NH" = "#3f78a6")

##### PERMANOVA ######
set.seed(48)

# Calculate bray curtis distance matrix
bray_pleo_early <- phyloseq::distance(pleo_early, method = "bray")
df_pleo_early <- data.frame(sample_data(pleo_early))
adon_pleo_early <- adonis2(bray_pleo_early ~ Location, data = df_pleo_early) #p=0.001

#check dispersions
disp_pleo_earl <- betadisper(bray_pleo_early, df_pleo_early$Location)
anova_pleo_earl <- anova(disp_pleo_earl) #p=2.749e-05. dispersions are different
plot(disp_pleo_earl)






### PCoA Late pleopods ###
min_lib_pleo_late <- min(sample_sums(pleo_late))

min_lib_pleo_late #53252
#subsample
set.seed(449)
pleo_rare_late <- rarefy_even_depth(pleo_late, sample.size = min_lib_pleo_late, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number

pleo_ord_late <- ordinate(physeq = pleo_rare_late, method ='PCoA', distance = 'bray')

pleo_PCOA_late <- plot_ordination(physeq = pleo_rare_late, ordination = pleo_ord_late, color='Location', shape = 'Location') +
  geom_point(aes(color = Location), alpha=1, size=4) +
  stat_ellipse(aes(group = Location), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals11, name = "") +
  scale_shape_manual(values = vals10, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.15, 0.8),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

pleo_PCOA_late

vals9 <- c("Early" = "#c17fee", "Late" = "#3f78a6")
vals10 <- c('DownEast' = 16, 'NH' = 17, 'SME' = 15)
vals11 <- c("DownEast" = "#c17fee", "NH" = "#3f78a6", "SME" = "#519c3d")
vals12 <- c('Early' = 16, 'Late' = 17)
vals13 <- c('DownEast' = 16, 'NH' = 17)
vals14 <- c("DownEast" = "#c17fee", "NH" = "#3f78a6")

##### PERMANOVA ######
set.seed(49)

# Calculate bray curtis distance matrix
pleo_late_nosme <- subset_samples(pleo_late, Location != "SME")
bray_pleo_late <- phyloseq::distance(pleo_late_nosme, method = "bray")

# make a data frame from the sample_data
df_pleo_late <- data.frame(sample_data(pleo_late_nosme))
adon_pleo_late <- adonis2(bray_pleo_late ~ Location, data = df_pleo_late) #p=0.032

#check dispersions
disp_pleo_late <- betadisper(bray_pleo_late, df_pleo_late$Location)
anova_pleo_late <- anova(disp_pleo_late) #p=0.8339. dispersions are not different
plot(disp_pleo_late)

#### Shared genera of pleopods early v late ####
# Remove taxa without genus assignment
pleo_ps_genus <- subset_taxa(phyob_pleo, !is.na(Genus) & Genus != "") #remove asvs where genera is na, 4672 asvs left
pleo_ps_genus <- tax_glom(pleo_ps_genus, taxrank = "Genus", NArm = TRUE) #agglomerate to genus, 744 genera

pleo_late_genus <- subset_samples(pleo_ps_genus, Embryogenesis == "Late")
pleo_late_genus <- prune_taxa(taxa_sums(pleo_late_genus) > 0, pleo_late_genus) #658 genera

pleo_early_genus <- subset_samples(pleo_ps_genus, Embryogenesis == "Early")
pleo_early_genus <- prune_taxa(taxa_sums(pleo_early_genus) > 0, pleo_early_genus) #445

pleo_early_taxa <- taxa_names(pleo_early_genus)
pleo_late_taxa <- taxa_names(pleo_late_genus)

asv_pleo_earl <- otu_table(pleo_early_genus)
asv_pleo_late <- otu_table(pleo_late_genus)

pleo_early_genus <- taxa_names(pleo_early_genus)[rowSums(asv_pleo_earl > 0) >= 2] #present in at least 2 samples, 452
pleo_late_genus <- taxa_names(pleo_late_genus)[rowSums(asv_pleo_late > 0) >= 2] #364


#### Make Venn Diagram
draw.pairwise.venn(
  area1 = length(pleo_early_genus),
  area2 = length(pleo_late_genus),
  cross.area = length(intersect(pleo_early_genus, pleo_late_genus)),
  category = c("Early Pleopods", "Late Pleopods"),
  fill = c("steelblue", "tomato"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.2)

## make a list of shared genera
genus_earl_list <- as.character(tax_table(emb_earl)[, "Genus"])
genus_late_list <- as.character(tax_table(emb_late)[, "Genus"])
intersect(genus_earl_list, genus_late_list)






### ANCOMBC to check LFC differences in pleopods early to late ###

ancom_pleo <- ancombc2(phyob_pleo, 
                               taxa_are_rows = FALSE, 
                               tax_level = "Genus", 
                               fix_formula = "Embryogenesis + Location",
                               group = "Embryogenesis", 
                               #global = TRUE, #only for >3 groups
                               #pairwise = TRUE, #only for >3 groups
                               alpha= 0.05, 
                               prv_cut = 0.1,
                               verbose = TRUE,
                               trend = TRUE)

res_pleo <- ancom_pleo$res
df_pleo <- res_pleo %>%
  dplyr::select(taxon, ends_with("Late")) 

save(df_pleo, file = "ancombc_pleo.csv")


## make LEfse plot
df_fig_pleo <- df_pleo %>%
  dplyr::filter(diff_EmbryogenesisLate == "TRUE",
                lfc_EmbryogenesisLate > 1.5 | lfc_EmbryogenesisLate < -1.5) %>% 
  dplyr::arrange(desc(lfc_EmbryogenesisLate)) %>%
  dplyr::mutate(direct = ifelse(lfc_EmbryogenesisLate > 0, "Enriched on Late", "Enriched on Early"),
                color = ifelse(diff_robust_EmbryogenesisLate, "aquamarine3", "black"))

df_fig_pleo$taxon <- factor(df_fig_pleo$taxon, levels = df_fig_pleo$taxon)
df_fig_pleo$direct <- factor(df_fig_pleo$direct, 
                            levels = c("Enriched on Late", "Enriched on Early"))

fig_pleo <- df_fig_pleo %>%
  ggplot(aes(x = taxon, y = lfc_EmbryogenesisLate, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_EmbryogenesisLate - se_EmbryogenesisLate, ymax = lfc_EmbryogenesisLate + se_EmbryogenesisLate), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Genera log fold change with > |1.5| on late-stage pleopods") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig_pleo$color))
fig_pleo









############## Embryos and Pleopod Comparison ####################
pleo_emb <- c("Pleopod", "Embryo")
ps_emb_pleo <- subset_samples(phyob, Sample_Kind == pleo_emb)
ps_emb_pleo <- prune_taxa(taxa_sums(ps_emb_pleo) >0, ps_emb_pleo)




## PCOA ##
ml_pleo_emb <- min(sample_sums(ps_emb_pleo))

ml_pleo_emb #50108
#subsample
set.seed(543)
pleo_emb_rare <- rarefy_even_depth(ps_emb_pleo, sample.size = ml_pleo_emb, verbose=FALSE, replace=TRUE) 

pleo_emb_ord <- ordinate(physeq = pleo_emb_rare, method ='PCoA', distance = 'bray')

pleo_emb_PCOA <- plot_ordination(physeq = pleo_emb_rare, ordination = pleo_emb_ord, color='Sample_Type', shape = 'Embryogenesis') +
  geom_point(aes(color = Sample_Type), alpha=1, size=4) +
  stat_ellipse(aes(group = Sample_Type), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals15, name = "") +
  scale_shape_manual(values = vals16, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.9, 0.15),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

vals15 <- c("Late_Embryo" = "#c17fee", "Early_Embryo" = "#3f78a6", "Pleo_Swab" ="#519c3d")
vals16 <- c('Early' = 16, 'Late' = 17)


pleo_emb_PCOA


### early pleo + embryo LFC x Location ###
early_pe <- subset_samples(ps_emb_pleo, Embryogenesis == "Early")
early_pe <- prune_taxa(taxa_sums(early_pe) >0, early_pe)

ancom_early_pe <- ancombc2(early_pe, 
                           taxa_are_rows = FALSE, 
                           tax_level = "Genus", 
                           fix_formula = "Sample_Kind",
                           pseudo = TRUE,
                           pseudo_sens = TRUE,
                           group = "Sample_Kind", 
                           global = FALSE, #only for >3 groups
                           pairwise = FALSE, #only for >3 groups
                           dunnet = FALSE,
                           alpha= 0.05, 
                           prv_cut = 0.1,
                           verbose = TRUE,
                           struc_zero = TRUE,
                           trend = TRUE)

res_early_pleo_pe <- ancom_early_pe$res

df_early_pleo_pe <- res_early_pleo_pe %>%
  dplyr::filter(diff_Sample_KindPleopod == "TRUE") %>%
  dplyr::select(taxon, ends_with("Pleopod"))  

df_type_early <- df_early_pleo_pe  %>%
  dplyr::select(taxon, ends_with("Pleopod")) %>%
  dplyr::arrange(desc(lfc_Sample_KindPleopod)) %>%
  dplyr::mutate(direct = ifelse(lfc_Sample_KindPleopod > 0, "Enriched on Pleopods", "Enriched on Embryos"))

df_type_early$taxon <- factor(df_type_early$taxon, levels = df_type_early$taxon)
df_type_early$direct <- factor(df_type_early$direct, 
                         levels = c("Enriched on Pleopods", "Enriched on Embryos"))

fig_type_early <- df_type_early %>%
  ggplot(aes(x = taxon, y = lfc_Sample_KindPleopod, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_Sample_KindPleopod - se_Sample_KindPleopod, ymax = lfc_Sample_KindPleopod + se_Sample_KindPleopod), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Early pleopods vs embryos") + 
  scale_fill_manual(values = direct_colors,
                    name = NULL) +
  #scale_fill_discrete(name = NULL) +
  #scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_early_pleo_pe$color))
fig_type_early














plot_df <- res_early_pleo_pe %>%
  filter(`q_Sample_KindPleopod:LocationNH` < 0.05) %>%
  transmute(taxon,
    Maine = lfc_Sample_KindPleopod,
    `New Hampshire` = lfc_Sample_KindPleopod + `lfc_Sample_KindPleopod:LocationNH`) %>%
  pivot_longer(
    cols = c(Maine, `New Hampshire`),
    names_to = "Location",
    values_to = "LFC")

early_pleo_pe_strat_LFC <- ggplot(plot_df, aes(x = Location, y = LFC, group = taxon, color = taxon)) +
  geom_line(alpha = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Location",
    y = "Log fold change (Pleopod vs Embryo)",
    title = "Location-dependent Early Pleopod–Embryo effects across genera") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )

##Caption: Colored interaction plot showing early pleopod–embryo log fold change estimates by location for genera exhibiting significant Sample_Kind × Location interactions (FDR < 0.05). Each line represents one genus.

## Interpretation: Seven genera exhibit strong pleopod enrichment relative to early embryos in Maine, but this contrast is substantially reduced or reversed in New Hampshire. One exception is the genus Alkanindiges, where pleopod enrichment is observed in New Hampshire relative to Maine.

##### PERMANOVA early Emb + pleo ######
set.seed(90)

# Calculate bray curtis distance matrix
bray_early_pe <- phyloseq::distance(early_pe, method = "bray")

# make a data frame from the sample_data
df_early_pe <- data.frame(sample_data(early_pe))
adon_early_pe <- adonis2(bray_early_pe ~ Sample_Type, data = df_early_pe) #p=0.001

#check dispersions
disp_early_pe <- betadisper(bray_early_pe, df_early_pe$Sample_Type)
anova_early_pe <- anova(disp_early_pe) #p =0.001783
plot(disp_early_pe) # dispersions don't look different









### LFC Late embryos + Pleopods ###
late_pe <- subset_samples(ps_emb_pleo, Embryogenesis == "Late" & Location != "SME")
late_pe <- prune_taxa(taxa_sums(late_pe) > 0, late_pe)
late_pe_sub <- subset_taxa(late_pe, Genus != "Unknown")

ancom_late_pe <- ancombc2(late_pe_sub, 
                           taxa_are_rows = FALSE, 
                           tax_level = "Genus", 
                           fix_formula = "Sample_Kind",
                           pseudo = TRUE,
                           pseudo_sens = TRUE,
                           group = "Sample_Kind",
                           global = FALSE, #only for >3 groups
                           pairwise = FALSE, #only for >3 groups
                           dunnet = FALSE,
                           alpha= 0.05, 
                           prv_cut = 0.1,
                           verbose = TRUE,
                           struc_zero = TRUE,
                           trend = TRUE)

res_late_pleo_pe <- ancom_late_pe$res

df_late_pleo_pe <- res_late_pleo_pe %>%
  dplyr::filter(diff_Sample_KindPleopod == "TRUE") %>%
  dplyr::select(taxon, ends_with("Pleopod"))  

df_type_late <- df_late_pleo_pe  %>%
  dplyr::select(taxon, ends_with("Pleopod")) %>%
  dplyr::arrange(desc(lfc_Sample_KindPleopod)) %>%
  dplyr::mutate(direct = ifelse(lfc_Sample_KindPleopod > 0, "Enriched on Pleopods", "Enriched on Embryos"))

df_type_late$taxon <- factor(df_type_late$taxon, levels = df_type_late$taxon)
df_type_late$direct <- factor(df_type_late$direct, 
                         levels = c("Enriched on Pleopods", "Enriched on Embryos"))

direct_colors <- c(
  "Enriched on Pleopods" = "#3f78a6",
  "Enriched on Embryos"  = "#519c3d")

fig_type_late <- df_type_late %>%
  ggplot(aes(x = taxon, y = lfc_Sample_KindPleopod, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_Sample_KindPleopod - se_Sample_KindPleopod, ymax = lfc_Sample_KindPleopod + se_Sample_KindPleopod), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Late pleopods vs embryos") + 
  scale_fill_manual(values = direct_colors,
                    name = NULL) +
  #scale_fill_discrete(name = NULL) +
  #scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_late_pleo_pe$color))
fig_type_late


library(ggpubr)
combined <- ggarrange(fig_type_early, fig_type_late,
                    ncol = 2, nrow = 1,  # Arrange in 1 row and 2 columns
                    labels = c("A", "B"), # Add labels (optional)
                    common.legend = TRUE, # Use a common legend for all plots
                    legend = "bottom"     # Position the common legend at the bottom
)


##### PERMANOVA late Emb + pleo ######
set.seed(93)

# Calculate bray curtis distance matrix
bray_late_pe <- phyloseq::distance(late_pe, method = "bray")

# make a data frame from the sample_data
df_late_pe <- data.frame(sample_data(late_pe))
adon_late_pe <- adonis2(bray_late_pe ~ Sample_Type, data = df_late_pe) #p=0.001

#check dispersions
disp_late_pe <- betadisper(bray_late_pe, df_late_pe$Sample_Type)
anova_late_pe <- anova(disp_late_pe) #p =0.007748
plot(disp_late_pe) # dispersions don't look different







##### PERMANOVA ######
set.seed(57)

# Calculate bray curtis distance matrix
bray_pleo_emb <- phyloseq::distance(ps_emb_pleo, method = "bray")

# make a data frame from the sample_data
df_pleo_emb <- data.frame(sample_data(ps_emb_pleo))
adon_pleo_emb <- adonis2(bray_pleo_emb ~ Sample_Type, data = df_pleo_emb) #p=0.001

#check dispersions
disp_pleo_emb <- betadisper(bray_pleo_emb, df_pleo_emb$Sample_Type)
anova_pleo_emb <- anova(disp_pleo_emb) #p is very small
plot(disp_pleo_emb) # dispersions don't look different

#### Shared genera of Pleopods and Embryos
# Remove taxa without genus assignment
pleo_ps_genus 
emb_ps_genus 

taxa_names_emb <- taxa_names(emb_ps_genus)
taxa_names_pleo <- taxa_names(pleo_ps_genus)

asv_emb <- otu_table(emb_ps_genus)
asv_pleo <- otu_table(pleo_ps_genus)

genera_emb <- taxa_names(emb_ps_genus)[rowSums(asv_emb > 0) >= 2] #present in at least 2 samples, 452
genera_pleo <- taxa_names(pleo_ps_genus)[rowSums(asv_pleo > 0) >= 2] #364


#### Make Venn Diagram

draw.pairwise.venn(
  area1 = length(genera_emb),
  area2 = length(genera_pleo),
  cross.area = length(intersect(genera_emb, genera_pleo)),
  category = c("Embryos", "Pleopods"),
  fill = c("steelblue", "tomato"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.2
)

## make a list of shared genera
genus_earl_list <- as.character(tax_table(emb_earl)[, "Genus"])
genus_late_list <- as.character(tax_table(emb_late)[, "Genus"])
intersect(genus_earl_list, genus_late_list)







################################################################################
###################### Shell ################################################
################################################################################
phyob_shell <- subset_samples(phyob, Sample_Type != "Early_Embryo" & Sample_Type != "Late_Embryo" & Sample_Type != "Hemolymph" & Sample_Type != "Pleo_Swab")
phyob_shellnosme <- subset_samples(phyob_shell, Location != "SME")
phyob_shell <- prune_taxa(taxa_sums(phyob_shell) > 0, phyob_shell)

shell_early <- subset_samples(phyob_shell, Embryogenesis == "Early")
shell_early <- prune_taxa(taxa_sums(shell_early) > 0, shell_early)

shell_late <- subset_samples(phyob_shell, Embryogenesis =="Late")
shell_late <- prune_taxa(taxa_sums(shell_late) >0, shell_late)

###Alpha Diversity Shell ~ Embryogenesis
rich.shell <- estimate_richness(phyob_shell, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(phyob_shell)), by = "Sample")

kruskal.test(Shannon ~ Embryogenesis, data = rich.shell) # not different by embryogenesis: chi-squared: 3.0155, p-value=0.082

kruskal.test(Shannon ~ Location, data = rich.shell) # not different by Location: chi-squared: 0.65572, p-value=0.7205

#Early Shell by location 
rich.shell.earl <- estimate_richness(shell_early, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(shell_early)), by = "Sample")

kruskal.test(Shannon ~ Location, data = rich.shell.earl) #not enough Downeast samples to reliably calculate

#Late shell by location
rich.shell.late <- estimate_richness(shell_late, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(shell_late)), by = "Sample")

kruskal.test(Shannon ~ Location, data = rich.shell.late) 

###### PCOA #######
min_lib_shell <- min(sample_sums(phyob_shellnosme))

min_lib_shell #3602
#subsample
set.seed(43)
shell_rare <- rarefy_even_depth(phyob_shellnosme, sample.size = min_lib_shell, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number
shell_rare #new ps object with rarefied and subsampled ASVs 

shell_ord <- ordinate(physeq = shell_rare, method ='PCoA', distance = 'bray')

shell_PCOA <- plot_ordination(physeq = shell_rare, ordination = shell_ord, color='Embryogenesis', shape = 'Location') +
  geom_point(aes(color = Embryogenesis), alpha=1, size=4) +
  stat_ellipse(aes(group = Embryogenesis), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals9, name = "") +
  scale_shape_manual(values = vals10, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

shell_PCOA

vals9 <- c("Early" = "#c17fee", "Late" = "#3f78a6")
vals10 <- c('DownEast' = 16, 'NH' = 17, 'SME' = 15)





##### PERMANOVA - All shell ######
set.seed(42)

# Calculate bray curtis distance matrix
bray_shell <- phyloseq::distance(phyob_shellnosme, method = "bray")

# make a data frame from the sample_data
df_shell <- data.frame(sample_data(phyob_shellnosme))
adon_shell <- adonis2(bray_shell ~ Embryogenesis, data = df_shell) #p=0.001

#Check dispersions
disp_shell <- betadisper(bray_shell, df_shell$Embryogenesis)
anova_shell <- anova(disp_shell) #p=0.01128. dispersions are different
plot(disp_shell)



### PCoA late shell only ###
ml_shell_late <- min(sample_sums(shell_late2))

ml_shell_late #14663
#subsample
set.seed(44)
late_shell_rare <- rarefy_even_depth(shell_late2, sample.size = ml_shell_late, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number
late_shell_rare #new ps object with rarefied and subsampled ASVs 

late_shell_ord <- ordinate(physeq = late_shell_rare, method ='PCoA', distance = 'bray')

shell_PCOA_late <- plot_ordination(physeq = late_shell_rare, ordination = late_shell_ord, color='Location', shape = 'ESD_Status') +
  geom_point(aes(color = Location), alpha=1, size=4) +
  #stat_ellipse(aes(group = Location), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals13, name = "") +
  scale_shape_manual(values = vals12, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

shell_PCOA_late

vals13 <- c("DownEast" = "#c17fee", "NH" = "#3f78a6")
vals12 <- c('Positive' = 16, 'Negative' = 17)

### PERMANOVA Late shell ###
#remove one SME sample - cannot perform stats on one 
shell_late2 <- subset_samples(shell_late, Sample != "L_ME1_cswab")
set.seed(45)
bray_shell_late <- phyloseq::distance(shell_late2, method = "bray")

df_shell_late <- data.frame(sample_data(shell_late2))
adon_shell_late <- adonis2(bray_shell_late ~ ESD_Status * Location, data = df_shell_late) #location pval = 0.027, esd_status pval = 0.003
library(pairwiseAdonis)
pairwise.adonis2(bray_shell_late ~ ESD_Status*Location, data=df_shell_late)


#Check dispersions
disp_shell_late <- betadisper(bray_shell_late, df_shell_late$Location)
anova_shell_late <- anova(disp_shell_late) #p=0.95. dispersions are not different
plot(disp_shell_late)

### ESD_Status * Location was significant, indicating the effect of ESD_Status differs by location.
#Check NH
shell_lateNH <- subset_samples(shell_late2, Location =="NH")
set.seed(87)
bray_shell_lateNH <- phyloseq::distance(shell_lateNH, method = "bray")
df_shell_lateNH <- data.frame(sample_data(shell_lateNH))
adon_shell_lateNH <- adonis2(bray_shell_lateNH ~ ESD_Status, data = df_shell_lateNH) # p = 0.002
#Check dispersions
disp_shell_lateNH <- betadisper(bray_shell_lateNH, df_shell_lateNH$ESD_Status)
anova_shell_lateNH <- anova(disp_shell_lateNH) #p=0.3194. dispersions are not different
plot(disp_shell_lateNH)

#Check ME
shell_lateME <- subset_samples(shell_late2, Location =="DownEast")
set.seed(87)
bray_shell_lateME <- phyloseq::distance(shell_lateME, method = "bray")
df_shell_lateME <- data.frame(sample_data(shell_lateME))
adon_shell_lateME <- adonis2(bray_shell_lateME ~ ESD_Status, data = df_shell_lateME) # p = 0.288
#Check dispersions
disp_shell_lateME <- betadisper(bray_shell_lateME, df_shell_lateME$ESD_Status)
anova_shell_lateME <- anova(disp_shell_lateME) #p=0.8374. dispersions are not different
plot(disp_shell_lateME)

###PCOA Late Shell NH ####
ml_shell_lateNH <- min(sample_sums(shell_lateNH))

ml_shell_lateNH #14663
#subsample
set.seed(77)
late_shell_rareNH <- rarefy_even_depth(shell_lateNH, sample.size = ml_shell_lateNH, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number
late_shell_rareNH #new ps object with rarefied and subsampled ASVs 

late_shell_ordNH <- ordinate(physeq = late_shell_rareNH, method ='PCoA', distance = 'bray')

shell_PCOA_lateNH <- plot_ordination(physeq = late_shell_rareNH, ordination = late_shell_ordNH, color='Location', shape = 'ESD_Status') +
  geom_point(aes(color = Location), alpha=1, size=4) +
  stat_ellipse(aes(group = ESD_Status), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals13, name = "") +
  scale_shape_manual(values = vals12, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.2, 0.85),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

shell_PCOA_lateNH

vals13 <- c("DownEast" = "#c17fee", "NH" = "#3f78a6")
vals12 <- c('Positive' = 16, 'Negative' = 17)



###PCOA Late Shell ME ####
ml_shell_lateME <- min(sample_sums(shell_lateME))

ml_shell_lateME #70319
#subsample
set.seed(78)
late_shell_rareME <- rarefy_even_depth(shell_lateME, sample.size = ml_shell_lateME, verbose=FALSE, replace=TRUE) #randomly subsample your ASVs so that all samples have the same number
late_shell_rareME #new ps object with rarefied and subsampled ASVs 

late_shell_ordME <- ordinate(physeq = late_shell_rareME, method ='PCoA', distance = 'bray')

shell_PCOA_lateME <- plot_ordination(physeq = late_shell_rareME, ordination = late_shell_ordME, color='Location', shape = 'ESD_Status') +
  geom_point(aes(color = Location), alpha=1, size=4) +
  stat_ellipse(aes(group = ESD_Status), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals13, name = "") +
  scale_shape_manual(values = vals12, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.2, 0.85),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

shell_PCOA_lateME

vals13 <- c("DownEast" = "#c17fee", "NH" = "#3f78a6")
vals12 <- c('Positive' = 16, 'Negative' = 17)




### ANCOMBC-2 of Shell early vs late ###
ancom_shell <- ancombc2(phyob_shellnosme, 
                       taxa_are_rows = FALSE, 
                       tax_level = "Genus", 
                       fix_formula = "Embryogenesis + Location + ESD_Status",
                       group = "Embryogenesis", 
                       #global = TRUE, #only for >3 groups
                       #pairwise = TRUE, #only for >3 groups
                       alpha= 0.05, 
                       prv_cut = 0.1,
                       verbose = TRUE,
                       trend = TRUE)

res_shell <- ancom_shell$res
df_shell_ancom <- res_shell %>%
  dplyr::select(taxon, ends_with("Late"), ends_with("NH"), ends_with("SME"), ends_with("Positive")) 

save(df_shell, file = "ancombc_shell.csv")

##make lefse
df_fig_shell <- df_shell %>%
  dplyr::filter(diff_EmbryogenesisLate == "TRUE")%>% 
  dplyr::arrange(desc(lfc_EmbryogenesisLate)) %>%
  dplyr::mutate(direct = ifelse(lfc_EmbryogenesisLate > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(diff_robust_EmbryogenesisLate, "aquamarine3", "black"))

df_fig_shell$taxon <- factor(df_fig_shell$taxon, levels = df_fig_shell$taxon)
df_fig_shell$direct <- factor(df_fig_shell$direct, 
                             levels = c("Positive LFC", "Negative LFC"))

fig_shell <- df_fig_shell %>%
  ggplot(aes(x = taxon, y = lfc_EmbryogenesisLate, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_EmbryogenesisLate - se_EmbryogenesisLate, ymax = lfc_EmbryogenesisLate + se_EmbryogenesisLate), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Genera log fold change on late-stage carapace controlled for location & ESD status") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig_pleo$color))
fig_shell

### ANCOMBC-2 of early shell x ESD_Status ###
early_shell_esd <- ancombc2(shell_early, 
                            taxa_are_rows = FALSE, 
                            tax_level = "Genus", 
                            fix_formula = "ESD_Status + Location",
                            group = "ESD_Status", 
                            #global = TRUE, #only for >3 groups
                            #pairwise = TRUE, #only for >3 groups
                            pseudo = TRUE,
                            pseudo_sens = TRUE,
                            alpha= 0.05, 
                            prv_cut = 0.1,
                            verbose = TRUE,
                            trend = TRUE)

early_shell_res <- early_shell_esd$res

early_shell_df <- early_shell_res %>%
  dplyr::filter(diff_ESD_StatusPositive == "TRUE") %>%
  select(taxon, lfc_ESD_StatusPositive, q_ESD_StatusPositive)


### ANCOMBC-2 of late shell x ESD_Status ###
late_shell_esd <- ancombc2(shell_late2, 
                           taxa_are_rows = FALSE, 
                           tax_level = "Genus", 
                           fix_formula = "ESD_Status * Location",
                           group = "ESD_Status", 
                           #global = TRUE, #only for >3 groups
                           #pairwise = TRUE, #only for >3 groups
                           pseudo = TRUE,
                           pseudo_sens = TRUE,
                           alpha= 0.05, 
                           prv_cut = 0.1,
                           verbose = TRUE,
                           trend = TRUE)

late_shell_esd <- late_shell_esd$res

late_shell_df <- late_shell_esd %>%
  dplyr::filter(diff_ESD_StatusPositive == "TRUE" | diff_LocationNH == "TRUE" | `diff_ESD_StatusPositive:LocationNH` == "TRUE") %>%
  select(taxon, lfc_ESD_StatusPositive, q_ESD_StatusPositive, lfc_LocationNH, q_LocationNH, `lfc_ESD_StatusPositive:LocationNH`, `q_ESD_StatusPositive:LocationNH`)

##check aquimarina
aqui <- late_shell_esd %>%
  filter(taxon == "Aquimarina")%>%
  select(taxon, lfc_ESD_StatusPositive, q_ESD_StatusPositive)










################################################################################
###################### Hemolymph ################################################
################################################################################

phyob_hem <- subset_samples(phyob, Sample_Type == "Hemolymph")
phyob_hem <- prune_taxa(taxa_sums(phyob_hem) > 0, phyob_hem)
phyob_hem2 <- subset_samples(phyob_hem, Location != "SME")

###Alpha Diversity Hem ~ Location
rich.hem <- estimate_richness(phyob_hem2, split = TRUE, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var= "Sample") %>%
  left_join(as.data.frame(sample_data(phyob_hem2)), by = "Sample")

kruskal.test(Shannon ~ Location, data = rich.hem) # not different by Location: chi-squared: 0.72985, p-value=0.6942

### PCoA Hemolymph ###
ml_hem <- min(sample_sums(phyob_hem2))

ml_hem #34574
#subsample
set.seed(47)
hem_rare <- rarefy_even_depth(phyob_hem2, sample.size = ml_hem, verbose=FALSE, replace=TRUE) 
hem_rare  

hem_ord <- ordinate(physeq = hem_rare, method ='PCoA', distance = 'bray')

hem_PCOA <- plot_ordination(physeq = hem_rare, ordination = hem_ord, color='Location', shape = 'Location') +
  geom_point(aes(color = Location), alpha=1, size=4) +
  stat_ellipse(aes(group = Location), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals13, name = "") +
  scale_shape_manual(values = vals12, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

hem_PCOA

vals13 <- c("DownEast" = "#c17fee", "NH" = "#3f78a6")
vals12 <- c('DownEast' = 16, 'NH' = 17)




### PERMANOVA Hemolymph ###
set.seed(49)
bray_hem <- phyloseq::distance(phyob_hem2, method = "bray")

df_hem <- data.frame(sample_data(phyob_hem2))
adon_hem <- adonis2(bray_hem ~ ESD_Status, data = df_hem) #pval = 0.095 

### Bar Plot
phyob_hem_bp <- phyob_hem %>%
  tax_glom(taxrank='Genus') %>% #select/combine at specified taxa level
  transform_sample_counts(function(x) {(x/sum(x))*100} ) %>% #convert to relative abundance
  psmelt() %>% #switch to ggplot usable format
  filter(Abundance > 5) %>% #filter by taxa with relative abundance >1%
  arrange(Genus) #sort alphabetically 

hem_bp <- ggplot(phyob_hem_bp, aes(x=Sample, y=Abundance, fill=Genus)) + 
  #facet_grid(~factor(Sample_Kind, levels=c("Embryo", "Pleopod", "Carapace", "Hemolymph")), scales = "free", space = "free")+ #makes panels on your plot
  guides(x = guide_axis(angle = 90)) + #sample names turn 90
  geom_bar(stat='identity', color = "black") + # black around bars
  scale_fill_manual(values = colors_125) + #colors set below
  guides(fill = guide_legend(reverse = FALSE, keywidth = 2, keyheight = 1, ncol = 1)) + 
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 11))+
  ylab('Relative Genus Abundance (>5%) \n') +
  xlab('Sample ID')

hem_bp

library(plotly)
ggplotly(hem_bp)


################### Shell and Hemolymph comparison ######################
ps_shell_hem <- merge_phyloseq(phyob_shell, phyob_hem)
ps_shell_hem <- prune_taxa(taxa_sums(ps_shell_hem)>0, ps_shell_hem)


## PCOA ##
ml_shell_hem <- min(sample_sums(ps_shell_hem))

ml_shell_hem #14663
#subsample
set.seed(543)
shell_hem_rare <- rarefy_even_depth(ps_shell_hem, sample.size = ml_shell_hem, verbose=FALSE, replace=TRUE) 

shell_hem_ord <- ordinate(physeq = shell_hem_rare, method ='PCoA', distance = 'bray')

shell_hem_PCOA <- plot_ordination(physeq = shell_hem_rare, ordination = shell_hem_ord, color='Type', shape = 'Type') +
  geom_point(aes(color = Type), alpha=1, size=4) +
  stat_ellipse(aes(group = Type), type = "norm", level = 0.95) +
  #ggtitle("Principal Coordinate Analysis", subtitle = "Early & late-stage embryos, larvae, and tank water") +
  scale_color_manual(values = vals15, name = "") +
  scale_shape_manual(values = vals16, name = "") + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size=20),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=13),
        legend.position = c(0.9, 0.15),
        legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95")) 

vals15 <- c("Carapace" = "#c17fee", "Hemolymph" = "#3f78a6")
vals16 <- c('Carapace' = 16, 'Hemolymph' = 17)



### Shared genera of carapace and hemolymph
tax_table_cara <- tax_table(phyob_shell)
tax_table_hem <- tax_table(phyob_hem)
tax_table_hem_cara <- tax_table(ps_shell_hem)

# Filter to get only genera (level 6 in the taxonomy table)
genus_cara <- tax_table_cara[, "Genus"]
genus_hem <- tax_table_hem[, "Genus"]
genus_shell_hem <- tax_table_hem_cara[, "Genus"]

# Get the taxa present in each group (removing any NA values)
genus_cara <- na.omit(genus_cara)
genus_hem <- na.omit(genus_hem)
genus_shell_hem <- na.omit(genus_shell_hem)

# Find shared genera
shared_genera <- intersect(genus_cara, genus_hem) #share 323 genera of 5771






### Carapace and Hemolymph LFC ###
ancom_cara_hem <- ancombc2(ps_shell_hem, 
                           taxa_are_rows = FALSE, 
                           tax_level = "Genus", 
                           fix_formula = "Type",
                           pseudo = TRUE,
                           pseudo_sens = TRUE,
                           group = "Type", 
                           global = FALSE, #only for >3 groups
                           pairwise = FALSE, #only for >3 groups
                           dunnet = FALSE,
                           alpha= 0.05, 
                           prv_cut = 0.1,
                           verbose = TRUE,
                           struc_zero = TRUE,
                           trend = TRUE)

res_cara_hem <- ancom_cara_hem$res

df_cara_hem <- res_cara_hem %>%
  dplyr::filter(diff_TypeHemolymph == "TRUE") %>%
  dplyr::select(taxon, ends_with("Hemolymph"))  

df_type_cara_hem <- df_cara_hem  %>%
  dplyr::select(taxon, ends_with("Hemolymph")) %>%
  dplyr::arrange(desc(lfc_TypeHemolymph)) %>%
  dplyr::mutate(direct = ifelse(lfc_TypeHemolymph > 0, "Enriched on Hemolymph", "Enriched on Carapace"))

df_type_cara_hem$taxon <- factor(df_type_cara_hem$taxon, levels = df_type_cara_hem$taxon)
df_type_cara_hem$direct <- factor(df_type_cara_hem$direct, 
                               levels = c("Enriched on Hemolymph", "Enriched on Carapace"))
direct_colors2 <- c(
  "Enriched on Hemolymph" = "#3f78a6",
  "Enriched on Carapace"  = "#519c3d")

fig_type_cara_hem <- df_type_cara_hem %>%
  ggplot(aes(x = taxon, y = lfc_TypeHemolymph, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_TypeHemolymph - se_TypeHemolymph, ymax = lfc_TypeHemolymph + se_TypeHemolymph), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Hemolymph vs Carapace") + 
  scale_fill_manual(values = direct_colors2,
                    name = NULL) +
  #scale_fill_discrete(name = NULL) +
  #scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_cara_hem$color))
fig_type_cara_hem








set_colors40 <- c("#95ae74",
                  "#5668e6",
                  "#6cd84e",
                  "#5762be",
                  "#c0dc47",
                  "#7d89e9",
                  "#73a72f",
                  "#366bb3",
                  "#59e186",
                  "#6196e2",
                  "#b7dd79",
                  "#646ea3",
                  "#49a65a",
                  "#98a6d7",
                  "#557e2c",
                  "#3f78a6",
                  "#8fdf9f",
                  "#646c89",
                  "#57e5c4",
                  "#3a718b",
                  "#d1e2ab",
                  "#55b0db",
                  "#69773d",
                  "#b6cded",
                  "#4a7f51",
                  "#63d6e0",
                  "#737b55",
                  "#c6e9db",
                  "#3a7883",
                  "#8fd3bb",
                  "#667560",
                  "#a3bfc6",
                  "#3f7a61",
                  "#7794ab",
                  "#40a488",
                  "#52756e",
                  "#8ba691",
                  "#59a1a0",
                  "#4071f3",
                  "#42a359")

colors_125 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#33A02C",
                "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3",
                "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99",
                "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
                "#CCEBC5", "#FFED6F", "#4E79A7", "#F28E2B", "#E15759",
                "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                "#9C755F", "#BAB0AC", "#003F5C", "#58508D", "#BC5090",
                "#FF6361", "#FFA600", "#2F4B7C", "#665191", "#A05195",
                "#D45087", "#F95D6A", "#FF7C43", "#FFA600", "#264653",
                "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", "#0072B2",
                "#D55E00", "#CC79A7", "#F0E442", "#56B4E9", "#009E73",
                "#999999", "#00429D", "#73A2C6", "#ABDDA4", "#FEE08B",
                "#F46D43", "#A50026", "#5E4FA2", "#3288BD", "#66C2A5",
                "#ABD9E9", "#E6F598", "#FFFFBF", "#FEE090", "#FDAE61",
                "#F46D43", "#D53E4F", "#9E0142", "#440154", "#482777",
                "#3E4989", "#31688E", "#26828E", "#1F9E89", "#35B779",
                "#6CCE59", "#B4DE2C", "#FDE725", "#0D0887", "#46039F",
                "#7201A8", "#9C179E", "#BD3786", "#D8576B", "#ED7953",
                "#FB9F3A", "#FDCA26", "#F0F921", "#1A9850", "#91CF60",
                "#D9EF8B", "#FEE08B", "#FC8D59", "#D73027", "#7F3B08",
                "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#D8DAEB",
                "#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00",
                "#CA6702", "#BB3E03", "#AE2012", "#9B2226", "#3A0F73")


