## Community analysis of in situ lobster embryo, pleopod, carapace, and hemolymph
## microbiomes (16S rRNA)
##
## Input:  data/phyob.rds (made by build_phyloseq.R)
## Output: results/ (tables) and results/figures/ (plots)
##
## Sections: 1 all samples | 2 embryos | 3 pleopods | 4 embryos vs. pleopods |
##           5 carapace | 6 hemolymph | 7 carapace vs. hemolymph

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggpubr)
library(ANCOMBC)
library(VennDiagram)
library(grid)

source("palettes.R")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

phyob <- readRDS("data/phyob.rds")

sample_order <- c("Early_Embryo", "Hemolymph", "Late_Embryo", "Pleo_Swab", "Shell", "Cara_Swab")
phyob@sam_data$Sample_Type <- factor(phyob@sam_data$Sample_Type, levels = sample_order)

# Collected test results, written to results/ at the end
kw_results       <- list()
perm_results     <- list()
rarefaction_log  <- list()


## ======================= Helper functions =======================

save_fig <- function(p, name, width = 7, height = 5) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(file.path("results/figures", paste0(name, ".png")), p,
         width = width, height = height, dpi = 300)
}

# Keep samples by metadata value(s), then drop taxa absent from the remaining samples
select_samples <- function(ps, column, include = NULL, exclude = NULL) {
  v <- as.character(data.frame(sample_data(ps))[[column]])
  keep <- rep(TRUE, length(v))
  if (!is.null(include)) keep <- keep & v %in% include
  if (!is.null(exclude)) keep <- keep & !(v %in% exclude)
  ps <- prune_samples(sample_names(ps)[keep], ps)
  prune_taxa(taxa_sums(ps) > 0, ps)
}

# Shannon and observed richness (non-rarefied) joined to the sample metadata
alpha_table <- function(ps) {
  rich <- estimate_richness(ps, split = TRUE, measures = c("Observed", "Shannon"))
  rich$Sample <- sample_names(ps)
  left_join(rich, as.data.frame(sample_data(ps)), by = "Sample")
}

kw_test <- function(df, group, label) {
  t <- kruskal.test(as.formula(paste("Shannon ~", group)), data = df)
  data.frame(analysis = label, test = paste("Kruskal-Wallis: Shannon ~", group),
             statistic = unname(t$statistic), df = unname(t$parameter), p = t$p.value)
}

# PERMANOVA on Bray-Curtis distances (non-rarefied counts), plus a dispersion test
permanova <- function(ps, formula, label, disp_group = NULL, seed = 42) {
  set.seed(seed)
  meta <- droplevels(data.frame(sample_data(ps)))
  dist_bray <- phyloseq::distance(ps, method = "bray")
  fit <- adonis2(as.formula(paste("dist_bray ~", formula)), data = meta)
  out <- data.frame(analysis = label, term = rownames(fit), df = fit$Df,
                    R2 = fit$R2, F = fit$F, p = fit$`Pr(>F)`, row.names = NULL)
  out <- out[!out$term %in% c("Residual", "Total"), ]
  if (!is.null(disp_group)) {
    disp <- betadisper(dist_bray, meta[[disp_group]])
    out$dispersion_group <- disp_group
    out$dispersion_p <- anova(disp)$`Pr(>F)`[1]
  }
  out
}

pcoa_theme <- function(legend_pos) {
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"),
        legend.key.size = unit(0.6, "cm"), legend.text = element_text(size = 13),
        legend.position = legend_pos, legend.background = element_rect(fill = "grey97"),
        legend.title = element_blank(),
        panel.grid = element_line(color = "grey85", linetype = 4),
        panel.background = element_rect(fill = "grey95"))
}

# Rarefy to the smallest library, run PCoA on Bray-Curtis, save the plot.
# Returns the rarefaction depth so it can be logged.
run_pcoa <- function(ps, name, color, shape = color, colors, shapes, seed,
                     ellipse_by = color, legend_pos = c(0.85, 0.2)) {
  set.seed(seed)
  depth <- min(sample_sums(ps))
  rare <- rarefy_even_depth(ps, sample.size = depth, replace = TRUE, verbose = FALSE)
  ord <- ordinate(rare, method = "PCoA", distance = "bray")
  p <- plot_ordination(rare, ord, color = color, shape = shape) +
    geom_point(aes(color = .data[[color]]), alpha = 1, size = 4)
  if (!is.null(ellipse_by)) {
    p <- p + stat_ellipse(aes(group = .data[[ellipse_by]]), type = "norm", level = 0.95)
  }
  p <- p + scale_color_manual(values = colors, name = "") +
    scale_shape_manual(values = shapes, name = "") + pcoa_theme(legend_pos)
  save_fig(p, name)
  data.frame(figure = name, rarefied_depth = depth, n_samples = nsamples(ps))
}

# Stacked relative-abundance barplot at a taxonomic rank
rel_abund_barplot <- function(ps, rank, min_abund, colors, facet_levels = NULL) {
  d <- ps %>%
    tax_glom(taxrank = rank) %>%
    transform_sample_counts(function(x) x / sum(x) * 100) %>%
    psmelt() %>%
    dplyr::filter(Abundance > min_abund) %>%
    dplyr::arrange(.data[[rank]])
  p <- ggplot(d, aes(x = Sample, y = Abundance, fill = .data[[rank]])) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = colors) +
    guides(x = guide_axis(angle = 90),
           fill = guide_legend(keywidth = 2, keyheight = 1, ncol = 1)) +
    labs(fill = rank, x = "Sample ID",
         y = paste0("Relative ", rank, " Abundance (>", min_abund, "%)\n")) +
    theme(axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 12), axis.line = element_line(color = "black"),
          axis.text.x = element_text(size = 10), strip.text = element_text(size = 11))
  if (!is.null(facet_levels)) {
    p <- p + facet_grid(~factor(Type, levels = facet_levels), scales = "free", space = "free")
  }
  p
}

# Genera / ASVs present in at least `min_samples` samples (orientation-safe)
presence_set <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  otu
}

asv_prevalence <- function(ps, min_samples = 2) {
  otu <- presence_set(ps)
  rownames(otu)[rowSums(otu > 0) >= min_samples]
}

genus_prevalence <- function(ps, min_samples = 2) {
  ps <- subset_taxa(ps, !is.na(Genus) & Genus != "")
  ps <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)
  otu <- presence_set(ps)
  genus <- as.character(tax_table(ps)[, "Genus"])
  unique(genus[rowSums(otu > 0) >= min_samples])
}

# Two-set Venn diagram, saved as PNG, with the shared set written to a text file
venn2 <- function(a, b, labels, fills, name) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  v <- draw.pairwise.venn(area1 = length(a), area2 = length(b),
                          cross.area = length(intersect(a, b)),
                          category = labels, fill = fills, alpha = 0.6,
                          cex = 1.5, cat.cex = 1.2, ind = FALSE)
  png(file.path("results/figures", paste0(name, ".png")), width = 1600, height = 1600, res = 300)
  grid.newpage()
  grid.draw(v)
  dev.off()
  writeLines(sort(intersect(a, b)), file.path("results", paste0(name, "_shared.txt")))
  data.frame(comparison = name, set1 = labels[1], n1 = length(a),
             set2 = labels[2], n2 = length(b), shared = length(intersect(a, b)))
}
venn_log <- list()

# ANCOM-BC2 at genus level. `pseudo` is the pseudo-count (R coerced TRUE/FALSE
# to 1/0 in the original calls, so those are written out explicitly here).
run_ancombc2 <- function(ps, fix_formula, group, pseudo = 0, struc_zero = FALSE, seed = 123) {
  set.seed(seed)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  sample_data(ps) <- sample_data(droplevels(data.frame(sample_data(ps))))
  ancombc2(ps, taxa_are_rows = FALSE, tax_level = "Genus", fix_formula = fix_formula,
           group = group, pseudo = pseudo, pseudo_sens = TRUE, prv_cut = 0.1, alpha = 0.05,
           struc_zero = struc_zero, global = FALSE, pairwise = FALSE, dunnet = FALSE,
           trend = TRUE, verbose = FALSE)$res
}

save_ancom <- function(res, name) {
  write.csv(res, file.path("results", paste0("ancombc2_", name, ".csv")), row.names = FALSE)
}

# Bar plot of log fold changes for genera with a significant `term`.
# Axis labels are colored by the ANCOM-BC2 sensitivity-analysis result if robust_labels = TRUE.
lfc_barplot <- function(res, term, pos_label, neg_label, fills, title,
                        min_abs_lfc = 0, robust_labels = FALSE) {
  lfc <- paste0("lfc_", term); se <- paste0("se_", term)
  diff <- paste0("diff_", term); robust <- paste0("diff_robust_", term)
  d <- res %>%
    dplyr::filter(.data[[diff]] %in% TRUE, abs(.data[[lfc]]) > min_abs_lfc) %>%
    dplyr::arrange(dplyr::desc(.data[[lfc]]))
  if (nrow(d) == 0) {
    message("No significant genera for: ", title)
    return(NULL)
  }
  label_colors <- if (robust_labels && robust %in% names(d)) {
    ifelse(d[[robust]] %in% TRUE, "aquamarine3", "black")
  } else "black"
  d <- d %>% dplyr::mutate(
    direct = factor(ifelse(.data[[lfc]] > 0, pos_label, neg_label), levels = c(pos_label, neg_label)),
    taxon = factor(taxon, levels = taxon))
  ggplot(d, aes(x = taxon, y = .data[[lfc]], fill = direct)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = .data[[lfc]] - .data[[se]], ymax = .data[[lfc]] + .data[[se]]),
                  width = 0.2, position = position_dodge(0.05), color = "black") +
    labs(x = NULL, y = "Log fold change", title = title) +
    scale_fill_manual(values = setNames(fills, c(pos_label, neg_label)), name = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1, color = label_colors))
}


## ======================= 1. All samples =======================

# Alpha diversity (non-rarefied)
richness_plot <- plot_richness(phyob, x = "Sample_Type", color = "Sample_Type",
                               measures = c("Observed", "Shannon"),
                               title = "Alpha Diversity (non-rarefied)", nrow = 2, scales = "free_y") +
  geom_point(aes(color = Sample_Type), alpha = 1, size = 7, shape = 20) +
  stat_summary(fun.data = mean_se, geom = "errorbar", colour = "gray10", linewidth = 0.9, width = 0.25) +
  scale_color_manual(values = pal_sample_type, name = "") +
  scale_x_discrete(labels = c(Early_Embryo = "Early Embryo", Late_Embryo = "Late Embryo",
                              Hemolymph = "Hemolymph", Cara_Swab = "Carapace Swab",
                              Shell = "Shell", Pleo_Swab = "Pleopod Swab")) +
  theme(axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 17),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 17, angle = 35, hjust = 0.9, vjust = 0.88),
        strip.text = element_text(size = 17), axis.ticks.length = unit(0.3, "cm"),
        plot.title = element_text(size = 19), legend.position = "none",
        panel.grid = element_line(color = "grey87", linetype = 4),
        panel.background = element_rect(fill = "grey95"))
save_fig(richness_plot, "alpha_diversity_by_sample_type", width = 8, height = 9)

alpha_all <- alpha_table(phyob)

shannon_box <- ggboxplot(alpha_all, x = "Sample_Type", y = "Shannon", color = "Sample_Type",
                         palette = unname(pal_sample_type[sample_order]), order = sample_order,
                         ylab = "Alpha Diversity", xlab = "Sample Type")
save_fig(shannon_box, "shannon_boxplot_by_sample_type", width = 8, height = 5)

kw_results$all_samples <- kw_test(alpha_all, "Sample_Type", "All samples")
pw <- pairwise.wilcox.test(alpha_all$Shannon, alpha_all$Sample_Type, p.adjust.method = "BH")
pw_table <- as.data.frame(as.table(pw$p.value)) %>%
  dplyr::filter(!is.na(Freq)) %>%
  dplyr::rename(group1 = Var1, group2 = Var2, p_adj_BH = Freq)
write.csv(pw_table, "results/alpha_pairwise_wilcox_sample_type.csv", row.names = FALSE)

# Mean alpha diversity by sample type and by sample kind
write.csv(alpha_all %>% group_by(Sample_Type) %>%
            summarise(n = n(), mean_shannon = mean(Shannon), mean_observed = mean(Observed)),
          "results/alpha_means_by_sample_type.csv", row.names = FALSE)
write.csv(alpha_all %>% group_by(Sample_Kind) %>%
            summarise(n = n(), mean_shannon = mean(Shannon), mean_observed = mean(Observed)),
          "results/alpha_means_by_sample_kind.csv", row.names = FALSE)

# Beta diversity
rarefaction_log$all_samples <- run_pcoa(phyob, "pcoa_all_samples", color = "Type",
                                        colors = pal_type, shapes = shape_type, seed = 3)
perm_results$all_samples <- permanova(phyob, "Sample_Type", "All samples",
                                      disp_group = "Sample_Type")

meta_all <- droplevels(data.frame(sample_data(phyob)))
disp_all <- betadisper(phyloseq::distance(phyob, method = "bray"), meta_all$Sample_Type)
write.csv(TukeyHSD(disp_all)$group, "results/dispersion_tukey_sample_type.csv")

# Phylum composition
phylum_bp <- rel_abund_barplot(phyob, "Phylum", 1, set_colors40,
                               facet_levels = c("Early Embryo", "Late Embryo", "Pleopod",
                                                "Carapace", "Hemolymph"))
save_fig(phylum_bp, "phylum_barplot_all_samples", width = 14, height = 7)


## ======================= 2. Embryos =======================

ps_early_emb <- select_samples(phyob, "Sample_Type", "Early_Embryo")
ps_late_emb  <- select_samples(phyob, "Sample_Type", "Late_Embryo")
ps_emb       <- select_samples(phyob, "Sample_Kind", "Embryo")
late_emb_nosme <- select_samples(ps_late_emb, "Location", exclude = "SME")
emb_nosme      <- select_samples(ps_emb, "Location", exclude = "SME")

# Beta diversity by location (the late-embryo PCoA includes SME; its PERMANOVA does not)
rarefaction_log$early_emb <- run_pcoa(ps_early_emb, "pcoa_early_embryos_location", color = "Location",
                                      colors = pal_location, shapes = shape_location, seed = 4)
rarefaction_log$late_emb  <- run_pcoa(ps_late_emb, "pcoa_late_embryos_location", color = "Location",
                                      colors = pal_location, shapes = shape_location, seed = 5)
perm_results$early_emb_location <- permanova(ps_early_emb, "Location", "Early embryos ~ Location",
                                             disp_group = "Location")
perm_results$late_emb_location  <- permanova(late_emb_nosme, "Location", "Late embryos (no SME) ~ Location",
                                             disp_group = "Location")

# Early vs. late embryos
rarefaction_log$emb_stage <- run_pcoa(ps_emb, "pcoa_embryos_stage", color = "Sample_Type",
                                      colors = pal_sample_type, shapes = shape_sample_type,
                                      seed = 6, legend_pos = c(0.25, 0.2))
perm_results$emb_stage <- permanova(ps_emb, "Embryogenesis", "All embryos ~ Embryogenesis",
                                    disp_group = "Embryogenesis")

# Alpha diversity
kw_results$early_emb_location <- kw_test(alpha_table(ps_early_emb), "Location", "Early embryos ~ Location")
kw_results$late_emb_location  <- kw_test(alpha_table(late_emb_nosme), "Location", "Late embryos (no SME) ~ Location")
kw_results$emb_stage          <- kw_test(alpha_table(ps_emb), "Sample_Type", "All embryos: early vs. late")

# Shared genera and ASVs, early vs. late embryos (present in >= 2 samples)
emb_early <- select_samples(ps_emb, "Embryogenesis", "Early")
emb_late  <- select_samples(ps_emb, "Embryogenesis", "Late")
venn_log$emb_genera <- venn2(genus_prevalence(emb_early), genus_prevalence(emb_late),
                             c("Early Embryos", "Late Embryos"), c("#519c3d", "tomato"),
                             "venn_genera_embryos_early_late")
venn_log$emb_asvs <- venn2(asv_prevalence(ps_early_emb), asv_prevalence(ps_late_emb),
                           c("Early Embryos", "Late Embryos"), c("#519c3d", "tomato"),
                           "venn_asvs_embryos_early_late")

# ANCOM-BC2, embryos without SME. Three models:
#   stage only | stage + location (for the location effect) | stage * location (interaction)
res_emb_stage <- run_ancombc2(emb_nosme, "Embryogenesis", group = "Embryogenesis", struc_zero = TRUE)
save_ancom(res_emb_stage, "embryos_stage")
save_fig(lfc_barplot(res_emb_stage, "EmbryogenesisLate", "Enriched on Late", "Enriched on Early",
                     c("#779fe5", "#519c3d"),
                     "Genera log fold change on late vs. early-stage embryos"),
         "ancombc2_embryos_stage", width = 10, height = 6)

res_emb_add <- run_ancombc2(emb_nosme, "Embryogenesis + Location", group = "Embryogenesis", struc_zero = TRUE)
save_ancom(res_emb_add, "embryos_stage_plus_location")
save_fig(lfc_barplot(res_emb_add, "LocationNH", "Enriched in NH", "Enriched in ME",
                     c(pal_location[["NH"]], pal_location[["DownEast"]]),
                     "Genera log fold change in NH vs. ME in all embryos", robust_labels = TRUE),
         "ancombc2_embryos_location", width = 10, height = 6)

res_emb_int <- run_ancombc2(emb_nosme, "Embryogenesis * Location", group = "Location",
                            pseudo = 1, struc_zero = TRUE)
res_emb_int <- res_emb_int %>%
  dplyr::mutate(lfc_NH_vs_Maine_Late = lfc_LocationNH + .data[["lfc_EmbryogenesisLate:LocationNH"]])
save_ancom(res_emb_int, "embryos_stage_x_location")
save_fig(lfc_barplot(res_emb_int, "EmbryogenesisLate", "Enriched on Late", "Enriched on Early",
                     c("#779fe5", "#519c3d"),
                     "Genera log fold change on late vs. early-stage embryos in Maine",
                     robust_labels = TRUE),
         "ancombc2_embryos_stage_in_maine", width = 10, height = 6)
save_fig(lfc_barplot(res_emb_int, "LocationNH", "Enriched in NH", "Enriched in ME",
                     c(pal_location[["NH"]], pal_location[["DownEast"]]),
                     "Genera log fold change in Maine vs. NH in early-stage embryos",
                     robust_labels = TRUE),
         "ancombc2_embryos_location_in_early", width = 10, height = 6)

# Genera whose early-to-late change depends on location
write.csv(res_emb_int %>%
            dplyr::filter(.data[["diff_EmbryogenesisLate:LocationNH"]] %in% TRUE) %>%
            dplyr::select(taxon, `lfc_EmbryogenesisLate:LocationNH`, `q_EmbryogenesisLate:LocationNH`,
                          lfc_NH_vs_Maine_Late),
          "results/embryos_location_dependent_stage_effect.csv", row.names = FALSE)


## ======================= 3. Pleopods =======================

phyob_pleo   <- select_samples(phyob, "Sample_Type", "Pleo_Swab")
pleo_early   <- select_samples(phyob_pleo, "Embryogenesis", "Early")
pleo_late    <- select_samples(phyob_pleo, "Embryogenesis", "Late")
pleo_late_nosme <- select_samples(pleo_late, "Location", exclude = "SME")

# Alpha diversity
alpha_pleo <- alpha_table(phyob_pleo)
kw_results$pleo_stage          <- kw_test(alpha_pleo, "Embryogenesis", "Pleopods: early vs. late")
kw_results$pleo_early_location <- kw_test(alpha_table(pleo_early), "Location", "Early pleopods ~ Location")
kw_results$pleo_late_location  <- kw_test(alpha_table(pleo_late), "Location", "Late pleopods ~ Location")
save_fig(ggboxplot(alpha_pleo, x = "Embryogenesis", y = "Shannon", color = "Embryogenesis",
                   palette = unname(pal_stage[c("Early", "Late")]), order = c("Early", "Late"),
                   ylab = "Alpha Diversity", xlab = "Embryogenesis stage"),
         "shannon_boxplot_pleopods_stage")

# Beta diversity
rarefaction_log$pleo <- run_pcoa(phyob_pleo, "pcoa_pleopods_stage", color = "Embryogenesis",
                                 shape = "Location", colors = pal_stage, shapes = shape_location,
                                 seed = 446, legend_pos = c(0.15, 0.8))
rarefaction_log$pleo_early <- run_pcoa(pleo_early, "pcoa_early_pleopods_location", color = "Location",
                                       colors = pal_location, shapes = shape_location,
                                       seed = 449, legend_pos = c(0.15, 0.8))
rarefaction_log$pleo_late <- run_pcoa(pleo_late, "pcoa_late_pleopods_location", color = "Location",
                                      colors = pal_location, shapes = shape_location,
                                      seed = 449, legend_pos = c(0.15, 0.8))
perm_results$pleo_stage <- permanova(phyob_pleo, "Embryogenesis", "Pleopods ~ Embryogenesis",
                                     disp_group = "Embryogenesis", seed = 47)
perm_results$pleo_early <- permanova(pleo_early, "Location", "Early pleopods ~ Location",
                                     disp_group = "Location", seed = 48)
perm_results$pleo_late  <- permanova(pleo_late_nosme, "Location", "Late pleopods (no SME) ~ Location",
                                     disp_group = "Location", seed = 49)

# Shared genera, early vs. late pleopods
venn_log$pleo_genera <- venn2(genus_prevalence(pleo_early), genus_prevalence(pleo_late),
                              c("Early Pleopods", "Late Pleopods"), c("steelblue", "tomato"),
                              "venn_genera_pleopods_early_late")

# ANCOM-BC2, early vs. late pleopods adjusting for location
res_pleo <- run_ancombc2(phyob_pleo, "Embryogenesis + Location", group = "Embryogenesis")
save_ancom(res_pleo, "pleopods_stage_plus_location")
save_fig(lfc_barplot(res_pleo, "EmbryogenesisLate", "Enriched on Late", "Enriched on Early",
                     c("#779fe5", "#519c3d"),
                     "Genera log fold change with > |1.5| on late-stage pleopods",
                     min_abs_lfc = 1.5, robust_labels = TRUE),
         "ancombc2_pleopods_stage", width = 10, height = 6)


## ======================= 4. Embryos vs. pleopods =======================

ps_emb_pleo <- select_samples(phyob, "Sample_Kind", c("Pleopod", "Embryo"))

rarefaction_log$emb_pleo <- run_pcoa(ps_emb_pleo, "pcoa_embryos_vs_pleopods", color = "Sample_Type",
                                     shape = "Embryogenesis", colors = pal_sample_type,
                                     shapes = shape_stage, seed = 543, legend_pos = c(0.9, 0.15))
perm_results$emb_pleo <- permanova(ps_emb_pleo, "Sample_Type", "Embryos + pleopods ~ Sample_Type",
                                   disp_group = "Sample_Type", seed = 57)

# Shared genera, all embryos vs. all pleopods
venn_log$emb_pleo_genera <- venn2(genus_prevalence(ps_emb), genus_prevalence(phyob_pleo),
                                  c("Embryos", "Pleopods"), c("steelblue", "tomato"),
                                  "venn_genera_embryos_pleopods")

pe_fills <- c("#3f78a6", "#519c3d")

# Early stage: pleopods vs. embryos, plus a location interaction model
early_pe <- select_samples(ps_emb_pleo, "Embryogenesis", "Early")
perm_results$early_pe <- permanova(early_pe, "Sample_Type", "Early embryos + pleopods ~ Sample_Type",
                                   disp_group = "Sample_Type", seed = 90)
res_early_pe <- run_ancombc2(early_pe, "Sample_Kind", group = "Sample_Kind", pseudo = 1, struc_zero = TRUE)
save_ancom(res_early_pe, "early_pleopods_vs_embryos")
fig_early_pe <- lfc_barplot(res_early_pe, "Sample_KindPleopod", "Enriched on Pleopods",
                            "Enriched on Embryos", pe_fills, "Early pleopods vs embryos")

res_early_pe_int <- run_ancombc2(early_pe, "Sample_Kind * Location", group = "Sample_Kind",
                                 pseudo = 1, struc_zero = TRUE)
save_ancom(res_early_pe_int, "early_pleopods_vs_embryos_x_location")
int_df <- res_early_pe_int %>%
  dplyr::filter(.data[["q_Sample_KindPleopod:LocationNH"]] < 0.05) %>%
  dplyr::transmute(taxon, Maine = lfc_Sample_KindPleopod,
                   `New Hampshire` = lfc_Sample_KindPleopod + .data[["lfc_Sample_KindPleopod:LocationNH"]]) %>%
  tidyr::pivot_longer(c(Maine, `New Hampshire`), names_to = "Location", values_to = "LFC")
if (nrow(int_df) > 0) {
  fig_int <- ggplot(int_df, aes(x = Location, y = LFC, group = taxon, color = taxon)) +
    geom_line() + geom_point(size = 2) + geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Location", y = "Log fold change (Pleopod vs Embryo)",
         title = "Location-dependent early pleopod-embryo effects across genera") +
    theme_classic() +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16),
          legend.title = element_text(size = 16), legend.text = element_text(size = 14))
  save_fig(fig_int, "ancombc2_early_pleopods_vs_embryos_interaction", width = 10, height = 6)
}

# Late stage (excluding SME and genera labeled "Unknown")
late_pe <- select_samples(ps_emb_pleo, "Embryogenesis", "Late")
late_pe <- select_samples(late_pe, "Location", exclude = "SME")
perm_results$late_pe <- permanova(late_pe, "Sample_Type", "Late embryos + pleopods ~ Sample_Type",
                                  disp_group = "Sample_Type", seed = 93)
late_pe_known <- subset_taxa(late_pe, !is.na(Genus) & Genus != "Unknown")
res_late_pe <- run_ancombc2(late_pe_known, "Sample_Kind", group = "Sample_Kind", pseudo = 1, struc_zero = TRUE)
save_ancom(res_late_pe, "late_pleopods_vs_embryos")
fig_late_pe <- lfc_barplot(res_late_pe, "Sample_KindPleopod", "Enriched on Pleopods",
                           "Enriched on Embryos", pe_fills, "Late pleopods vs embryos")

if (!is.null(fig_early_pe) && !is.null(fig_late_pe)) {
  save_fig(ggarrange(fig_early_pe, fig_late_pe, ncol = 2, labels = c("A", "B"),
                     common.legend = TRUE, legend = "bottom"),
           "ancombc2_pleopods_vs_embryos_early_late", width = 14, height = 6)
}


## ======================= 5. Carapace (Shell and Cara_Swab samples) =======================

phyob_shell <- select_samples(phyob, "Sample_Type", c("Shell", "Cara_Swab"))
shell_nosme <- select_samples(phyob_shell, "Location", exclude = "SME")
shell_early <- select_samples(phyob_shell, "Embryogenesis", "Early")
shell_late  <- select_samples(phyob_shell, "Embryogenesis", "Late")
# One SME sample is removed because tests need more than one sample per group
shell_late2 <- select_samples(shell_late, "Sample", exclude = "L_ME1_cswab")

# Alpha diversity
alpha_shell <- alpha_table(phyob_shell)
kw_results$shell_stage    <- kw_test(alpha_shell, "Embryogenesis", "Carapace: early vs. late")
kw_results$shell_location <- kw_test(alpha_shell, "Location", "Carapace ~ Location")
kw_results$shell_early_location <- tryCatch(
  kw_test(alpha_table(shell_early), "Location", "Early carapace ~ Location"),
  error = function(e) { message("Early carapace location test skipped: ", conditionMessage(e)); NULL })
kw_results$shell_late_location <- kw_test(alpha_table(shell_late), "Location", "Late carapace ~ Location")

# Beta diversity, all carapace samples (no SME)
rarefaction_log$shell <- run_pcoa(shell_nosme, "pcoa_carapace_stage", color = "Embryogenesis",
                                  shape = "Location", colors = pal_stage, shapes = shape_location, seed = 43)
perm_results$shell_stage <- permanova(shell_nosme, "Embryogenesis", "Carapace (no SME) ~ Embryogenesis",
                                      disp_group = "Embryogenesis", seed = 42)

# Late-stage carapace: ESD status x location
rarefaction_log$shell_late <- run_pcoa(shell_late2, "pcoa_late_carapace_esd", color = "Location",
                                       shape = "ESD_Status", colors = pal_location, shapes = shape_esd,
                                       seed = 44, ellipse_by = NULL)
perm_results$shell_late <- permanova(shell_late2, "ESD_Status * Location",
                                     "Late carapace ~ ESD_Status * Location",
                                     disp_group = "Location", seed = 45)

# Because the interaction was significant, test ESD status within each location
shell_lateNH <- select_samples(shell_late2, "Location", "NH")
shell_lateME <- select_samples(shell_late2, "Location", "DownEast")
perm_results$shell_late_NH <- permanova(shell_lateNH, "ESD_Status", "Late carapace, NH ~ ESD_Status",
                                        disp_group = "ESD_Status", seed = 87)
perm_results$shell_late_ME <- permanova(shell_lateME, "ESD_Status", "Late carapace, ME ~ ESD_Status",
                                        disp_group = "ESD_Status", seed = 87)
rarefaction_log$shell_late_NH <- run_pcoa(shell_lateNH, "pcoa_late_carapace_NH", color = "Location",
                                          shape = "ESD_Status", colors = pal_location, shapes = shape_esd,
                                          seed = 77, ellipse_by = "ESD_Status", legend_pos = c(0.2, 0.85))
rarefaction_log$shell_late_ME <- run_pcoa(shell_lateME, "pcoa_late_carapace_ME", color = "Location",
                                          shape = "ESD_Status", colors = pal_location, shapes = shape_esd,
                                          seed = 78, ellipse_by = "ESD_Status", legend_pos = c(0.2, 0.85))

# ANCOM-BC2: stage (adjusting for location and ESD status)
res_shell <- run_ancombc2(shell_nosme, "Embryogenesis + Location + ESD_Status", group = "Embryogenesis")
save_ancom(res_shell, "carapace_stage_location_esd")
save_fig(lfc_barplot(res_shell, "EmbryogenesisLate", "Positive LFC", "Negative LFC",
                     c("#779fe5", "#519c3d"),
                     "Genera log fold change on late-stage carapace controlled for location & ESD status",
                     robust_labels = TRUE),
         "ancombc2_carapace_stage", width = 10, height = 6)

# ANCOM-BC2: ESD status in early and late carapace
res_shell_early_esd <- run_ancombc2(shell_early, "ESD_Status + Location", group = "ESD_Status", pseudo = 1)
save_ancom(res_shell_early_esd, "early_carapace_esd")
res_shell_late_esd <- run_ancombc2(shell_late2, "ESD_Status * Location", group = "ESD_Status", pseudo = 1)
save_ancom(res_shell_late_esd, "late_carapace_esd_x_location")


## ======================= 6. Hemolymph =======================

phyob_hem  <- select_samples(phyob, "Sample_Type", "Hemolymph")
hem_nosme  <- select_samples(phyob_hem, "Location", exclude = "SME")

kw_results$hem_location <- kw_test(alpha_table(hem_nosme), "Location", "Hemolymph (no SME) ~ Location")
rarefaction_log$hem <- run_pcoa(hem_nosme, "pcoa_hemolymph_location", color = "Location",
                                colors = pal_location, shapes = shape_location, seed = 47)
perm_results$hem_esd <- permanova(hem_nosme, "ESD_Status", "Hemolymph (no SME) ~ ESD_Status",
                                  disp_group = "ESD_Status", seed = 49)

save_fig(rel_abund_barplot(phyob_hem, "Genus", 5, colors_125), "genus_barplot_hemolymph",
         width = 14, height = 7)


## ======================= 7. Carapace vs. hemolymph =======================

ps_shell_hem <- select_samples(phyob, "Sample_Type", c("Shell", "Cara_Swab", "Hemolymph"))

rarefaction_log$shell_hem <- run_pcoa(ps_shell_hem, "pcoa_carapace_vs_hemolymph", color = "Type",
                                      colors = pal_type, shapes = shape_type, seed = 543,
                                      legend_pos = c(0.9, 0.15))

# Shared genera (any presence)
venn_log$shell_hem_genera <- venn2(genus_prevalence(phyob_shell, 1), genus_prevalence(phyob_hem, 1),
                                   c("Carapace", "Hemolymph"), c("#c17fee", "#3f78a6"),
                                   "venn_genera_carapace_hemolymph")

res_shell_hem <- run_ancombc2(ps_shell_hem, "Type", group = "Type", pseudo = 1, struc_zero = TRUE)
save_ancom(res_shell_hem, "carapace_vs_hemolymph")
save_fig(lfc_barplot(res_shell_hem, "TypeHemolymph", "Enriched on Hemolymph", "Enriched on Carapace",
                     c("#3f78a6", "#519c3d"), "Hemolymph vs Carapace"),
         "ancombc2_carapace_vs_hemolymph", width = 10, height = 6)


## ======================= Write collected results =======================

write.csv(dplyr::bind_rows(kw_results),   "results/alpha_diversity_tests.csv", row.names = FALSE)
write.csv(dplyr::bind_rows(perm_results), "results/permanova_tests.csv", row.names = FALSE)
write.csv(dplyr::bind_rows(rarefaction_log), "results/rarefaction_depths.csv", row.names = FALSE)
write.csv(dplyr::bind_rows(venn_log),     "results/venn_set_sizes.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
