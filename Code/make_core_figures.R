#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
package <- normalizePath(file.path(dirname(script_arg), ".."), mustWork=FALSE)
source_data <- file.path(package, "Source_Data")
out <- file.path(package, "Figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

blue <- "#0072B2"; orange <- "#D55E00"; green <- "#009E73"
purple <- "#CC79A7"; sky <- "#56B4E9"; gray <- "#777777"; light <- "#D9D9D9"
layer_cols <- c("L1_contrib"=blue, "L2_contrib"=orange, "L3_contrib"=green,
                "L4_contrib"=purple, "L5_contrib"=sky)

theme_safe <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.25),
        legend.key.size = unit(0.35, "cm"), plot.margin = margin(5, 5, 5, 5))

save_plot <- function(p, stem, w=7.2, h=3.2) {
  ggsave(file.path(out, paste0(stem, ".png")), p, width=w, height=h, dpi=300, bg="white")
  ggsave(file.path(out, paste0(stem, ".tiff")), p, width=w, height=h, dpi=300,
         compression="lzw", bg="white")
  ggsave(file.path(out, paste0(stem, ".pdf")), p, width=w, height=h, bg="white")
}

# Main Figure 4: descriptive proximity analysis
d <- read.csv(file.path(source_data, "Fig4_GWAS_3D_proximity.csv"))
d$direction <- ifelse(grepl("B_to_A|Opening", d$switch_type), "B to A", "A to B")
d$gene <- factor(d$gene, levels=d$gene[order(d$colocalization_score)])
p2a <- ggplot(d, aes(colocalization_score, gene, fill=direction)) +
  geom_col(width=.68) + scale_fill_manual(values=c("B to A"=orange, "A to B"=blue)) +
  labs(x="Proximity-derived score", y=NULL, fill=NULL, tag="A") + theme_safe
p2b <- ggplot(d, aes(gene, dist_to_switch_kb, fill=direction)) +
  geom_col(width=.68) + geom_hline(yintercept=500, linetype="dashed", color=orange) +
  scale_fill_manual(values=c("B to A"=orange, "A to B"=blue), guide="none") +
  labs(x=NULL, y="Distance to nearest switch (kb)", tag="B") +
  theme_safe + theme(axis.text.x=element_text(angle=50, hjust=1))
p2c <- ggplot(d, aes(chromatin_score, OB_expr, color=colocalization_score, label=gene)) +
  geom_point(size=2.5) + geom_text(nudge_x=.003, nudge_y=.04, hjust=0, size=2.2, check_overlap=TRUE) +
  scale_color_viridis_c(name="Proximity\nscore") +
  scale_x_continuous(expand=expansion(mult=c(.05,.28))) +
  labs(x="Compartment-change score", y="Osteoblast expression", tag="C") + theme_safe
fig2 <- (p2a | p2b | p2c)
save_plot(fig2, "Figure4_GWAS_3D_Proximity")

# Figure 3: four-candidate proteomic intersection, not a proteome-wide volcano
d <- read.csv(file.path(source_data, "SuppFig1_proteomic_candidate_intersection.csv"))
d$gene <- factor(d$gene, levels=d$gene[order(-log10(d$P.Value))])
d$comparison_label <- ifelse(grepl("AD_vs_BMMSC", d$comparison),
                             "AD and OB vs BMMSC", "OB vs AD and BMMSC")
d$switch_label <- ifelse(grepl("A_to_B", d$switch_type), "A→B", "B→A")
p3a <- ggplot(d, aes(-log10(P.Value), gene)) + geom_col(fill=purple, width=.68) +
  labs(x="Proteomic association, −log10(P)", y=NULL, tag="A") + theme_safe
long <- rbind(data.frame(gene=d$gene, component="3D chromatin", value=abs(d$chromatin_score)),
              data.frame(gene=d$gene, component="Proteomics (scaled)", value=-log10(d$P.Value)/10))
p3b <- ggplot(long, aes(value, gene, fill=component)) + geom_col(position="dodge", width=.65) +
  scale_fill_manual(values=c("3D chromatin"=blue, "Proteomics (scaled)"=purple)) +
  labs(x="Evidence component", y=NULL, fill=NULL, tag="B") + theme_safe
p3c <- ggplot(d, aes(comparison_label, gene, fill=switch_label)) + geom_tile(color="white") +
  scale_fill_manual(values=c("A→B"=blue, "B→A"=orange), name="Switch") +
  labs(x=NULL, y=NULL, tag="C") + theme_safe +
  theme(axis.text.x=element_text(angle=25, hjust=1, size=7), axis.line=element_blank(), axis.ticks=element_blank())
fig3 <- (p3a | p3b | p3c)
save_plot(fig3, "Supplementary_Figure1_3D_Proteomics_Candidates", h=3.5)

# Supplementary Figure 2: expression-product communication screen and data-derived heatmap
d <- read.csv(file.path(source_data, "SuppFig2_ligand_receptor_scores.csv"))
d$label <- paste0(d$ligand, "–", d$receptor, " (", substr(d$sender,1,5), "→", substr(d$receiver,1,5), ")")
top <- head(d[order(-d$interaction_score),], 12)
top$label <- factor(top$label, levels=rev(top$label))
ob <- d[d$sender == "Osteoblast",]
ob <- head(ob[order(-ob$interaction_score),], 10)
ob$label <- factor(ob$label, levels=rev(ob$label))
agg <- aggregate(interaction_score ~ sender + receiver, d, sum)
abbr <- c("Adipocyte"="Adip.", "Chondrocyte"="Chond.", "Endothelial"="Endo.",
          "Erythroid"="Ery.", "Fibroblast"="Fibro.", "Lymphocyte"="Lymph.",
          "Myeloid"="Myel.", "Osteoblast"="Osteo.", "Pericyte"="Peric.")
agg$receiver_short <- unname(abbr[agg$receiver])
p4a <- ggplot(top, aes(interaction_score, label)) + geom_col(fill=orange, width=.68) +
  scale_x_continuous(breaks=c(0,4,8)) +
  labs(x="Expression-product score", y=NULL, tag="A") + theme_safe
p4b <- ggplot(ob, aes(interaction_score, label)) + geom_col(fill=blue, width=.68) +
  scale_x_continuous(breaks=c(0,1,2)) +
  labs(x="Osteoblast-originating score", y=NULL, tag="B") + theme_safe
p4c <- ggplot(agg, aes(receiver_short, sender, fill=interaction_score)) + geom_tile(color="white", linewidth=.15) +
  scale_fill_viridis_c(name="Summed\nscore") +
  labs(x="Receiver cell class", y="Sender", tag="C") + theme_safe +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=7.5),
        axis.line=element_blank(), axis.ticks=element_blank(),
        plot.margin=margin(5, 8, 18, 8))
fig4 <- (p4a | p4b | p4c) + plot_layout(widths=c(1.25, 1.15, 1.65))
save_plot(fig4, "Supplementary_Figure2_Ligand_Receptor_Screen", w=8.8, h=4.5)

# Figure 5: coverage-aware ranking with explicit assay availability
d <- read.csv(file.path(source_data, "Fig5_3D_ORS_scores.csv"), check.names=FALSE)
d$available_score <- d[["3D_ORS_available"]]
d$strict_score <- d[["3D_ORS_strict_5layer"]]
eligible <- d[!is.na(d[["3D_ORS_available"]]),]
four <- eligible[eligible$evidence_layers_measured == 4,]
four <- head(four[order(-four[["3D_ORS_available"]]),], 10)
four$coverage_group <- "Four measured layers"
five <- eligible[eligible$evidence_layers_measured == 5,]
five <- head(five[order(-five[["3D_ORS_strict_5layer"]]),], 5)
five$coverage_group <- "All five layers"
top <- rbind(four, five)
top$coverage_group <- factor(top$coverage_group,
                             levels=c("Four measured layers", "All five layers"))
top$gene_display <- factor(top$gene, levels=rev(top$gene))
for (flag in c("L1_measured","L2_measured","L3_measured","L4_measured","L5_measured")) {
  top[[flag]] <- top[[flag]] %in% c(TRUE, "TRUE", "True")
}
measured_weight <- 0.30*top$L1_measured + 0.25*top$L2_measured +
  0.20*top$L3_measured + 0.15*top$L4_measured + 0.10*top$L5_measured
top$L1_contrib <- ifelse(top$L1_measured, 0.30*top$L1_3d_score/measured_weight, 0)
top$L2_contrib <- ifelse(top$L2_measured, 0.25*top$L2_gwas_score/measured_weight, 0)
top$L3_contrib <- ifelse(top$L3_measured, 0.20*top$L3_expr_score/measured_weight, 0)
top$L4_contrib <- ifelse(top$L4_measured, 0.15*top$L4_chip_score/measured_weight, 0)
top$L5_contrib <- ifelse(top$L5_measured, 0.10*top$L5_prot_score/measured_weight, 0)
long <- pivot_longer(top, cols=all_of(names(layer_cols)), names_to="layer", values_to="value")
long$layer <- factor(long$layer, levels=names(layer_cols))
p5a <- ggplot(long, aes(value, gene_display, fill=layer)) + geom_col(width=.72) +
  scale_fill_manual(values=layer_cols, labels=c("3D", "GWAS", "scRNA", "ChIP", "Proteomics")) +
  facet_grid(coverage_group ~ ., scales="free_y", space="free_y") +
  labs(x="Within-group coverage-normalized 3D-ORS", y=NULL, fill=NULL, tag="A") + theme_safe +
  theme(legend.position="bottom")
coverage <- as.data.frame(table(d$evidence_layers_measured))
names(coverage) <- c("layers", "n")
coverage$layers <- factor(coverage$layers, levels=c("3", "4", "5"))
p5b <- ggplot(coverage, aes(layers, n, fill=layers)) + geom_col(width=.62) +
  geom_text(aes(label=n), vjust=-.4, fontface="bold", size=3) +
  scale_fill_manual(values=c("3"=light, "4"=sky, "5"=green), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0, .12))) +
  labs(x="Measured evidence layers", y="Candidate genes", tag="B") + theme_safe
eligible$coverage_group <- factor(eligible$evidence_layers_measured, levels=c(4,5),
                                  labels=c("Four layers", "Five layers"))
p5c <- ggplot(eligible, aes(coverage_group, available_score, fill=coverage_group)) +
  geom_boxplot(width=.48, outlier.shape=NA, alpha=.7) +
  geom_jitter(width=.12, size=1.2, alpha=.65) +
  scale_fill_manual(values=c("Four layers"=sky, "Five layers"=green), guide="none") +
  labs(x=NULL, y="Coverage-normalized 3D-ORS", tag="C") + theme_safe
sens_raw <- read.csv(file.path(source_data, "Supplementary_3D_ORS_leave_one_layer_out_gene_ranks.csv"))
sens <- aggregate(leave_one_layer_out_rank ~ gene, sens_raw,
                  function(x) c(mean=mean(x), min=min(x), max=max(x)))
sens <- data.frame(gene=sens$gene, mean_rank=sens$leave_one_layer_out_rank[,"mean"],
                   min_rank=sens$leave_one_layer_out_rank[,"min"],
                   max_rank=sens$leave_one_layer_out_rank[,"max"])
sens <- sens[order(sens$mean_rank),]
sens$gene <- factor(sens$gene, levels=rev(sens$gene))
p5d <- ggplot(sens, aes(mean_rank, gene)) +
  geom_errorbar(aes(xmin=min_rank, xmax=max_rank), width=.2, color=gray, orientation="y") +
  geom_point(color=blue, size=2) + scale_x_reverse() +
  labs(x="Rank after omitting one layer (lower is better)", y=NULL, tag="D") + theme_safe
fig5 <- (p5a | p5b) / (p5c | p5d)
save_plot(fig5, "Figure5_3DORS_Coverage_Aware", h=5.8)

message("Wrote revised figures to: ", out)
