library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(cowplot)
library(dplyr)
library(biovizBase)
library(hdf5r)
library(EnhancedVolcano)
set.seed(1234)

ctrl <- readRDS("~/hpchome/projects/CTRL_only/CTRL_only/Loadable Data/ctrl.rds")

##########Figure 2##########
FeaturePlot(ctrl, features = c("Egfr"), split.by = "orig.ident", pt.size = 0.5, reduction = "wnn.umap")
ggsave("Egfr_expression.tiff", width = 8, height = 8, units = c("in"), dpi = 300)

FeaturePlot(ctrl, features = c("Vim"), split.by = "orig.ident", pt.size = 0.5, reduction = "wnn.umap")
ggsave("Vim_expression.tiff", width = 8, height = 8, units = c("in"), dpi = 300)

#Subset clusters
epidermis <- subset(ctrl, idents = c('2', '5', '3', '10', '7', '8'))

dermis <- subset(ctrl, idents = c('18', '0', '1', '4', '6', '9'))

# Reorder clusters
epidermis.reorder <- epidermis
levels(x = epidermis.reorder) <- c('8', '7', '10', '3', '5', '2')

dermis.reorder <- dermis
levels(x = dermis.reorder) <- c('9', '6', '4', '1', '0', '18')

#Define markers
epidermis.markers <-c("Krt14", "Sox6", "Krt10", "Dsg1a", "Lor", "Hrnr", "Ptch2", "Edar")

dermis.markers <- c("Glis1", "Pde4d", "Twist2", "Crabp1", "Dlk1", "Agtr2")

#Dotplot of marker genes
dotplot_epidermis_markers <- DotPlot(epidermis.reorder, features = epidermis.markers, cols = c("lightgrey", "royalblue3"), dot.scale = 10,  col.min = -1.5, col.max = 2.5)
dotplot_epidermis_markers
ggsave("dotplot_epidermis_markers.tiff", width = 12, height = 5, units = c("in"), dpi = 300)

dotplot_dermis_markers <- DotPlot(dermis.reorder, features = dermis.markers, cols = c("lightgrey", "red2"), dot.scale = 10, col.min = -1.5, col.max = 2.5)
dotplot_dermis_markers
ggsave("dotplot_dermis_markers.tiff", width = 12, height = 5, units = c("in"), dpi = 300)

##########Figures 3-5##########

###Example Volcano
#Set Default Assay
DefaultAssay(ctrl) <- "ATAC"

#Find differentially expressed genes
deg.basalvspinous.ATAC <- FindMarkers(ctrl, ident.1 = c(2,5), ident.2 = c(3,10))
head(deg.basalvspinous.ATAC)
write.csv(deg.basalvspinous.ATAC, file = "deg.basalvspinous.ATAC.csv")

#Volcano Plot Visualization
volcano_dac_basalvspinous <- EnhancedVolcano(deg.basalvspinous.ATAC, 
                                             rownames(deg.basalvspinous.ATAC),
                                             x ="avg_log2FC", 
                                             y ="p_val", 
                                             FCcutoff = 0.1, 
                                             labSize = 5)
volcano_dac_basalvspinous
ggsave("volcano_dac_basalvspinous.tiff", width = 10, height = 10, units = c("in"), dpi = 300)

deg.basalvspinous.ATAC$diffexpressed <- "NO"
deg.basalvspinous.ATAC$diffexpressed[deg.basalvspinous.ATAC$avg_log2FC < -0.1 & deg.basalvspinous.ATAC$p_val < 0.05] <- "DOWN"
deg.basalvspinous.ATAC$diffexpressed[deg.basalvspinous.ATAC$avg_log2FC > 0.1 & deg.basalvspinous.ATAC$p_val < 0.05] <- "UP"
deg.basalvspinous.ATAC$delabel <- NA
deg.basalvspinous.ATAC <-tibble::rownames_to_column(deg.basalvspinous.ATAC, "Gene")
deg.basalvspinous.ATAC$delabel[deg.basalvspinous.ATAC$Gene %in% c("chrX-1000-1000", "chrX-1000-1000")] <- deg.basalvspinous.ATAC$Gene[deg.basalvspinous.ATAC$Gene %in% c("chrX-1000-1000", "chrX-1000-1000")] 

deg.basalvspinous.ATAC.volcano <- ggplot(deg.basalvspinous.ATAC, aes(x=avg_log2FC, y=-log10(p_val), col=diffexpressed, label=delabel)) + geom_point(alpha = 0.6, size = 2.5) + scale_color_manual(values=c("#d8952f", "gray80", "#449cd8")) + theme_classic() + geom_text_repel(aes(), size = 4, box.padding = unit(1, "lines"), point.padding = unit(0.3, "lines"), segment.linetype = 2, max.overlaps = 50, fontface = 2, colour = "#4b4d4e")+ geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed")
deg.basalvspinous.ATAC.volcano
ggsave("deg.basalvspinous.ATAC.volcano.tiff", deg.basalvspinous.ATAC.volcano, width = 8, height = 7, units = c("in"), dpi = 300)

###Example TF Enrichment
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)

DefaultAssay(ctrl) <- "ATAC"

main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- which(as.character(seqnames(granges(ctrl[["ATAC"]]))) %in% main.chroms)
ctrl[["ATAC"]] <- subset(ctrl[["ATAC"]], features = rownames(ctrl[["ATAC"]])[keep.peaks])
head(rownames(ctrl[["ATAC"]]))

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Add motif information
ctrl <- AddMotifs(
  object = ctrl,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

da_motif <- FindMarkers(
  object = ctrl,
  ident.1 = c(2,5),
  ident.2 = c(3,10),
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC'
)

sig_peak <- rownames(da_motif[da_motif$p_val_adj <.05,])

# find peaks open
open.peaks <- AccessiblePeaks(ctrl, idents = c(2, 5, 3, 10))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(ctrl, assay = "ATAC", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[sig_peak, ],
  n = 50000
)

da_motif <- FindMotifs(ctrl, features = sig_peak, background = peaks.matched)

da_motif.plot <- MotifPlot(
  object = ctrl,
  motifs = head(rownames(da_motif))
)
da_motif.plot

ggsave("basalvsspinous.motif.plot.tiff", da_motif.plot, height=10, width=10)

###Example View Peaks
combined_auto <- CoveragePlot(
  object = ctrl,
  region = "chrX-1000-1000",
  features = "gene",
  expression.assay = "SCT",
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 10000,
  extend.downstream = 10000,
  ymax = 16
)
combined_auto

color <- c("#00aefa", "#ff6c92","#00bdd1", "#f265e7", "#be80ff", "#00b7e8")
recolor <- combined_auto & scale_fill_manual(values = color)
recolor
ggsave("gene_promoter_combined_auto.tiff", recolor, height = 7, width = 10)
