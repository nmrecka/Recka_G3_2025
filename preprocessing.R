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

PATH_ATAC_PEAKS = "~/hpchome/projects/2023_10X_multiomics_E15.5_murine_skin/raw_data/dnacore454.healthcare.uiowa.edu/CTRL_count/outs/ctrl_atac_peaks.bed"
PATH_BARCODE_METRICS = "~/hpchome/projects/2023_10X_multiomics_E15.5_murine_skin/raw_data/dnacore454.healthcare.uiowa.edu/CTRL_count/outs/per_barcode_metrics.csv"
PATH_ATAC_FRAGMENTS = "~/hpchome/projects/2023_10X_multiomics_E15.5_murine_skin/raw_data/dnacore454.healthcare.uiowa.edu/CTRL_count/outs/atac_fragments.tsv.gz"
PATH_RNA = "~/hpchome/projects/2023_10X_multiomics_E15.5_murine_skin/raw_data/dnacore454.healthcare.uiowa.edu/CTRL_count/outs/filtered_feature_bc_matrix.h5"
OUTPUT_RDS = "~/hpchome/projects/CTRL_only/CTRL_only/Loadable Data/ctrl.rds"

#Read in peak sets
ctrl.peaks <- read.table(file = PATH_ATAC_PEAKS, col.names = c("chr", "start", "end"))

#Convert to genomic ranges
ctrl.gr <- makeGRangesFromDataFrame(ctrl.peaks)

#Filter out bad peaks based on length
peakwidths <- width(ctrl.gr)
ctrl.gr <- ctrl.gr[peakwidths < 10000 & peakwidths > 20]

#Create fragment objects
ctrl.md <- read.table(file = PATH_BARCODE_METRICS, stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)

#Add fragpath
ctrl.fragpath <- PATH_ATAC_FRAGMENTS

#Create fragment object
ctrl.fragobj <- CreateFragmentObject(path = ctrl.fragpath)

#Add ctrl.gr 
ctrl.counts <- FeatureMatrix(fragments = ctrl.fragobj, features = ctrl.gr, cells = rownames(ctrl.md))

#Get gene annotations for genome (in this case mm10), set to correct annotation style
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

#Change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

#Create chromatin assay
ctrl.assay <- CreateChromatinAssay(ctrl.counts, fragments = ctrl.fragobj, annotation = annotations, genome = "mm10")

#Create seurat object
ctrl <- CreateSeuratObject(ctrl.assay, assay = "ATAC", meta.data = ctrl.md, project = "CTRL")

#Add the gene information to the object
Annotation(ctrl) <- annotations

#Add RNA data
ctrl.rna <- Read10X_h5(PATH_RNA)

#Extra step from seurat issues 7631
ctrl_rna_names <- colnames(ctrl.rna$`Gene Expression`)
ctrl_RNA_names <- colnames(ctrl.assay)
ctrl_intersect <- intersect(ctrl_RNA_names, ctrl_rna_names)
ctrl <- subset(ctrl, cells = ctrl_intersect)

#Add RNA data
ctrl[['RNA']] <- CreateAssayObject(counts = ctrl.rna$`Gene Expression`[ , ctrl_intersect])

#CTRL ATAC
DefaultAssay(ctrl) <- "ATAC"
ctrl <- NucleosomeSignal(ctrl)
ctrl <- TSSEnrichment(ctrl)

#Add blacklist ratio and fraction of reads in peaks
ctrl$pct_reads_in_peaks <- ctrl$atac_peak_region_fragments / ctrl$atac_fragments * 100

ctrl$blacklist_fraction <- FractionCountsInRegion(
  object = ctrl, 
  assay = 'ATAC',
  regions = blacklist_mm10
)

#Plot variables 
VlnPlot(object = ctrl, features = c("nCount_ATAC", "nCount_RNA", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks", "blacklist_fraction", "atac_peak_region_fragments"), ncol = 4, pt.size = 0.5)
ggsave("ctrl_QC_VlnPlot.tiff", width = 12, height = 8, units = c("in"), dpi = 300)

#Filter out low quality cells (change these numbers to suit your data visualized in the above VlnPlot)
ctrl <- subset(x = ctrl, subset = nCount_ATAC < 40000 & nCount_RNA < 30000 & nCount_ATAC > 100 & nCount_RNA > 100 & nucleosome_signal < 4 & TSS.enrichment > 1 & TSS.enrichment < 30 & pct_reads_in_peaks > 15 & blacklist_fraction < 0.1 & atac_peak_region_fragments <40000)

#Check variables after 
VlnPlot(object = ctrl, features = c("nCount_ATAC", "nCount_RNA", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks", "blacklist_fraction", "atac_peak_region_fragments"), ncol = 4, pt.size = 0.5)
ggsave("ctrl_QC_after_VlnPlot.tiff", width = 12, height = 8, units = c("in"), dpi = 300)

#CTRL RNA
DefaultAssay(ctrl) <- "RNA"
ctrl <- SCTransform(ctrl)
ctrl <- RunPCA(ctrl)
# You can change the number of neighbors using n.neighbors = X in the line below (usually 5-50, default is 30)
ctrl <- RunUMAP(ctrl, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

#More CTRL ATAC
DefaultAssay(ctrl) <- "ATAC"
ctrl <- FindTopFeatures(ctrl, min.cutoff = 5)
ctrl <- RunTFIDF(ctrl)
ctrl <- RunSVD(ctrl)
# You can change the number of neighbors using n.neighbors = X in the line below (usually 5-50, default is 30)
ctrl <- RunUMAP(ctrl, reduction = 'lsi', dims = 2:50, reduction.name = "umap.ATAC", reduction.key = "atacUMAP_")

#Weighted nearest neighbors
ctrl <- FindMultiModalNeighbors(ctrl, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), k.nn = 30)
ctrl <- RunUMAP(ctrl, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ctrl <- FindClusters(ctrl, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.9)

#Visualization
ctrl.dimplot <- DimPlot(object = ctrl, reduction = 'wnn.umap', label = T, label.size = 4.5, pt.size = 0.35)
ctrl.dimplot
ggsave("ctrl.tiff", ctrl.dimplot, height=7, width=10)

#Save files
saveRDS(ctrl, file = OUTPUT_RDS)