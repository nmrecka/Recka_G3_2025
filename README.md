# Recka_G3_2025
Code for Recka et al. 2025 published in G3

## Requirements
To run the code in this repository, you will need the following:

### R Libraries
- Signac
- Seurat
- GenomicRanges
- rtracklayer
- dplyr
- data.table
- purrr
- future
- future.apply
- Matrix
- pbapply
- ggplot2
- cowplot
- biovizBase
- hdf5r
- EnhancedVolcano
- EnsDb.Mmusculus.v79
- JASPAR2020
- TFBSTools
- BSgenome.Mmusculus.UCSC.mm10
- patchwork

### Command-line Tools
- samtools
- bamCoverage (from deepTools)

## Instructions

### Preprocessing
1. **Preprocessing**:
   - Use the `preprocessing.R` script to prepare the data. This script processes ATAC and RNA data, filters low-quality cells, and creates a Seurat object.

### Generating Plots
1. **Gene Expression Visualization**:
   - Use the `figure_generation.R` script to generate various plots to observe gene expression and chromatin accessibility data (feature plot, dot plot, volcano plot, etc.). 
   - Example: Subset clusters and define markers for visualization. Use `DotPlot` for marker visualization and save the plots.
   - Example: Perform differential analysis using the `FindMarkers` function. Visualize results with volcano plots using the `EnhancedVolcano` package.

### Generating BigWig Files for RNA and ATAC
1. **Exporting BigWig Files**:
   - Use the `generate_RNA_bw_and_cluster_barcodes.R` script to export BigWig files for RNA and ATAC data. The `ExportGroupBW` function is used for this purpose.
   - For ATAC data, ensure the default assay is set to "ATAC" and use the `SplitFragments` function to split fragments by clusters before exporting.

2. **Normalizing BAM to BigWig**:
   - Use the `normalize_bam_to_bw.sh` script to normalize BAM files and generate BigWig files. Example usage:
     ```bash
     ./normalize_bam_to_bw.sh --input-bam sample.bam --output-bw sample.bw --scale-factor 1 --bin-size 25
     ```

### Additional Notes
- Ensure all required libraries and tools are installed before running the scripts.
- Modify file paths in the scripts to match your local setup.
