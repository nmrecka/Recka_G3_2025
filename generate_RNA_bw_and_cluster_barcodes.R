library(Signac)
library(Seurat)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)
library(purrr)
library(future)
library(future.apply)
library(Matrix)
library(pbapply)

# Define the function to export bigwig files 
ExportGroupBW  <- function(
    object,
    assay = NULL,
    group.by = NULL,
    idents = NULL,
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    cutoff = NULL,
    chromosome = NULL,
    outdir = NULL,
    verbose=TRUE
) {
  # Check if temporary directory exist
  if (!dir.exists(outdir)){
    dir.create(outdir)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) { 
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/") 
    return(NULL) 
  }
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  DefaultAssay(object = object) <- assay
  group.by <- SetIfNull(x = group.by, y = 'ident')
  Idents(object = object) <- group.by
  idents <- SetIfNull(x = idents, y = levels(x = object))
  GroupsNames <- names(x = table(object[[group.by]])[table(object[[group.by]]) > minCells])
  GroupsNames <- GroupsNames[GroupsNames %in% idents]
  # Check if output files already exist
  lapply(X = GroupsNames, FUN = function(x) {
    fn <- paste0(outdir, .Platform$file.sep, x, ".bed")
    if (file.exists(fn)) {
      message(sprintf("The group \"%s\" is already present in the destination folder and will be overwritten !",x))
      file.remove(fn)
    }
  })      
  # Splitting fragments file for each idents in group.by
  SplitFragments(
    object = object,
    assay = assay,
    group.by = group.by,
    idents = idents,
    outdir = outdir,
    file.suffix = "",
    append = TRUE,
    buffer_length = 256L,
    verbose = verbose
  )
  # Column to normalized by
  if(!is.null(x = normMethod)) {
    if (tolower(x = normMethod) %in% c('rc', 'ncells', 'none')){
      normBy <- normMethod
    } else{
      normBy <- object[[normMethod, drop = FALSE]]
    }
  }
  # Get chromosome information
  if(!is.null(x = chromosome)){
    seqlevels(object) <- chromosome
  }
  availableChr <- names(x = seqlengths(object))
  chromLengths <- seqlengths(object)
  chromSizes <- GRanges(
    seqnames = availableChr,
    ranges = IRanges(
      start = rep(1, length(x = availableChr)),
      end = as.numeric(x = chromLengths)
    )
  )
  if (verbose) {
    message("Creating tiles")
  }
  # Create tiles for each chromosome, from GenomicRanges
  tiles <- unlist(
    x = slidingWindows(x = chromSizes, width = tileSize, step = tileSize)
  )
  if (verbose) {
    message("Creating bigwig files at ", outdir)
  }
  # Run the creation of bigwig for each cellgroups
  if (nbrOfWorkers() > 1) { 
    mylapply <- future_lapply 
  } else { 
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply) 
  }
  
  covFiles <- mylapply(
    GroupsNames,
    FUN = CreateBWGroup,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    tileSize,
    normMethod,
    cutoff,
    outdir
  )
  return(covFiles)
}

CreateBWGroup <- function(
    groupNamei,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    tileSize,
    normMethod,
    cutoff,
    outdir
) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) { 
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/") 
    return(NULL) 
  }
  normMethod <- tolower(x = normMethod)
  # Read the fragments file associated to the group
  fragi <- rtracklayer::import(
    paste0(outdir, .Platform$file.sep, groupNamei, ".bed"), format = "bed"
  )
  cellGroupi <- unique(x = fragi$name)
  # Open the writing bigwig file
  covFile <- file.path(
    outdir,
    paste0(groupNamei, "-TileSize-",tileSize,"-normMethod-",normMethod,".bw")
  )
  
  covList <- lapply(X = seq_along(availableChr), FUN = function(k) {
    fragik <- fragi[seqnames(fragi) == availableChr[k],]
    tilesk <- tiles[BiocGenerics::which(S4Vectors::match(seqnames(tiles), availableChr[k], nomatch = 0) > 0)]
    if (length(x = fragik) == 0) {
      tilesk$reads <- 0
      # If fragments
    } else {
      # N Tiles
      nTiles <- chromLengths[availableChr[k]] / tileSize
      # Add one tile if there is extra bases
      if (nTiles%%1 != 0) {
        nTiles <- trunc(x = nTiles) + 1
      }
      # Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(fragik)$name, cellGroupi)
      
      # For each tiles of this chromosome, create start tile and end tile row,
      # set the associated counts matching with the fragments
      mat <- sparseMatrix(
        i = c(trunc(x = start(x = fragik) / tileSize),
              trunc(x = end(x = fragik) / tileSize)) + 1,
        j = as.vector(x = c(matchID, matchID)),
        x = rep(1, 2*length(x = fragik)),
        dims = c(nTiles, length(x = cellGroupi))
      )
      
      # Max count for a cells in a tile is set to cutoff
      if (!is.null(x = cutoff)){
        mat@x[mat@x > cutoff] <- cutoff
      }
      # Sums the cells
      mat <- rowSums(x = mat)
      tilesk$reads <- mat
      # Normalization
      if (!is.null(x = normMethod)) {
        if (normMethod == "rc") {
          tilesk$reads <- tilesk$reads * 10^4 / length(fragi$name)
        } else if (normMethod == "ncells") {
          tilesk$reads <- tilesk$reads / length(cellGroupi)
        } else if (normMethod == "none") {
        } else {
          if (!is.null(x = normBy)){
            tilesk$reads <- tilesk$reads * 10^4 / sum(normBy[cellGroupi, 1])
          }
        }
      }
    }
    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]
    tilesk
  })
  
  names(covList) <- availableChr
  covList <- as(object = covList, Class = "RleList")
  rtracklayer::export.bw(object = covList, con = covFile)
  return(covFile)
}

SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

# Following is for exporting ATAC data bigwig files for individual clusters and combined clusters

TILE_SIZE <- 25
MIN_CELLS <- 5

# load the CTRL object
ctrl <- readRDS("~/hpchome/projects/CTRL_only/CTRL_only/Loadable Data/ctrl.rds")
output_dir <- "~/hpchome/projects/CTRL_only"

# set the default assay ATAC
DefaultAssay(ctrl) <- "ATAC"

# export ATAC pseudo-bulk bigwig files for individual clusters
SplitFragments(
  object = ctrl,
  assay = "ATAC",
  group.by = "seurat_clusters",
  idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),
  outdir = output_dir,
  file.suffix = "",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE
)

ExportGroupBW(object = ctrl,
              assay = 'ATAC',
              group.by = ',seurat_clusters', # or 'combined_cluster' for combined clusters,
              idents = NULL, #NULL [for all clusters] or list(seurat_clusters=c(21)) [for selected clusters],
              normMethod = "RC", #"TSS.enrichment", #'ncells', 'none' or any quantitative values from @meta.data
              tileSize = TILE_SIZE,
              minCells = MIN_CELLS,
              cutoff = NULL,
              chromosome = standardChromosomes(ctrl),
              outdir = output_dir,
              verbose=TRUE
)


# save the barcodes of each cluster
barcode_list <- list()
individual_clusters <- unique(ctrl$seurat_clusters)
for (cluster in individual_clusters) {
  barcode_list[[as.character(cluster)]] <- Cells(ctrl)[ctrl$seurat_clusters == cluster]
}
for (group in names(barcode_list)) {
  write.table(
    barcode_list[[group]],
    file = paste0(output_dir, group, "_barcodes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}


# export ATAC pseudo-bulk bigwig files for combined clusters

ctrl$combined_cluster <- as.character(ctrl$seurat_clusters)

# define the combined clusters
ctrl$combined_cluster[ctrl$seurat_clusters %in% c(2,5)] <- "basal"
ctrl$combined_cluster[ctrl$seurat_clusters %in% c(18)] <- "dermal_condensate"
ctrl$combined_cluster[ctrl$seurat_clusters %in% c(3,10)] <- "spinous"
ctrl$combined_cluster[ctrl$seurat_clusters %in% c(7)] <- "differentiated"
ctrl$combined_cluster[ctrl$seurat_clusters %in% c(8)] <- "hair_follicle"
ctrl$combined_cluster[ctrl$seurat_clusters %in% c(2,3,5,7,8,10)] <- "epidermis"
ctrl$combined_cluster[ctrl$seurat_clusters %in% c(0,1,4,6,9)] <- "dermis"
ctrl$combined_cluster[ctrl$seurat_clusters %in% c(11,12,13,14,15,16,17,19,20,21,22)] <- "other"

ExportGroupBW(object = ctrl,
              assay = 'ATAC',
              group.by = 'combined_cluster', # or 'seurat_clusters' for individual clusters,
              idents = NULL, #NULL [for all clusters] or list(seurat_clusters=c(21)) [for selected clusters],
              normMethod = "RC", #"TSS.enrichment", #'ncells', 'none' or any quantitative values from @meta.data
              tileSize = TILE_SIZE,
              minCells = MIN_CELLS,
              cutoff = NULL,
              chromosome = standardChromosomes(ctrl),
              outdir = output_dir,
              verbose=TRUE
)

# save the barcodes of each combined cluster
barcode_list <- list()
combined_clusters <- unique(ctrl$combined_cluster)
for (cluster in combined_clusters) {
  barcode_list[[cluster]] <- Cells(ctrl)[ctrl$combined_cluster == cluster]
}
for (group in names(barcode_list)) {
  write.table(
    barcode_list[[group]],
    file = paste0(output_dir, group, "_barcodes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}
