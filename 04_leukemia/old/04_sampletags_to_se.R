#!/usr/bin/env R
##
## Parses the CB's leukemia sampletags into a summarizedexp R object
##
## 23rd Nov 2023
## Izaskun Mallona

library(SingleCellExperiment)
library(Matrix)

read_matrix <- function(mtx, cells, features, cell.column = 1, feature.column = 1, modality, wta_whitelist) {
  cell.barcodes <- read.table(
    file = cells,
    header = FALSE,
    row.names = cell.column)


  feature.names <- read.table(
    file = features,
    header = FALSE,
    row.names = feature.column)

  d <- as(readMM(mtx), 'CsparseMatrix')

  ## if (modality == 'wta') {
  ## colnames(d) <- gsub('_', '', rownames(cell.barcodes))
  ## } else if (modality == 'tso') {
  ##     ## remove the fixed parts of the TSO CBs
  ##     colnames(d) <- paste0(
  ##         substr(rownames(cell.barcodes), 1, 9),
  ##         substr(rownames(cell.barcodes), 9+4+1, 9+4+9),
  ##         substr(rownames(cell.barcodes), 9+4+9+4+1, 9+4+9+4+9))
  ## }

  colnames(d) <- rownames(cell.barcodes)
  rownames(d) <- rownames(feature.names)
  
  if (modality == 'tso') {
      d <- d[,wta_whitelist]
  }

  return(d)
}

## sampletag stuff
wd <- '/home/imallona/ebrunner_spectral/leukemia_sampletags'
raw <- file.path(wd, 'patients_wta_vs_sampletags', 'Solo.out', 'Gene', 'raw')

se <- SummarizedExperiment(read_matrix(file.path(raw, 'matrix.mtx'),
                                       file.path(raw, 'barcodes.tsv'),
                                       file.path(raw, 'features.tsv'),
                                       modality = 'wta'))

saveRDS(se, file.path(wd, 'sampletags_summarizedexperiment.rds'))
