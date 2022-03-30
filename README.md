# SCHNAPPsContributions #


some contributions to the SCHNAPPs project (https://c3bi-pasteur-fr.github.io/UTechSCB-SCHNAPPs/)

## Trajectory related contributions

### Tempora: cell trajectory inference using time-series single-cell RNA sequencing data

https://github.com/BaderLab/Tempora
https://www.baderlab.org/Software/Tempora

Tempora is a novel cell trajectory inference method that orders cells using time information from time-series scRNAseq data. Tempora uses biological pathway information to help identify cell type relationships and can identify important time-dependent pathways to help interpret the inferred trajectory.

### scorpius

https://github.com/rcannood/SCORPIUS
https://www.biorxiv.org/content/10.1101/079509v2

SCORPIUS an unsupervised approach for inferring linear developmental chronologies from single-cell RNA sequencing data. 

### ELPIGraph Elastic principal graphs

https://github.com/Albluca/ElPiGraph.R

The mapping of a principal graph into multidimensional data space is regularized by minimizing the stretching of the graph edges and the deviation from harmonicity for the graph stars.


## CelliD

https://github.com/RausellLab/CelliD

[Gene signature extraction and cell identity recognition at the single-cell level with Cell-ID, Nature Biotechnology 2021](https://rdcu.be/cjFWE)

CelliD is a robust statistical method that performs gene signature extraction and functional annotation for each individual cell in a single-cell RNA-seq dataset.

## Self-organizing Maps (SOM) based on the Rsomoclu, and kohonen packages.

Unsupervised clustering of genes.


## scDEA 

https://github.com/Zhangxf-ccnu/scDEA

https://pubmed.ncbi.nlm.nih.gov/34571530/

applies ensemble learning for single-cell differential expression analysis on single-cell RNA-seq dataset.

This is actually implemented in the main up, but only available if scDEA is installed. Since it is very heavy computationally it is documented only here.

# Prepare data for SCHNAPPs

## read data from 10X experiments

load some utility functions. 

```{R}
source("functions.R")


seuratList = list() # list of count objects coming from the modified Seurat H5 reader

sampleNames <- c("Buffercells", "GFPcells")

# PROJECTFOLDER = 10X project folder that holds the 'counts' directory
# we are reading filtered H5 files.

for (sIdx in 1:2){
  h5file <- paste0("PROJECTFOLDER/counts/",sampleNames[sIdx],"/outs/filtered_feature_bc_matrix.h5")
  seuratList[[sIdx]] <- myRead10X_h5(h5file, use.names = TRUE, unique.features = TRUE, sampleName = sampleNames[sIdx])
}
scexSeurat = seuratList[[1]]
for (sIdx in 2){
  scexSeurat <- cbind(scexSeurat, seuratList[[sIdx]])
}
# print the result
scexSeurat

# filter genes with no expression
scexSeurat = scexSeurat[rowSums(as.matrix(scexSeurat))>0,]
# make sure that no genes exist with the same overall expression. This would cause problems for some calculations.
scexSeurat = unique(as.matrix(scexSeurat))



scEx <- SingleCellExperiment(
  assay = list(counts = scexSeurat),
  colData = pd,
  rowData = featuredata
)





featureData_summary <- data.frame(
  "Description" = NA,
  "gene_id" = my_gene$gene_id,
  "Chromosome.Name" = my_gene$seqid,
  "Associated.Gene.Name" = my_gene$gene_name,
  stringsAsFactors = F
)
# featureData_summary <- featureData_summary[, -2]
featuredata <- data.frame(featureData_summary, stringsAsFactors = F)
featuredata$id <- featuredata$Associated.Gene.Name
featuredata$symbol <- make.unique(featuredata$Associated.Gene.Name)
featuredata <- featuredata[which(featuredata$symbol %in% rownames(scexSeurat)), ]
featuredata <- unique(featuredata)
featuredata[which(duplicated(featuredata$symbol)), ]

# rownames(featuredata) = 1:nrow(featuredata)
# featureDatatmp <- featuredata[, -2]
# 
# featureDatatmp <- unique(featureDatatmp)

#same symbol, different ensg numbers:
# featuredata = featuredata[rownames(featureDatatmp), ]
rownames(featuredata) <- featuredata$symbol
nrow(featuredata)
nrow(scexSeurat)
nrFD = nrow(featuredata)
newRowNames = rownames(scexSeurat)[which(!rownames(scexSeurat) %in% featuredata$Associated.Gene.Name)]

if (length(newRowNames) > 0) {
  featuredata[(nrFD + 1):(nrFD + length(newRowNames)),] <- NA
  rownames(featuredata[(nrFD + 1):(nrFD + length(newRowNames)),]) <- newRowNames
  featuredata[(nrFD + 1):(nrFD + length(newRowNames)),"symbol"] <- newRowNames
  featuredata[(nrFD + 1):(nrFD + length(newRowNames)),"Description"] <- newRowNames
  featuredata[(nrFD + 1):(nrFD + length(newRowNames)),"Associated.Gene.Name"] <- newRowNames
  featuredata[(nrFD + 1):(nrFD + length(newRowNames)),"id"] <- newRowNames
}

rownames(scexSeurat)[which(!rownames(scexSeurat) %in% featuredata$Associated.Gene.Name)]
rownames(scexSeurat)[which(!featuredata$Associated.Gene.Name %in% rownames(scexSeurat)) ]
# scexSeurat <- scexSeurat[featuredata$Associated.Gene.Name, ]
featuredata <- featuredata[rownames(scexSeurat), ]
# featuredata$Description = featureData_summary93[featuredata$gene_id,"Description"]
nrow(featuredata)
nrow(scexSeurat)


pd <- data.frame(
  barcode = sub("(.*)-(.*)", "\\1", colnames(scexSeurat)),
  sampleNames = sub(".*-(.*)", "\\1", colnames(scexSeurat))
)
pd$barcode <- as.character(pd$barcode)
rownames(pd) <- colnames(scexSeurat)

hs.pairs <- readRDS("mouse_cycle_markers.rds")
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(featuredata), keytype="SYMBOL", column="ENSEMBL")
assignments <- cyclone(scexSeurat, hs.pairs, gene.names=ensembl, BPPARAM = MulticoreParam())
pd$phases = assignments$phases
pd$G1score = assignments$normalized.scores$G1
pd$Sscore = assignments$normalized.scores$S
pd$G2Mscore = assignments$normalized.scores$G2M

pd$phases[is.na(pd$phases)] = "NA"
pd$G1score[is.na(pd$G1score)] = "NA"
pd$Sscore[is.na(pd$Sscore)] = "NA"
pd$G2Mscore[is.na(pd$G2Mscore)] = "NA"


# scExTemp <- applySingleR(scEx, DatabaseImmuneCellExpressionData(), "cellTypes")

# colData(scEx)$cellTypes = as.factor(colData(scExTemp)$cellTypes)

fname = "Ferdinand-scRNAseq"
outfile <- paste0(fname, ".RData")
save(file = outfile, list = c("scEx"))

# set.seed(1)
# colIdx <- sample(1:ncol(scexSeurat), 2000, replace = FALSE)
# scEx <- SingleCellExperiment(
#   assay = list(counts = scexSeurat[, colIdx]),
#   colData = pd[colIdx, ],
#   rowData = featuredata
# )
# outfile <- paste0(fname, ".sml.RData")
# save(file = outfile, list = c("scEx"))

require(rmarkdown)
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  echo = FALSE,
  include = TRUE
)
getwd()
rm("params")
render("gbmReport.Rmd",
       output_file = paste0(fname, ".report.html"),
       output_format = "html_document",
       params = list(
         fileN = paste0(fname, ".RData"),
         min.genes = 2,
         min.cells = 3,
         low.thres1 = 2,
         low.thres2 = -Inf,
         high.thres1 = 2500,
         high.thres2 = 2000
       )
)

```
