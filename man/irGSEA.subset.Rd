% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irGSEA.subset.R
\name{irGSEA.subset}
\alias{irGSEA.subset}
\title{Subset the enrichment score matrix}
\usage{
irGSEA.subset(
  object = NULL,
  features = NULL,
  invert = F,
  method = c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper")
)
}
\arguments{
\item{object}{Seurat object (V4 or V5).}

\item{features}{Vector. The geneset want to subset.}

\item{invert}{If filp the geneset want to subset.}

\item{method}{A vector. Default c("AUCell", "UCell", "singscore", "ssgsea",
"JASMINE", "viper").}
}
\value{
Seurat object including score matrix.
}
\description{
Subset the enrichment score matrix
}
\examples{
\dontrun{
# load PBMC dataset by R package SeuratData
library(Seurat)
library(SeuratData)
# download 3k PBMCs from 10X Genomics
InstallData("pbmc3k")
library(Seurat)

library(RcppML)
library(irGSEA)
library(tidyverse)
library(clusterProfiler)
data("pbmc3k.final")
pbmc3k.final <- SeuratObject::UpdateSeuratObject(pbmc3k.final)

# download gmt file
gmt_url1 <- "https://data.broadinstitute.org/"
gmt_url2 <- "gsea-msigdb/msigdb/release/2023.2.Hs/",
gmt_url3 <- "h.all.v2023.2.Hs.symbols.gmt"
gmt_url <- paste0(gmt_url1, gmt_url2, gmt_url3)
local_gmt <- "./h.all.v2023.2.Hs.symbols.gmt"
download.file(gmt_url , local_gmt)
msigdb <- clusterProfiler::read.gmt("./h.all.v2023.2.Hs.symbols.gmt")

# convert to list[hallmarker] required by irGSEA package
msigdb$term <- factor(msigdb$term)
msigdb <- msigdb \%>\%
  dplyr::group_split(term, .keep = F) \%>\%
  purrr::map( ~.x \%>\% dplyr::pull(gene) \%>\% unique(.)) \%>\%
  purrr::set_names(levels(msigdb$term))

pbmc3k.final1 <- irGSEA.score(object = pbmc3k.final, assay = "RNA", slot = "data",
                              custom = T, geneset = msigdb[1:25],
                              method = c("AUCell", "UCell", "singscore", "ssgsea"),
                              kcdf = 'Gaussian')
pbmc3k.final2 <- irGSEA.score(object = pbmc3k.final, assay = "RNA", slot = "data",
                              custom = T, geneset = msigdb[26:50],
                              method = c("AUCell", "UCell", "singscore", "ssgsea"),
                              kcdf = 'Gaussian')

pbmc3k.final3 <- irGSEA.merge(object.x = pbmc3k.final1,
                              object.y = pbmc3k.final2,
                              method = c("AUCell", "UCell", "singscore", "ssgsea"),
                              overwrite = T)
pbmc3k.final4 <- irGSEA.subset(object = pbmc3k.final3,
                               features = rownames(pbmc3k.final3[["AUCell"]])[1:25],
                               method = c("AUCell", "UCell", "singscore", "ssgsea"),
                               invert = F)
}

}
