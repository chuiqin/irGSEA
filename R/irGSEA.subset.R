#' Subset the enrichment score matrix
#'
#' Subset the enrichment score matrix
#'
#' @param object Seurat object (V4 or V5).
#' @param features Vector. The geneset want to subset.
#' @param invert If filp the geneset want to subset.
#' @param method A vector. Default c("AUCell", "UCell", "singscore", "ssgsea",
#' "JASMINE", "viper").
#'
#'
#' @return Seurat object including score matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' # load PBMC dataset by R package SeuratData
#  # devtools::install_github('satijalab/seurat-data')
#' library(Seurat)
#' library(SeuratData)
#' # download 3k PBMCs from 10X Genomics
#' InstallData("pbmc3k")
#' library(Seurat)
#'
#' library(RcppML)
#' library(irGSEA)
#' library(tidyverse)
#' library(clusterProfiler)
#' data("pbmc3k.final")
#' pbmc3k.final <- SeuratObject::UpdateSeuratObject(pbmc3k.final)
#'
#' # download gmt file
#' gmt_url1 <- "https://data.broadinstitute.org/"
#' gmt_url2 <- "gsea-msigdb/msigdb/release/2023.2.Hs/",
#' gmt_url3 <- "h.all.v2023.2.Hs.symbols.gmt"
#' gmt_url <- paste0(gmt_url1, gmt_url2, gmt_url3)
#' local_gmt <- "./h.all.v2023.2.Hs.symbols.gmt"
#' download.file(gmt_url , local_gmt)
#' msigdb <- clusterProfiler::read.gmt("./h.all.v2023.2.Hs.symbols.gmt")
#'
#' # convert to list[hallmarker] required by irGSEA package
#' msigdb$term <- factor(msigdb$term)
#' msigdb <- msigdb %>%
#'   dplyr::group_split(term, .keep = F) %>%
#'   purrr::map( ~.x %>% dplyr::pull(gene) %>% unique(.)) %>%
#'   purrr::set_names(levels(msigdb$term))
#'
#' pbmc3k.final1 <- irGSEA.score(object = pbmc3k.final, assay = "RNA", slot = "data",
#'                               custom = T, geneset = msigdb[1:25],
#'                               method = c("AUCell", "UCell", "singscore", "ssgsea"),
#'                               kcdf = 'Gaussian')
#' pbmc3k.final2 <- irGSEA.score(object = pbmc3k.final, assay = "RNA", slot = "data",
#'                               custom = T, geneset = msigdb[26:50],
#'                               method = c("AUCell", "UCell", "singscore", "ssgsea"),
#'                               kcdf = 'Gaussian')
#'
#' pbmc3k.final3 <- irGSEA.merge(object.x = pbmc3k.final1,
#'                               object.y = pbmc3k.final2,
#'                               method = c("AUCell", "UCell", "singscore", "ssgsea"),
#'                               overwrite = T)
#' pbmc3k.final4 <- irGSEA.subset(object = pbmc3k.final3,
#'                                features = rownames(pbmc3k.final3[["AUCell"]])[1:25],
#'                                method = c("AUCell", "UCell", "singscore", "ssgsea"),
#'                                invert = F)
#' }
#'
irGSEA.subset <- function(object = NULL,
                          features = NULL,
                          invert = F,
                          method = c("AUCell", "UCell", "singscore",
                                    "ssgsea","JASMINE", "viper")
                         ){

  # work
  for (k in method) {
    # if v5 or v3
    if (class(object)[1] == "Assay5") {
      options(Seurat.object.assay.version = "v5")
      version <- "v5"
    }else{
      options(Seurat.object.assay.version = "v3")
      version <- "v3"
    }
    if (invert) {
      features <- rownames(object[[k]])[!rownames(object[[k]]) %in% features]
    }

    object[[k]] <- subset(object[[k]], features = features)

  }


  return(object)

}


