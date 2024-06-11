#' Merge the enrichment score assay
#'
#' Merge the enrichment score assay among various Seurat objects
#'
#' @param object.x Seurat object (V4 or V5).
#' @param object.y The list includes various Seurat objects. The verison of
#' Seurat object should be the same as the version of object.x.
#' @param method A vector. Default c("AUCell", "UCell", "singscore", "ssgsea",
#' "JASMINE", "viper").
#' @param overwrite Default True. The same geneset name exists in two gene scoring
#' matrices, the newly added geneset will overwrite the previous geneset if
#' the overwrite is true. The newly added geneset will be forcibly renamed
#' as "geneset's name + serial number" if the overwrite is false.
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
#'
#' }
#'
irGSEA.merge <- function(object.x = NULL, object.y = NULL,
                         method = c("AUCell", "UCell", "singscore",
                                    "ssgsea","JASMINE", "viper"),
                         overwrite = T
                         ){
  # create a list
  if (length(object.y) == 1) {
    object.y <- list(object.y)
  }

  # work
  for (k in method) {

    for (i in seq_along(object.y)) {

      # acts.x
      acts.x <- tryCatch({
        acts.x <- SeuratObject::GetAssayData(object.x, assay = k, slot = "scale.data")
      }, error = function(e) {
        acts.x <- NULL
      })

      # acts.y
      acts.y <- tryCatch({
        acts.y <- SeuratObject::GetAssayData(object.y[[i]], assay = k, slot = "scale.data")
      }, error = function(e) {
        acts.y <- NULL
      })

      # next
      if (is.null(acts.x) & is.null(acts.y)) {
        next
      }



      # merge, if overwrite, or renames the same geneset
      if (overwrite) {
        index.intersect <- !rownames(acts.x) %in% rownames(acts.y)
        acts <- rbind(acts.x[index.intersect,], acts.y)


      }else{
        index.intersect <- match(rownames(acts.x), rownames(acts.y))
        if (any(!is.na(index.intersect))) {
          index.intersect2 <- index.intersect[!is.na(index.intersect)]
          rownames(acts.y)[index.intersect2] <- paste0(rownames(acts.y)[index.intersect2], "-1")
          acts <- rbind(acts.x, acts.y)
        }else{
          acts <- rbind(acts.x, acts.y)
        }


      }

      if (is.null(acts)) { next }


      # if version
      version <- tryCatch({
        # if v5 or v3
        if (class(object.x[[k]])[1] == "Assay5") {
          options(Seurat.object.assay.version = "v5")
          version <- "v5"
        }else{
          options(Seurat.object.assay.version = "v3")
          version <- "v3"
        }
      }, error = function(e) {
        version <- tryCatch({
          # if v5 or v3
          if (class(object.y[[i]][[k]])[1] == "Assay5") {
            options(Seurat.object.assay.version = "v5")
            version <- "v5"
          }else{
            options(Seurat.object.assay.version = "v3")
            version <- "v3"
          }
        }, error = function(e) {
          options(Seurat.object.assay.version = "v3")
          version <- "v3"
        })
      })

      # if verison of Assay same
      tryCatch({
        if (class(object.x[[k]])[1] != class(object.y[[i]][[k]])[1]) {
          message(paste0("The Seurat/Assay versions of the object.x and object.y: ",
                         i,
                         " are inconsistent. \n",
                         "We convert object.y: ",
                         i,
                         " from ",
                         class(object.y[[i]][[k]])[1],
                         " to ",
                         class(object.x[[k]])[1],
                         " ."))

          # convert
          if (class(object.y[[i]][[k]])[1] == "Assay5") {
            object.y[[i]][[k]] <- methods::as(object = object.y[[i]][[k]], Class = "Assay")
          }else{
            object.y[[i]][[k]] <- methods::as(object = object.y[[i]][[k]], Class = "Assay5")
          }
        }
      }, error = function(e) {
        # print("")
      })




      # meta.features
      # object.x

      if (is.null(acts.x)) {
        object.x.meta.features <- NULL
      }else{
        if (purrr::is_empty(object.x[[k]]@meta.features)) {
          if (version == "v3") {
            object.x.meta.features <- data.frame(geneset = rownames(object.x[[k]]),
                                                 target.gene = "")
          }else{
            object.x.meta.features <- data.frame(geneset = rownames(object.x[[k]]),
                                                 target.gene = "") %>%
              tibble::column_to_rownames(var = "geneset")
          }

        }else{
          if (version == "v3") {
            object.x.meta.features <- object.x[[k]]@meta.features %>%
              tibble::rownames_to_column(var = "geneset")
          }else{
            object.x.meta.features <- object.x[[k]]@meta.features
          }
        }
        if (overwrite) {
          index.intersect <- !rownames(acts.x) %in% rownames(acts.y)
          object.x.meta.features <- object.x.meta.features[index.intersect, ]

        }
      }



      # object.y

      if (is.null(acts.y)) {
        object.y.meta.features <- NULL
      }else{
        if (purrr::is_empty(object.y[[i]][[k]]@meta.features)) {

          if (version == "v3") {
            object.y.meta.features <- data.frame(geneset = rownames(object.y[[i]][[k]]),
                                                 target.gene = "")
          }else{
            object.y.meta.features <- data.frame(geneset = rownames(object.y[[i]][[k]]),
                                                 target.gene = "") %>%
              tibble::column_to_rownames(var = "geneset")
          }
          object.y.meta.features$geneset <- rownames(acts.y)

        }else{
          if (version == "v3") {
            object.y.meta.features <- object.y[[i]][[k]]@meta.features %>%
              tibble::rownames_to_column(var = "geneset")
          }else{
            object.y.meta.features <- object.y[[i]][[k]]@meta.features
          }
          object.y.meta.features$geneset <- rownames(acts.y)
        }

      }




      if (is.null(object.x.meta.features) & is.null(object.y.meta.features) ) {
        next
      }



      # add matrix meta.features

      object.x[[k]] <- SeuratObject::CreateAssayObject(counts = acts)
      object.x <- SeuratObject::SetAssayData(object.x, slot = "scale.data",
                                             new.data = acts,
                                             assay = k)
      if (utils::packageVersion("Seurat") >= "5.0.0") {
        object.x[[k]]$counts <- NULL}



      if (version == "v3") {
        object.x[[k]]@meta.features <- as.data.frame(rbind(object.x.meta.features,
                                                           object.y.meta.features)) %>%
          tibble::column_to_rownames(var = "geneset")

      }else{
        object.x[[k]]@meta.features <- rbind(object.x.meta.features,
                                             object.y.meta.features)

      }

      rm(acts)
      gc()


      message(paste0("Finish merge ", k, " scores of object.y: ", i))

    }

  }


  return(object.x)

}


