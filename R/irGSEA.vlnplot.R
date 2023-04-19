#' Half vlnplot
#'
#' Easy to show the data distribution by half vlnplot
#'
#' @param object A Seurat after perform \code{\link{irGSEA.score}}
#' @param method A character. It should be one or more of the followling : AUCell,
#' UCell, singscore, ssgsea, JASMINE, VAM, scSE, VISION, gficf, GSVA, zscore,
#' plage, wmean, wsum, mdt, viper,  GSVApy, AddModuleScore, pagoda2.
#' @param show.geneset A character. It should be one of the rownames of
#' enrichment score matrix.
#' @param group.by  Default ident when it is set to NULL. You can specify other
#' column of metadata.
#' @param color.cluster A vector. Default "ggsci::pal_igv()(the number of colnames
#' of enrichment score matrix)" when it is set to NULL.
#' @param cluster.levels A vector equal to the number of clusters.
#'
#' @return vlnplot
#' @export
#'
#' @examples
#' \dontrun{
#' # load PBMC dataset by R package SeuratData
#' library(Seurat)
#' library(SeuratData)
#' # download 3k PBMCs from 10X Genomics
#' InstallData("pbmc3k")
#' data("pbmc3k.final")
#' pbmc3k.final <- SeuratObject::UpdateSeuratObject(pbmc3k.final)
#'
#' # Seurat object
#' pbmc3k.final <- irGSEA.score(object = pbmc3k.final, assay = "RNA",
#' slot = "data", msigdb = T, species = "Homo sapiens",
#' category = "H", geneid = "symbol",
#' method = c("AUCell", "UCell", "singscore", "ssgsea"), kcdf = 'Gaussian')
#'
#' irGSEA.vlnplot.plot1 <- irGSEA.vlnplot(object = pbmc3k.final,
#' method = c("AUCell", "UCell", "singscore", "ssgsea"),
#' show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
#' irGSEA.vlnplot.plot2 <- irGSEA.vlnplot(object = pbmc3k.final,
#' method = c("AUCell", "UCell", "singscore", "ssgsea"),
#' show.geneset = "HALLMARK-IL6-JAK-STAT3-SIGNALING")
#'
#'
#' }
#'
irGSEA.vlnplot <- function(object = NULL, method = NULL,
                               show.geneset = NULL, group.by = NULL,
                               color.cluster = NULL, cluster.levels = NULL){
  # pretreatment
  ident <- NULL
  geneset <- NULL
  if ((! all(method %in% Seurat::Assays(object))) |(purrr::is_null(method))) {
    stop("`method` should be one or more of the followling : AUCell, UCell, singscore, ssgsea, JASMINE, VAM, scSE, VISION, gficf, GSVA, zscore, plage, wmean, wsum, mdt, viper,  GSVApy, AddModuleScore, pagoda2.")
  }

  # group
  if (purrr::is_null(group.by)) {
    anno.ident <- SeuratObject::Idents(object)
  }else{
    object <- SeuratObject::SetIdent(object, value = group.by)
    anno.ident <- SeuratObject::Idents(object)
  }
  # factors are sorted alphabetically
  anno.ident <- as.factor(as.character(anno.ident))
  # set levels of cluster
  if (! purrr::is_null(cluster.levels)) {
    anno.ident <- factor(anno.ident, levels = cluster.levels)
  }
  SeuratObject::Idents(object) <- anno.ident

  # set colors
  if (purrr::is_null(color.cluster)) {
    color.cluster <- ggsci::pal_igv()(length(levels(object)))
  }

  # geneset
  if (purrr::is_null(show.geneset)) {
    stop("`show.geneset` can not be empty.")
  }

  # plot
  for (i in method) {
    object@meta.data <- object@meta.data %>%
      dplyr::mutate(!!rlang::sym(i):= as.numeric(object[[i]]@counts[show.geneset,]))
  }
  scores.vlnplot <- Seurat::VlnPlot(object = object,
                                    assay = "RNA",
                                    combine = T,
                                    stack = T,
                                    flip = T,
                                    group.by = group.by,
                                    cols = color.cluster,
                                    features = method,
                                    pt.size = 0,
                                    fill.by = "ident")+
    ggplot2::ggtitle(show.geneset)+
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   title = ggplot2::element_text(size =12),
                   axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 0.5, angle = 45))+
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Cluster"))


  return(scores.vlnplot)

}


