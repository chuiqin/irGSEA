#' Ridge plot
#' Easy to show the data distribution by ridge plot
#' @param object A Seurat after perform \code{\link{irGSEA.score}}
#' @param method A character. It should be one of the followling : AUCell,
#' UCell, singscore, ssgsea.
#' @param show.geneset A character. It should be one of the rownames of
#' enrichment score matrix.
#' @param group.by  Default ident when it is set to NULL. You can specify other
#' column of metadata.
#' @param color.cluster A vector. Default "ggsci::pal_igv()(the number of colnames
#' of enrichment score matrix)" when it is set to NULL.
#' @param cluster.levels A vector equal to the number of clusters.
#' @param ... More parameters pass to \code{\link[ggridges]{geom_density_ridges}}
#'
#' @return ridge plot
#' @export
#'
#' @examples
#'
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
#' irGSEA.ridgeplot.plot1 <- irGSEA.ridgeplot(object = pbmc3k.final,
#' method = "UCell", show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
#' irGSEA.ridgeplot.plot2 <- irGSEA.ridgeplot(object = pbmc3k.final,
#' method = "ssgsea", show.geneset = "HALLMARK-IL6-JAK-STAT3-SIGNALING")
#'
#' }
#'
irGSEA.ridgeplot <- function(object = NULL, method = NULL,
                             show.geneset = NULL, group.by = NULL,
                             color.cluster = NULL, cluster.levels = NULL, ...){
  # pretreatment
  ident <- NULL
  geneset <- NULL
  if ((! all(method %in% Seurat::Assays(object))) | (length(method) > 1) | (purrr::is_null(method))) {
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea, JASMINE, VAM, scSE, VISION, gficf, GSVA, zscore, plage, wmean, wsum, mdt, viper,  GSVApy, AddModuleScore, pagoda2.")
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
  }else{
    custom.geneset <- show.geneset[show.geneset %in% rownames(object[[method]])]
    if (purrr::is_null(custom.geneset)) {
      stop("All genesets of `show.geneset` are not in the `method`.")
    }
    if (! all(show.geneset %in% rownames(object[[method]]))) {
      a <- show.geneset[! show.geneset %in% rownames(object[[method]])]
      message(paste0("Following genesets of `show.geneset` are not in such `method` : ",a))
    }
  }

  # plot
  scores.ridgeplot <- Seurat::RidgePlot(object = object,
                                        assay = method,
                                        slot = "scale.data",
                                        group.by = group.by,
                                        cols = color.cluster,
                                        features = custom.geneset,
                                        fill.by = "ident") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.title.x = ggplot2::element_text(hjust = 0.5),
                   panel.grid = ggplot2::element_blank())+
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Cluster"))

  scores.ridgeplot <- scores.ridgeplot$data %>%
    dplyr::rename(c("geneset" = tidyselect::all_of(custom.geneset))) %>%
    ggplot2::ggplot(ggplot2::aes(x = geneset, y = ident, fill = ident)) +
    ggridges::geom_density_ridges(jittered_points=TRUE, scale = .95, rel_min_height = .01,
                                  point_shape = "|", point_size = 3, size = 0.25,
                                  position = ggridges::position_points_jitter(height = 0), ...) +
    ggplot2::scale_y_discrete(expand = c(.01, 0), name = "Cluster") +
    ggplot2::scale_x_continuous(expand = c(0, 0), name = paste0(method, " scores")) +
    ggplot2::scale_fill_manual(values = color.cluster) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      title="Cluster",
      override.aes = list(fill = color.cluster,  color = NA, point_color = NA))) +
    ggplot2::ggtitle(custom.geneset) +
    ggridges::theme_ridges(center = TRUE)+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                   axis.title = ggplot2::element_text(size = 12))


  return(scores.ridgeplot)

}

