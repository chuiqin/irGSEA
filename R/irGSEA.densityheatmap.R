#' Density heatmap
#'
#' Easy to show the data distribution by density heatmap
#'
#' @param object A Seurat after perform \code{\link{irGSEA.score}}
#' @param method A character. It should be one of the followling : AUCell,
#' UCell, singscore, ssgsea.
#' @param show.geneset A character. It should be one of the rownames of
#' enrichment score matrix.
#' @param group.by Default ident when it is set to NULL. You can specify other
#' column of metadata.
#' @param heatmap_width Width of the whole heatmap (including heatmap
#' components), default 12.
#' @param heatmap_height Heigh of the whole heatmap (including heatmap
#' components), default 12.
#' @param cluster.levels A vector equal to the number of clusters.
#' @param ... More parameters pass to \code{\link[ComplexHeatmap]{densityHeatmap}}
#'
#' @return density heatmap
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
#' irGSEA.densityheatmap.plot1 <- irGSEA.densityheatmap(object = pbmc3k.final,
#' method = "UCell", show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
#' irGSEA.densityheatmap.plot2 <- irGSEA.densityheatmap(object = pbmc3k.final,
#' method = "ssgsea", show.geneset = "HALLMARK-IL6-JAK-STAT3-SIGNALING")
#'
#' }
#'
irGSEA.densityheatmap <- function(object = NULL, method = NULL,
                                  show.geneset = NULL, group.by = NULL,
                                  heatmap_width = 12, heatmap_height = 12,
                                  cluster.levels = NULL, ...){
  # pretreatment
  ident <- NULL
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
  scores.densityheatmap <- Seurat::VlnPlot(object = object,
                                           assay = method,
                                           slot = "data",
                                           group.by = group.by,
                                           features = custom.geneset,
                                           pt.size = 0,
                                           fill.by = "ident")
  scores.densityheatmap <- scores.densityheatmap$data %>%
    dplyr::group_split(ident, .keep = F) %>%
    purrr::set_names(levels(scores.densityheatmap$data$ident)) %>%
    purrr::map(~ .x %>% dplyr::pull(.))


  scores.densityheatmap <- ComplexHeatmap::densityHeatmap(scores.densityheatmap,
                                                          title = paste0(method,"'s Distribution"),
                                                          ylab = custom.geneset,
                                                          title_gp = grid::gpar(fontsize = 9),
                                                          ylab_gp = grid::gpar(fontsize = 8),
                                                          tick_label_gp = grid::gpar(fontsize = 9),
                                                          column_names_gp = grid::gpar(fontsize = 9),
                                                          column_names_rot = 45,
                                                          show_heatmap_legend = T,
                                                          heatmap_width = grid::unit(heatmap_width, "cm"),
                                                          heatmap_height = grid::unit(heatmap_height, "cm"),
                                                          ...)
  scores.densityheatmap <- ggplotify::as.ggplot(
    grid::grid.grabExpr(ComplexHeatmap::draw(scores.densityheatmap)))

  return(scores.densityheatmap)

}


