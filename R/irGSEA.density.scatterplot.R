#' Density Scatter plot
#'
#' Easy to the data distribution by density scatter plot
#'
#' @param object A Seurat after perform \code{\link{irGSEA.score}}
#' @param method A character. It should be one of the followling : AUCell,
#' UCell, singscore, ssgsea.
#' @param show.geneset A character. It should be one of the rownames of
#' enrichment score matrix.
#' @param reduction A character. It can not be empty and should be calculated
#' in advance.
#' @param ... More parameters pass to \code{\link[Nebulosa]{plot_density}}
#'
#' @return density scatter plot
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
#' irGSEA.density.scatterplot1 <- irGSEA.density.scatterplot(object = pbmc3k.final,
#' method = "UCell", show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE",
#' reduction = "umap")
#' irGSEA.density.scatterplot2 <- irGSEA.density.scatterplot(object = pbmc3k.final,
#' method = "ssgsea", show.geneset = "HALLMARK-IL6-JAK-STAT3-SIGNALING",
#' reduction = "umap")
#'
#' }
#'
#'
irGSEA.density.scatterplot <- function(object = NULL, method = NULL,
                                       show.geneset = NULL, reduction = "umap",
                                       ...){
  # pretreatment
  if ((! all(method %in% Seurat::Assays(object))) | (length(method) > 1) | (purrr::is_null(method))) {
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea, JASMINE, VAM, scSE, VISION, gficf, GSVA, zscore, plage, wmean, wsum, mdt, viper,  GSVApy, AddModuleScore, pagoda2.")
  }

  if ((! all(reduction %in% Seurat::Reductions(object))) | (length(reduction) > 1) | (purrr::is_null(reduction))) {
    stop("`reductions` can not be empty and should be calculated in advance.")
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
  SeuratObject::DefaultAssay(object) <- method
  scores.scatterplot  <- Nebulosa::plot_density(object,
                                                features = custom.geneset,
                                                slot = "scale.data",
                                                reduction = reduction,
                                                method = "wkde",
                                                joint = T,
                                                ...) +
    ggplot2::ggtitle(paste0(method, ": ", custom.geneset))+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
                   axis.title = ggplot2::element_text(size = 10))

  return(scores.scatterplot)

}
