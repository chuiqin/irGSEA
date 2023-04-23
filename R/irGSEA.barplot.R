
#' Stacked bar plot
#'
#' Easy to show analysis results by stacked bar plot
#'
#' @param object A list after perform \code{\link{irGSEA.integrate}}
#' @param method A vector. It should be one or more of the followling : AUCell,
#' UCell, singscore, ssgsea or RRA. Default all methods when it is set to NULL.
#' @param significance.color A vector. Default "c("#D0DFE6FF","#f87669")" when
#' it is set to NULL.
#' @param color.cluster A vector. Default "ggsci::pal_igv()(the number of colnames
#' of enrichment score matrix)" when it is set to NULL.
#' @param color.method A vector. Default "ggsci::pal_igv()(the number of methods
#' )" when it is set to NULL.
#' @param cluster.levels A vector equal to the number of clusters.
#'
#'
#' @return stacked bar plot
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
#' # Integrated analysis
#' result.dge <- irGSEA.integrate(object = pbmc3k.final,
#' group.by = "seurat_annotations", metadata = NULL, col.name = NULL,
#' method = c("AUCell","UCell","singscore","ssgsea"))
#'
#' irGSEA.barplot.plot1 <- irGSEA.barplot(object = result.dge)
#'
#' }
#'
irGSEA.barplot <- function(object = NULL, method = NULL,
                           significance.color = NULL,
                           color.cluster = NULL, color.method = NULL,
                           cluster.levels = NULL){
  # pretreatment
  if (! purrr::is_list(object)) {
    stop("object should be a list.")
  }
  if (! all(method %in% names(object))) {
    stop("`method` should be one or more of the followling : AUCell, UCell, singscore, ssgsea, JASMINE, VAM, scSE, VISION, gficf, GSVA, zscore, plage, wmean, wsum, mdt, viper,  GSVApy, AddModuleScore, pagoda2, RRA.")
  }
  cluster <- NULL
  direction <- NULL
  pvalue <- NULL
  cell <- NULL
  Name <- NULL
  geneset <- NULL
  proportion <- NULL
  anno.cluster <- NULL
  anno.method <- NULL

  if (purrr::is_null(method)){ method <- names(object) }

  for (i in method) {
    if (i != "RRA"){
      object[i] <- object[i] %>% purrr::map( ~.x %>% dplyr::rename(pvalue = p_val_adj))
    }
  }


  # matrix
  sig.genesets.barplot <- list()
  for (i in seq_along(names(object[method]))) {
    sig.genesets.barplot[[i]] <- object[method][[names(object[method])[i]]] %>%
      dplyr::mutate(cell = stringr::str_c(cluster, direction, sep = "_")) %>%
      dplyr::select(c("Name", "pvalue", "cell")) %>%
      dplyr::mutate(pvalue = dplyr::if_else(pvalue < 0.05, "significant","no significant")) %>%
      tidyr::spread(cell, pvalue, fill = "no significant") %>%
      tidyr::gather(cell, pvalue, -Name) %>%
      dplyr::mutate(direction = stringr::str_extract(cell, pattern = "up|down"),
                    cluster = stringr::str_remove(cell, pattern = "_up|_down"),
                    method = names(object[method])[i]) %>%
      dplyr::mutate(geneset = dplyr::if_else(pvalue == "no significant", "no significant", direction))
  }

  sig.genesets.barplot <- do.call(rbind, sig.genesets.barplot) %>%
    dplyr::group_by(cluster, method, geneset) %>%
    dplyr::summarise(proportion = dplyr::n()) %>%
    dplyr::mutate(cell = paste0(cluster, "_", method),
                  method = factor(method, levels = names(object))) %>%
    dplyr::arrange(cluster, method) %>%
    dplyr::mutate(cell = factor(cell, levels = unique(cell))) %>%
    dplyr::filter(method %in% tidyselect::all_of(method))

  # set levels of cluster
  if (! purrr::is_null(cluster.levels)) {
    sig.genesets.barplot <- sig.genesets.barplot %>%
      dplyr::mutate(cluster = factor(cluster, levels = cluster.levels)) %>%
      dplyr::arrange(cluster) %>%
      dplyr::mutate(cell = factor(cell, levels = unique(cell)))
  }

  # set color
  if (purrr::is_null(significance.color)) {
    significance.color <- c("#4575B4","#D0DFE6FF","#D73027")
  }
  # middle plot
  barplot.middle <- ggplot2::ggplot(sig.genesets.barplot , ggplot2::aes(cell, proportion, fill = geneset))+
    ggplot2::geom_bar(stat = "identity", position = "fill")+
    ggplot2::theme_bw()+
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Significance"))+
    ggplot2::scale_fill_manual(values = significance.color)+
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())+
    ggplot2::ylab("Proportion") +
    ggplot2::scale_y_continuous(expand = c(0,0))
  # set color
  significance.down.up.color <- c(significance.color[1],significance.color[3])

  # above color
  barplot.above <- sig.genesets.barplot %>%
    dplyr::filter(geneset != "no significant") %>%
    ggplot2::ggplot(ggplot2::aes(cell, proportion, fill = geneset)) +
    ggplot2::geom_bar(stat = "identity", colour = "white") +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::ylab("Count") +
    ggplot2::scale_fill_manual(values = significance.down.up.color) +
    ggplot2::geom_text(ggplot2::aes(label = proportion), color = "white",
                       position = ggplot2::position_stack(vjust = 0.5), size = 3)
  # set color
  if (purrr::is_null(color.cluster)) {
    color.cluster <- ggsci::pal_igv()(length(unique(sig.genesets.barplot$cluster)))
  }
  # below picture
  labels.cluster <- sig.genesets.barplot %>%
    dplyr::mutate(anno.cluster = "Cluster") %>%
    ggplot2::ggplot(ggplot2::aes(cell, y = anno.cluster, fill = cluster)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = color.cluster, name = "Cluster") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = NULL)
  # set color
  if (purrr::is_null(color.method)) {
    color.method <- ggsci::pal_igv()(length(unique(sig.genesets.barplot$method)))
  }
  # below picture
  labels.method <- sig.genesets.barplot %>%
    dplyr::mutate(anno.method = "Method") %>%
    ggplot2::ggplot(ggplot2::aes(cell, y = anno.method, fill = method)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = color.method, name = "Method") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = NULL)
  # combine all plots
  sig.genesets.bar.plot <- barplot.middle %>%
    aplot::insert_top(barplot.above, height = 1) %>%
    aplot::insert_bottom(labels.method, height = .1) %>%
    aplot::insert_bottom(labels.cluster, height = .1)
  sig.genesets.bar.plot <- ggplotify::as.ggplot(sig.genesets.bar.plot)
  return(sig.genesets.bar.plot)

}



