#' Bubble plot
#'
#' Easy to show analysis results by bubble plot
#'
#' @param object A list after perform \code{\link{irGSEA.integrate}}
#' @param method A character. It should be one of the followling : AUCell, UCell,
#' singscore, ssgsea or RRA. Default RRA.
#' @param top The top gene sets. Default 50.
#' @param show.geneset A vector including special gene sets. Default NULL.
#' @param cluster.color A vector. Default "ggsci::pal_igv()(the number of colnames
#' of enrichment score matrix)" when it is set to NULL.
#' @param direction.color A vector. Default "c("#4575B4","#D73027")" when it
#' is set to NULL.
#' @param significance.color A vector. Default "c("#D0DFE6FF","#f87669")" when
#' it is set to NULL.
#' @param cluster_rows Whether to make cluster on rows. Defaul True.
#' @param cluster.levels A vector equal to the number of clusters.
#'
#' @return bubble plot
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
#' irGSEA.bubble.plot1 <- irGSEA.bubble(object = result.dge,
#' method = "RRA", top = 50)
#' irGSEA.bubble.plot2 <- irGSEA.bubble(object = result.dge,
#' method = "ssgsea", top = 50)
#'
#' }
#'
#'
#'
irGSEA.bubble <- function(object = NULL, method = "RRA",
                          top = 50, show.geneset = NULL,
                          cluster.color = NULL, direction.color = NULL,
                          significance.color = NULL,
                          cluster_rows = T, cluster.levels = NULL){
  # pretreatment
  cluster <- NULL
  direction <- NULL
  cell <- NULL
  pvalue <- NULL
  Name <- NULL
  significance <- NULL
  anno.cluster <- NULL
  anno.direction <- NULL
  value <- NULL

  if (! purrr::is_list(object)) {
    stop("object should be a list.")
  }
  if ((! all(method %in% names(object))) | (length(method) > 1) | (purrr::is_null(method))) {
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea, JASMINE, VAM, scSE, VISION, gficf, GSVA, zscore, plage, wmean, wsum, mdt, viper,  GSVApy, AddModuleScore, pagoda2, RRA.")
  }
  pvalue <- NULL
  if (method %in% names(object)[! names(object) == "RRA"]) {
    object[method] <- object[method] %>% purrr::map( ~.x %>% dplyr::rename(pvalue = p_val_adj))
  }

  # matrix

  sig.genesets.bubble <- object[[method]] %>%
    dplyr::mutate(cell = stringr::str_c(cluster, direction, sep = "_")) %>%
    dplyr::select(c("Name", "pvalue", "cell")) %>%
    dplyr::mutate(pvalue = dplyr::case_when(  pvalue < 0.0001 ~ "< 0.0001",
                                              pvalue < 0.001 ~ "< 0.001",
                                              pvalue < 0.01 ~ "< 0.01",
                                              pvalue < 0.05 ~ "< 0.05",
                                              pvalue >= 0.05 ~ ">= 0.05",
                                              TRUE ~ NA_character_)) %>%
    dplyr::mutate(Name = factor(Name, levels = unique(Name))) %>%
    tidyr::spread(cell, pvalue, fill = ">= 0.05") %>%
    tibble::column_to_rownames(var = "Name")


  if (length(unique(object[[method]]$cluster)) != 0.5*ncol(sig.genesets.bubble)) {
    cell.name <- c(stringr::str_c(unique(object[[method]]$cluster), c("up"), sep = "_"),
                   stringr::str_c(unique(object[[method]]$cluster), c("down"), sep = "_"))
    cell.name <- cell.name[!cell.name %in% colnames(sig.genesets.bubble)]
    for (i in cell.name) {
      sig.genesets.bubble <- sig.genesets.bubble %>%
        dplyr::mutate(!!rlang::sym(i):= ">= 0.05")
      sig.genesets.bubble <- sig.genesets.bubble[, sort(colnames(sig.genesets.bubble))]
    }

  }

  # top rows or custom genesets
  if (purrr::is_null(show.geneset)) {
    sig.genesets.bubble <- sig.genesets.bubble %>% dplyr::slice_head(n = top)
  }else{
    sig.genesets.bubble <- sig.genesets.bubble[rownames(sig.genesets.bubble) %in% show.geneset, ]
    sig.genesets.bubble <- sig.genesets.bubble[intersect(show.geneset, rownames(sig.genesets.bubble)), ]

    if (purrr::is_null(sig.genesets.bubble)) {
      stop("All genesets of `show.geneset` are not in the `method`.")
    }
    if (! all(show.geneset %in% rownames(sig.genesets.bubble))) {
      a <- show.geneset[! show.geneset %in% rownames(sig.genesets.bubble)]
      message(paste0("Some genesets of `show.geneset` are not in the `method` : ",a))
    }
  }

  # continue to edit matrix
  sig.genesets.bubble <- sig.genesets.bubble %>%
    tibble::rownames_to_column(var = "Name") %>%
    tidyr::gather(cell, pvalue, -Name) %>%
    dplyr::mutate(direction = stringr::str_extract(cell, pattern = "up|down"),
                  anno.direction = "Direction",
                  cluster = stringr::str_remove(cell, pattern = "_up|_down"),
                  anno.cluster = "Cluster",
                  pvalue = factor(pvalue, levels = rev(levels(factor(pvalue)))),
                  significance = dplyr::if_else(pvalue == ">= 0.05", "no significant", "significant")) %>%
    dplyr::mutate(value = dplyr::if_else(significance == "no significant", 0, 1))

  # set levels of cluster
  if (! purrr::is_null(cluster.levels)) {
    sig.genesets.bubble <- sig.genesets.bubble %>%
      dplyr::mutate(cluster = factor(cluster, levels = cluster.levels)) %>%
      dplyr::arrange(cluster) %>%
      dplyr::mutate(cell = factor(cell, levels = unique(cell)))
  }



  # set color
  if (purrr::is_null(cluster.color)) {
    cluster.color <- ggsci::pal_igv()(length(unique(sig.genesets.bubble$cluster)))
  }
  # above picture
  labels.cluster <- sig.genesets.bubble %>%
    ggplot2::ggplot(ggplot2::aes(cell, y = anno.cluster, fill = cluster)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = cluster.color, name = "Cluster") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = NULL)

  # set color
  if (purrr::is_null(direction.color)) {
    direction.color <- c("#4575B4","#D73027")
  }
  # above picture
  labels.direction <- sig.genesets.bubble %>%
    ggplot2::ggplot(ggplot2::aes(cell, y = anno.direction, fill = direction)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = direction.color, name = "Direction") +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = NULL)

  # set color
  if (purrr::is_null(significance.color)) {
    significance.color <- c("#D0DFE6FF","#f87669")
  }
  # middle picture
  bubble.plot <- ggplot2::ggplot(sig.genesets.bubble, ggplot2::aes(x = cell, y = Name))+
    ggplot2::geom_point(ggplot2::aes(size = pvalue, color = significance))+
    ggplot2::scale_color_manual(values = significance.color, name = method)+
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 8))+
    ggplot2::labs(x=NULL,y=NULL)+
    ggplot2::guides(size = ggplot2::guide_legend(title = "P Value")) +
    ggplot2::scale_size_manual(values = c(1,1.5,2,2.5,3))

  # bulid ggtree matrix
  sig.genesets.bubble.matrix <- sig.genesets.bubble %>%
    dplyr::select(c(Name, cell, value)) %>%
    dplyr::mutate(Name = factor(Name, levels = unique(Name))) %>%
    tidyr::spread(cell, value) %>%
    tibble::column_to_rownames(var = "Name")

  # row tree
  phr <- ggtree::ggtree(stats::hclust(stats::dist(sig.genesets.bubble.matrix)))
  # combine tree
  if (cluster_rows == T) {
    sig.genesets.bubble.plot <- bubble.plot %>%
      aplot::insert_left(phr, width=.1) %>%
      aplot::insert_top(labels.direction, height = .05) %>%
      aplot::insert_top(labels.cluster, height = .05)
  }else{
    sig.genesets.bubble.plot <- bubble.plot %>%
      aplot::insert_top(labels.direction, height = .05) %>%
      aplot::insert_top(labels.cluster, height = .05)
  }
  sig.genesets.bubble.plot <- ggplotify::as.ggplot(sig.genesets.bubble.plot)
  return(sig.genesets.bubble.plot)

}

