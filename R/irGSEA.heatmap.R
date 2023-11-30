#' Heatmap plot
#'
#' Easy to show analysis results by heatmap plot
#'
#' @param object A list after perform \code{\link{irGSEA.integrate}}
#' @param method A character. It should be one of the followling : AUCell, UCell,
#' singscore, ssgsea or RRA. Default RRA.
#' @param top The top gene sets. Default 50.
#' @param show.geneset A vector including special gene sets. Default NULL.
#' @param cluster_rows Whether to make cluster on rows. Defaul True.
#' @param significance.color A vector. Default "c("#D0DFE6FF","#f87669")" when
#' it is set to NULL.
#' @param cluster.color A vector. Default "ggsci::pal_igv()(the number of colnames
#' of enrichment score matrix)" when it is set to NULL.
#' @param direction.color A vector. Default "c("#4575B4","#D73027")" when it
#' is set to NULL.
#' @param rowname.fointsize The fointsize of rownames. Default 7.
#' @param heatmap.width Width of the whole heatmap (including heatmap
#' components), default 17.
#' @param heatmap.heigh Height of the whole heatmap (including heatmap
#' components), default 13.
#' @param cluster.levels A vector equal to the number of clusters.
#' @param ... More parameters pass to \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @return heatmap plot
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
#' irGSEA.heatmap.plot1 <- irGSEA.heatmap(object = result.dge, method = "RRA",
#' top = 50, show.geneset = NULL)
#'
#' irGSEA.heatmap.plot2 <- irGSEA.heatmap(object = result.dge, method = "ssgsea",
#' top = 50, show.geneset = NULL)
#'
#' }
#'
#'
#'
irGSEA.heatmap <- function(object = NULL, method = "RRA",
                           top = 50, show.geneset = NULL,
                           cluster_rows = T,
                           significance.color = NULL, cluster.color = NULL,
                           direction.color = NULL, rowname.fointsize = 7,
                           heatmap.width = 17, heatmap.heigh = 13,
                           cluster.levels = NULL, ...
){
  # pretreatment
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
  cluster <- NULL
  direction <- NULL
  pvalue <- NULL
  cell <- NULL
  value <- NULL
  Name <- NULL
  sig.genesets.heatmap <- object[[method]] %>%
    dplyr::mutate(cell = stringr::str_c(cluster, direction, sep = "_")) %>%
    dplyr::select(c("Name", "pvalue", "cell")) %>%
    dplyr::mutate(pvalue = dplyr::if_else(pvalue < 0.05, "significant","no significant")) %>%
    tidyr::spread(cell, pvalue, fill = "no significant") %>%
    tibble::column_to_rownames(var = "Name")

  if (length(unique(object[[method]]$cluster)) != 0.5*ncol(sig.genesets.heatmap)) {
    cell.name <- c(stringr::str_c(unique(object[[method]]$cluster), c("up"), sep = "_"),
                   stringr::str_c(unique(object[[method]]$cluster), c("down"), sep = "_"))
    cell.name <- cell.name[!cell.name %in% colnames(sig.genesets.heatmap)]
    for (i in cell.name) {
      sig.genesets.heatmap <- sig.genesets.heatmap %>%
        dplyr::mutate(!!rlang::sym(i):= "no significant")
      sig.genesets.heatmap <- sig.genesets.heatmap[, sort(colnames(sig.genesets.heatmap))]
    }

  }

  sig.genesets.heatmap.text <- object[[method]] %>%
    dplyr::mutate(cell = stringr::str_c(cluster, direction, sep = "_")) %>%
    dplyr::select(c("Name", "pvalue", "cell")) %>%
    dplyr::mutate(pvalue = dplyr::case_when(  pvalue < 0.0001 ~ "****",
                                              pvalue < 0.001 ~ "***",
                                              pvalue < 0.01 ~ "**",
                                              pvalue < 0.05 ~ "*",
                                              pvalue >= 0.05 ~ " ",
                                              TRUE ~ NA_character_)) %>%
    tidyr::spread(cell, pvalue, fill = " ") %>%
    tibble::column_to_rownames(var = "Name")

  if (length(unique(object[[method]]$cluster)) != 0.5*ncol(sig.genesets.heatmap.text)) {
    cell.name <- c(stringr::str_c(unique(object[[method]]$cluster), c("up"), sep = "_"),
                   stringr::str_c(unique(object[[method]]$cluster), c("down"), sep = "_"))
    cell.name <- cell.name[!cell.name %in% colnames(sig.genesets.heatmap.text)]
    for (i in cell.name) {
      sig.genesets.heatmap.text <- sig.genesets.heatmap.text %>%
        dplyr::mutate(!!rlang::sym(i):= " ")
      sig.genesets.heatmap.text <- sig.genesets.heatmap.text[, sort(colnames(sig.genesets.heatmap.text))]
    }

  }


  # set levels
  if (! purrr::is_null(cluster.levels)) {
    cluster.direction <- NULL
    heatmap.levels <- data.frame(cluster.direction = colnames(sig.genesets.heatmap)) %>%
      dplyr::mutate(cluster = stringr::str_remove(cluster.direction, pattern = "_up|_down")) %>%
      dplyr::arrange(factor(cluster, levels = cluster.levels)) %>%
      dplyr::pull(cluster.direction)
    sig.genesets.heatmap <- sig.genesets.heatmap %>% dplyr::select(heatmap.levels)
    sig.genesets.heatmap.text <- sig.genesets.heatmap.text %>% dplyr::select(heatmap.levels)
  }

  # top rows
  if (purrr::is_null(show.geneset)) {
    sig.genesets.heatmap <- sig.genesets.heatmap %>% dplyr::slice_head(n = top)
    sig.genesets.heatmap.text <- sig.genesets.heatmap.text %>% dplyr::slice_head(n = top)
  }else{
    sig.genesets.heatmap <- sig.genesets.heatmap[rownames(sig.genesets.heatmap) %in% show.geneset, ]
    sig.genesets.heatmap.text <- sig.genesets.heatmap.text[rownames(sig.genesets.heatmap.text) %in% show.geneset, ]
    if (purrr::is_null(sig.genesets.heatmap)) {
      stop("All genesets of `show.geneset` are not in the `method`.")
    }
    if (! all(show.geneset %in% rownames(sig.genesets.heatmap))){
      a <- show.geneset[! show.geneset %in% rownames(sig.genesets.heatmap)]
      message(paste0("Some genesets of `show.geneset` are not in the `method` : ",a))
    }
  }

  # top annotation
  sig.genesets.heatmap.cluster <- stringr::str_remove(colnames(sig.genesets.heatmap), pattern = "_up|_down")
  sig.genesets.heatmap.direction <- stringr::str_extract(colnames(sig.genesets.heatmap), pattern = "up|down")
  if (purrr::is_null(cluster.color)) {
    cluster.color <- ggsci::pal_igv()(length(unique(sig.genesets.heatmap.cluster)))
  }

  if (purrr::is_null(direction.color)) {
    direction.color <- c("#4575B4","#D73027")
  }
  heatmap.top.anno <- ComplexHeatmap::HeatmapAnnotation(Cluster = sig.genesets.heatmap.cluster,
                                                        Direction = sig.genesets.heatmap.direction,
                                                        show_legend = F,
                                                        show_annotation_name = T,
                                                        gap = grid::unit(1, "mm"),
                                                        annotation_name_gp= grid::gpar(fontsize = 8),
                                                        col = list(Cluster = c(structure(cluster.color,
                                                                                         names = unique(sig.genesets.heatmap.cluster))),
                                                                   Direction = c(structure(direction.color,
                                                                                           names = c("down","up")))))
  # heatmap body
  if (purrr::is_null(significance.color)) {
    significance.color <- structure(c("#D0DFE6FF","#f87669"), names = c(0,1))
  }

  sig.genesets.heatmap <- sig.genesets.heatmap %>%
    tibble::rownames_to_column(var = "Name") %>%
    tidyr::gather(cell, value, -Name) %>%
    dplyr::mutate(value = dplyr::if_else(value=="no significant", 0, 1)) %>%
    tidyr::spread(cell, value) %>%
    tibble::column_to_rownames(var = "Name")

  # set levels
  if (! purrr::is_null(cluster.levels)) {
    cluster.direction <- NULL
    heatmap.levels <- data.frame(cluster.direction = colnames(sig.genesets.heatmap)) %>%
      dplyr::mutate(cluster = stringr::str_remove(cluster.direction, pattern = "_up|_down")) %>%
      dplyr::arrange(factor(cluster, levels = cluster.levels)) %>%
      dplyr::pull(cluster.direction)
    sig.genesets.heatmap <- sig.genesets.heatmap %>% dplyr::select(heatmap.levels)
    sig.genesets.heatmap.text <- sig.genesets.heatmap.text %>% dplyr::select(heatmap.levels)
  }

  sig.genesets.heatmap <- as.matrix(sig.genesets.heatmap)

  heatmap.body <- ComplexHeatmap::Heatmap(sig.genesets.heatmap,
                                          heatmap_width = grid::unit(heatmap.width, "cm"),
                                          heatmap_height = grid::unit(heatmap.heigh, "cm"),
                                          name = method,
                                          col = significance.color,
                                          cluster_rows = cluster_rows,
                                          row_names_max_width = ComplexHeatmap::max_text_width(
                                            rownames(sig.genesets.heatmap),
                                            gp = grid::gpar(fontsize = rowname.fointsize)
                                          ),
                                          cluster_columns = F,
                                          top_annotation = heatmap.top.anno,
                                          color_space = "RGB",
                                          show_column_names = F,
                                          row_names_side="right",
                                          row_names_gp = grid::gpar(fontsize = rowname.fointsize),
                                          rect_gp = grid::gpar(col = "white", lwd = 2),
                                          show_heatmap_legend = F,
                                          cell_fun = function(j, i, x, y, width, height, fill){
                                            grid::grid.text(sig.genesets.heatmap.text[i, j], x, y, gp = grid::gpar(fontsize = 10))
                                          }, ...)

  # legend
  # represent sigficant
  lgd1 <- ComplexHeatmap::Legend(labels = c("no significant","significant"),
                                 title = method, legend_gp = grid::gpar(fill = significance.color))
  # represent p value
  lgd2 <- ComplexHeatmap::Legend(pch = c("*","**","***","****"),
                                 type = "points", labels = c("< 0.05","< 0.01","< 0.001","< 0.0001"),
                                 title = "P Value")
  # represent cluster
  lgd3 <- ComplexHeatmap::Legend(labels = unique(sig.genesets.heatmap.cluster),
                                 title = "Cluster",
                                 legend_gp = grid::gpar(fill = cluster.color))
  # represent direction
  lgd4 <- ComplexHeatmap::Legend(labels = c("down","up"),
                                 title = "Direction",
                                 legend_gp = grid::gpar(fill = direction.color),
                                 labels_gp = grid::gpar(fill = direction.color))

  # merge all legend
  heatmap.legend <- ComplexHeatmap::packLegend(lgd3, lgd4, lgd1, lgd2,
                                               direction = "vertical",
                                               column_gap = grid::unit(1, "cm"))

  # plot
  heatmap.plot <- grid::grid.grabExpr(ComplexHeatmap::draw(heatmap.body,
                                                           annotation_legend_list = heatmap.legend)) %>% ggplotify::as.ggplot()
  return(heatmap.plot)
}
