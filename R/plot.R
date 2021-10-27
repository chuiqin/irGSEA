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
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea, RRA.")
  }
  if (method %in% c("AUCell", "UCell", "singscore", "ssgsea")) {
    object[1:4] <- object[1:4] %>% purrr::map( ~.x %>% dplyr::rename(pvalue = p_val_adj))
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
    tidyr::spread(cell, pvalue, fill = ">= 0.05") %>%
    tibble::column_to_rownames(var = "Name")

  # top rows or custom genesets
  if (purrr::is_null(show.geneset)) {
    sig.genesets.bubble <- sig.genesets.bubble %>% dplyr::slice_head(n = top)
  }else{
    sig.genesets.bubble <- sig.genesets.bubble[rownames(sig.genesets.bubble) %in% show.geneset, ]
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
    cell.levels <- unlist(lapply(cluster.levels, function(x){paste(x, c("up","down"), sep = "_")}))
    sig.genesets.bubble <- sig.genesets.bubble %>%
      dplyr::mutate(cluster = factor(cluster, levels = cluster.levels)) %>%
      dplyr::mutate(cell = factor(cell, levels = cell.levels))
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
    ggplot2::scale_fill_manual(values = direction.color, name = method) +
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
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea, RRA.")
  }
  if (method %in% c("AUCell", "UCell", "singscore", "ssgsea")) {
    object[1:4] <- object[1:4] %>% purrr::map( ~.x %>% dplyr::rename(pvalue = p_val_adj))
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

  if (! purrr::is_null(cluster.levels)) {
    heatmap.levels <- unlist(lapply(cluster.levels, function(x){paste(x, c("up","down"), sep = "_")}))
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
  if (! purrr::is_null(cluster.levels)) {
    heatmap.levels <- unlist(lapply(cluster.levels, function(x){paste(x, c("up","down"), sep = "_")}))
    sig.genesets.heatmap <- sig.genesets.heatmap %>% dplyr::select(heatmap.levels)
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
  heatmap.legend <- ComplexHeatmap::packLegend(lgd1, lgd2, lgd3, lgd4,
                                               direction = "vertical",
                                               column_gap = grid::unit(1, "cm"))

  # plot
  heatmap.plot <- grid::grid.grabExpr(ComplexHeatmap::draw(heatmap.body,
                                                           annotation_legend_list = heatmap.legend)) %>% ggplotify::as.ggplot()
  return(heatmap.plot)
}




#' Upset plot
#'
#' Easy to show analysis results by upset plot
#'
#' @param object A list after perform \code{\link{irGSEA.integrate}}
#' @param method A character. It should be one of the followling : AUCell, UCell,
#' singscore, ssgsea or RRA. Default RRA.
#' @param upset.width Width of the whole upset plot. Default 13.
#' @param upset.height Height of the whole upset plot. Default 7.
#' @param title.size The fointsize of rownames. Default 10.
#' @param text.size The fointsize of rownames. Default 9.
#' @param cluster.color A vector. Default "ggsci::pal_igv()(the number of colnames
#' of enrichment score matrix)" when it is set to NULL.
#' @param bar.color A character. Default "black" when it is set to NULL.
#' @param cluster.levels A vector equal to the number of clusters.
#' @param mode A character. It should be one of the followling : distinct,
#' intersect, or union. Default distinct. It represents the mode for forming
#' the combination set. see Mode section \url{https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html}
#' for details.
#' @param set.size The minimal combination set size. Default 1.
#' @param set.degree A vector. Show all combination sets when it set to NULL.
#' It would show different combination set when it is set to different number.
#' For example, it only show the interactions among two cluster or three cluster
#' when it's set to 2 or 3.
#' @param table.generate Deault FALSE. It will output a list including all
#' combination sets and their gene sets when it set to TRUE.
#' @param ... More parameters pass to \code{\link[ComplexHeatmap]{UpSet}}
#'
#' @return upset plot or list
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
#' irGSEA.upset.plot1 <- irGSEA.upset(object = result.dge, method = "RRA")
#' irGSEA.upset.plot2 <- irGSEA.upset(object = result.dge, method = "ssgsea")
#'
#' }
#'
#'

irGSEA.upset <- function(object = NULL, method = "RRA",
                         upset.width = 13, upset.height = 7,
                         title.size = 10, text.size = 9,
                         cluster.color = NULL, bar.color = "black",
                         cluster.levels = NULL, mode = "distinct",
                         set.size = 1, set.degree = NULL, table.generate = F,
                         ...){
  # pretreatment
  if (! purrr::is_list(object)) {
    stop("object should be a list.")
  }
  if ((! all(method %in% names(object))) | (length(method) > 1) | (purrr::is_null(method))) {
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea, RRA.")
  }
  pvalue <- NULL
  if (method %in% c("AUCell", "UCell", "singscore", "ssgsea")) {
    object[1:4] <- object[1:4] %>% purrr::map( ~.x %>% dplyr::rename(pvalue = p_val_adj))
  }
  # matrix
  cluster <- NULL
  sig.genesets.upset <- object[[method]] %>%
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::select(c("Name","cluster"))
  sig.genesets.upset.names <- levels(factor(as.character(sig.genesets.upset$cluster)))
  if (purrr::is_null(cluster.color)) {
    cluster.color <- ggsci::pal_igv()(length(sig.genesets.upset.names))
  }

  sig.genesets.upset <- sig.genesets.upset %>%
    dplyr::group_split(cluster,.keep = F) %>%
    purrr::set_names(sig.genesets.upset.names)
  sig.genesets.upset <- lapply(sig.genesets.upset, function(x){x <- x$Name})
  # set levels of cluster
  if (! purrr::is_null(cluster.levels)) {
    sig.genesets.upset <- sig.genesets.upset[cluster.levels]
  }

  # matrix
  m <- ComplexHeatmap::make_comb_mat(sig.genesets.upset, mode = mode)
  # degree
  if (purrr::is_null(set.degree)) {
    m <- m[ComplexHeatmap::comb_degree(m) > 0]
  }else{
    m <- m[ComplexHeatmap::comb_degree(m) %in% c(set.degree)]
  }
  # set size
  m <- m[ComplexHeatmap::comb_size(m) >= set.size]

  # set levels of cluster
  ss <- ComplexHeatmap::set_size(m)
  if (! purrr::is_null(cluster.levels)) {
    cluster.order <- cluster.levels
  }else{
    cluster.order <- order(ss)
  }

  # table
  if (table.generate == T) {
    comb.list <- lapply(ComplexHeatmap::comb_name(m), function(i){ComplexHeatmap::extract_comb(m, i)})
    names(comb.list) <- ComplexHeatmap::comb_name(m, readable = TRUE)
    return(comb.list)
  }
  # plot
  cs <- ComplexHeatmap::comb_size(m)
  ht <- ComplexHeatmap::UpSet(m,
                              bg_col = "#D0DFE6FF",
                              heatmap_width = grid::unit(upset.width, "cm"),
                              heatmap_height = grid::unit(upset.height, "cm"),
                              column_title = method,
                              column_title_gp = grid::gpar(fontsize = title.size),
                              set_order = cluster.order,
                              comb_order = order(ComplexHeatmap::comb_degree(m), -cs),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(
                                "Intersections" = ComplexHeatmap::anno_barplot(cs,
                                                                               ylim = c(0, max(cs)*1.1),
                                                                               border = FALSE,
                                                                               gp = grid::gpar(color = bar.color),
                                                                               height = grid::unit(2, "cm")
                                ),
                                annotation_name_side = "left",
                                annotation_name_rot = 90,
                                annotation_name_gp = grid::gpar(fontsize = text.size)),
                              left_annotation = ComplexHeatmap::rowAnnotation(
                                "Counts" = ComplexHeatmap::anno_barplot(-ss, baseline = 0,
                                                                        axis_param = list(at = c(0, -500, -1000, -1500),
                                                                                          labels = c(0, 500, 1000, 1500),
                                                                                          labels_rot = 0),
                                                                        border = FALSE,
                                                                        gp = grid::gpar(fill = cluster.color),
                                                                        width = grid::unit(1.5, "cm")), annotation_name_gp = grid::gpar(fontsize = text.size),
                                set_name = ComplexHeatmap::anno_text(ComplexHeatmap::set_name(m),
                                                                     location = 0.5,
                                                                     just = "center",
                                                                     gp = grid::gpar(fontsize = text.size),
                                                                     width = ComplexHeatmap::max_text_width(text = ComplexHeatmap::set_name(m),
                                                                                                            gp = grid::gpar(fontsize = text.size)) + grid::unit(1, "mm"))),
                              right_annotation = NULL,
                              show_row_names = FALSE)

  plot_upset <- function(){
    ht <- ComplexHeatmap::draw(ht)
    od <- ComplexHeatmap::column_order(ht)
    ComplexHeatmap::decorate_annotation("Intersections", {
      grid::grid.text(cs[od], x = seq_along(cs), y = grid::unit(cs[od], "native") + grid::unit(2, "pt"),
                      default.units = "native", just = c("left", "bottom"),
                      gp = grid::gpar(fontsize = 6, col = "black"), rot = 45)
    })

  }

  upset.plot <- grid::grid.grabExpr(plot_upset()) %>%
    ggplotify::as.ggplot()
  return(upset.plot)
}




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
    stop("`method` should be one or more of the followling : AUCell, UCell, singscore, ssgsea, RRA.")
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
  object[1:4] <- object[1:4] %>% purrr::map( ~.x %>% dplyr::rename(pvalue = p_val_adj))
  if (purrr::is_null(method)){ method <- names(object) }
  # matrix
  sig.genesets.barplot <- list()
  for (i in seq_along(names(object))) {
    sig.genesets.barplot[[i]] <- object[[names(object)[i]]] %>%
      dplyr::mutate(cell = stringr::str_c(cluster, direction, sep = "_")) %>%
      dplyr::select(c("Name", "pvalue", "cell")) %>%
      dplyr::mutate(pvalue = dplyr::if_else(pvalue < 0.05, "significant","no significant")) %>%
      tidyr::spread(cell, pvalue, fill = "no significant") %>%
      tidyr::gather(cell, pvalue, -Name) %>%
      dplyr::mutate(direction = stringr::str_extract(cell, pattern = "up|down"),
                    cluster = stringr::str_remove(cell, pattern = "_up|_down"),
                    method = names(object)[i]) %>%
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
    cell.levels <- unlist(lapply(cluster.levels, function(x){paste(x, method, sep = "_")}))
    sig.genesets.barplot <- sig.genesets.barplot %>%
      dplyr::mutate(cluster = factor(cluster, levels = cluster.levels)) %>%
      dplyr::mutate(cell = factor(cell, levels = cell.levels))
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
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea.")
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
                                                ...) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
                   axis.title = ggplot2::element_text(size = 10))

  return(scores.scatterplot)

}



#' Half vlnplot
#'
#' Easy to show the data distribution by half vlnplot
#'
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
#'
#' @return half vlnplot
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
#' irGSEA.halfvlnplot.plot1 <- irGSEA.halfvlnplot(object = pbmc3k.final,
#' method = "UCell", show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
#' irGSEA.halfvlnplot.plot2 <- irGSEA.halfvlnplot(object = pbmc3k.final,
#' method = "ssgsea", show.geneset = "HALLMARK-IL6-JAK-STAT3-SIGNALING")
#'
#'
#' }
#'
irGSEA.halfvlnplot <- function(object = NULL, method = NULL,
                               show.geneset = NULL, group.by = NULL,
                               color.cluster = NULL, cluster.levels = NULL){
  # pretreatment
  ident <- NULL
  geneset <- NULL
  if ((! all(method %in% Seurat::Assays(object))) | (length(method) > 1) | (purrr::is_null(method))) {
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea.")
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
  scores.vlnplot <- Seurat::VlnPlot(object = object,
                                    assay = method,
                                    slot = "scale.data",
                                    group.by = group.by,
                                    cols = color.cluster,
                                    features = custom.geneset,
                                    pt.size = 0,
                                    fill.by = "ident") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())+
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Cluster"))

  scores.vlnplot <- scores.vlnplot$data
  scores.vlnplot <- scores.vlnplot %>%
    dplyr::rename(c("geneset" = tidyselect::all_of(custom.geneset))) %>%
    ggplot2::ggplot(ggplot2::aes(x = ident, y = geneset, fill = ident), colour = "white") +
    gghalves::geom_half_boxplot(side = "r", outlier.color = NA,errorbar.draw = TRUE)+
    gghalves::geom_half_violin() +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_manual(values = color.cluster)+
    ggplot2::ylab(paste0(method, " scores")) +
    ggplot2::xlab("") +
    ggplot2::ggtitle(custom.geneset) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   title = ggplot2::element_text(size =12),
                   axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 0.5, angle = 45))+
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Cluster"))

  return(scores.vlnplot)

}


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
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea.")
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
                                        fill.by = "ident", ...) +
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
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea.")
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
                                           slot = "scale.data",
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


