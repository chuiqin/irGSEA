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
    stop("`method` should be one of the followling : AUCell, UCell, singscore, ssgsea, JASMINE, VAM, scSE, VISION, gficf, GSVA, zscore, plage, wmean, wsum, mdt, viper,  GSVApy, AddModuleScore, pagoda2, RRA.")
  }
  pvalue <- NULL
  if (method %in% names(object)[! names(object) == "RRA"]) {
    object[method] <- object[method] %>% purrr::map( ~.x %>% dplyr::rename(pvalue = p_val_adj))
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
                              show_row_names = FALSE, ...)

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

