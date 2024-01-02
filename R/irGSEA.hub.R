#' Calculate the hub gene of the geneset
#'
#' Easy to calculate the hub gene of the geneset based on the correlation between
#' the geneset's score and the expression or rank of gene included in the geneset
#'
#' @param object A Seurat after perform \code{\link{irGSEA.score}}
#' @param assay Name of assay to calculate the correction. Default RNA.
#' @param slot Default data.
#' @param method A character.
#' @param show.geneset A character.
#' @param ncores Default 4.
#' @param type  expression or rank. Default rank. Calculate the correlation between
#' the geneset's score and the  rank of gene included in the geneset while the
#' type is rank. Calculate the correlation between the geneset's score and the
#' expression of gene included in the geneset while the type is expression.
#' @param maxRank Maximum number of genes to rank per cell; above this rank,
#' a given gene is considered as not expressed. Default 2000.
#' @param top The number of top-ranked positively correlated genes in each method
#' is displayed in a heatmap.
#' @param correlation.color A vector.
#' @param method.color A vector. Default "ggsci::pal_png()(the number of
#' methods)" when it is set to NULL.
#' @return list includes hub_result and hub_plot
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
#' hub.result <- irGSEA.hub(object = pbmc3k.final, assay = "RNA", slot = "data",
#' method = c("AUCell","UCell","singscore", "ssgsea"),
#' show.geneset = c("HALLMARK-INFLAMMATORY-RESPONSE", "HALLMARK-APOPTOSIS"),
#' ncores = 4, type = "rank", maxRank = 2000, top = 5,
#' correlation.color = c("#0073c2","white","#efc000"), method.color = NULL)
#'
#' head(hub.result$hub_result)
#' hub.result$hub_plot$`HALLMARK-APOPTOSIS`
#' hub.result$hub_plot$`HALLMARK-INFLAMMATORY-RESPONSE`
#'
#'
#' }
#'
irGSEA.hub <- function(object = NULL, assay = "RNA", slot = "data",
                       method = NULL, show.geneset = NULL,
                       ncores = 4, type = "rank", maxRank = 2000, top = 5,
                       correlation.color = c("#0073c2","white","#efc000"),
                       method.color = NULL){

  method.data <- lapply(method, function(i){
    message(i)
    # method
    tryCatch({

      # All the genesets are not on the assay
      if (all(show.geneset %in% rownames(object[[i]])) == F) {
        message(paste0("All the genesets: ", show.geneset, " are not on the assay: ",
                       i, "."))
        return(NULL)
      }

      show.geneset.before <- length(show.geneset)
      show.geneset.new <- show.geneset[show.geneset %in% rownames(object[[i]])]
      show.geneset.after <- length(show.geneset)

      # Some genesets are not on the assay
      if (show.geneset.before > show.geneset.after) {
        message(paste0("Some genesets: ", show.geneset[!show.geneset %in% show.geneset.new],
                       " are not on the assay: ", i, "."))
      }

      # geneset matrix
      geneset.data <- SeuratObject::GetAssayData(object[[i]], slot = "scale.data")[show.geneset.new, , drop=F]

      cor.geneset <- lapply(show.geneset.new, function(x){
        message(x)

        # the target gene of geneset
        if (class(object[[i]])[1] == "Assay5") {
          gene <- object[[i]]@meta.data[x, "target.gene"]
        }else{
          gene <- object[[i]]@meta.features[x, "target.gene"]
        }

        # if target gene of geneset is null
        if (gene == "") {
          message(paste0("Please check the column `target.gene` in meta.data or meta.features of",
                         " the geneset `", x, "` in the assay `", i, "`. It maybe empty."))
          return(NULL)
        }

        gene <- stringr::str_remove_all(gene, pattern = "\\+$|-$")
        gene <- stringr::str_split(gene, pattern = ", ")[[1]]

        # gene expression matrix or gene rank matrix
        if (type == "rank") {
          expression.data <- UCell::StoreRankings_UCell(matrix = SeuratObject::GetAssayData(object[[assay]], slot = slot),
                                                        maxRank = maxRank,
                                                        ncores = ncores)
          expression.data <- as.matrix(expression.data)[gene, , drop=F]
        }else{
          expression.data <- SeuratObject::GetAssayData(object[[assay]], slot = slot)[gene, , drop=F]
        }

        # calcute correlation
        cor.list <- lapply(gene, function(k){
          #message(k)
          test <- stats::cor.test(as.numeric(expression.data[k, ]),
                           as.numeric(geneset.data[x, ]),
                           type = "spearman")
          cor.data <- data.frame(method = i,
                                 geneset = x,
                                 gene = k,
                                 correlation = test$estimate,
                                 p.value = test$p.value)
          return(cor.data)
        })
        cor.list <- do.call(rbind, cor.list)
        return(cor.list)

      })
      cor.geneset <- do.call(rbind, cor.geneset)
      return(cor.geneset)

    }, error = function(e) {
      cat("Error: ", conditionMessage(e), "\n")
    })

  })
  method.data <- do.call(rbind, method.data)
  method.data$geneset <- factor(method.data$geneset, levels = show.geneset)
  method.data$method <- factor(method.data$method, levels = method)
  geneset <- NULL
  method.data2 <- method.data %>%
    dplyr::group_split(geneset, .keep = T) %>%
    purrr::set_names(levels(method.data$geneset))

  # draw
  method.data2 <- lapply(method.data2, function(method.data3){


    method <- NULL
    Module <- NULL
    correlation <- NULL
    p.value <- NULL
    gene <- NULL
    pvalue <- NULL
    Variable <- NULL
    Value <- NULL
    text <- NULL
    Method <- NULL

    method.data3 <- method.data3 %>%
      dplyr::group_by(method) %>%
      dplyr::arrange(dplyr::desc(correlation)) %>%
      dplyr::slice_head(n = top)

    method.data3 <- method.data3 %>%
      dplyr::mutate(pvalue = p.value, Value = correlation,
                    Variable = gene, Module = method) %>%
      dplyr::mutate(pvalue = dplyr::case_when(pvalue < 1e-04 ~ "****",
                                              pvalue < 0.001 ~ "***",
                                              pvalue < 0.01 ~ "**",
                                              pvalue <= 0.05 ~ "*",
                                              pvalue > 0.05 ~ "",
                                              TRUE ~ NA_character_))

    method.data3$Value <- round(method.data3$Value, 2)
    method.data3$text <- NULL
    # text
    for (i in 1:nrow(method.data3)) {
      if (method.data3$pvalue[i]=="") {
        method.data3$text[i] <- method.data3$Value[i]
      }else{
        method.data3$text[i] <- paste0(method.data3$Value[i], "\n", "(", method.data3$pvalue[i], ")" )
      }

    }


    # plot
    p.heatmap <- ggplot2::ggplot(method.data3, ggplot2::aes(x = Module, y = Variable, fill = Value)) +
      ggplot2::geom_tile(color = "black") +
      ggplot2::geom_text(ggplot2::aes(label = text), color = "black", show.legend = T)+
      ggplot2::scale_fill_gradientn(colours = grDevices::colorRampPalette(correlation.color)(100))+
      ggplot2::theme(panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.text.x.bottom = ggplot2::element_blank())+
      ggplot2::xlab("")+
      ggplot2::ylab("")+
      ggplot2::labs(fill = "Correlation")+
      Seurat::NoLegend()

    # label
    if (is.null(method.color)) {
      method.color <- ggsci::pal_npg()(length(unique(method.data3$Module)))
    }
    p.label <- method.data3 %>%
      dplyr::mutate(Method = "Method") %>%
      ggplot2::ggplot(ggplot2::aes(Module, y = Method, fill = Module)) +
      ggplot2::ylab("Method")+
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(values = method.color, name = "Method") +
      ggplot2::scale_y_discrete(position = "left") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(angle = 0, size =10, vjust = 0.5),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank())+
      Seurat::NoLegend()

    p <- p.heatmap %>% aplot::insert_top(p.label, height = 0.1)
    p <- ggplotify::as.ggplot(p) +
      ggplot2::ggtitle(unique(method.data3$geneset))+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    # make legend
    col_fun <- circlize::colorRamp2(c(stats::fivenum(method.data3$Value)[1],
                                      stats::fivenum(method.data3$Value)[3],
                                      stats::fivenum(method.data3$Value)[5]),
                                    correlation.color)

    lgd.method <- ComplexHeatmap::Legend(labels = levels(factor(method.data3$Module)),
                                         title = "Method",
                                         legend_gp = grid::gpar(fill = method.color))


    lgd.cor <- ComplexHeatmap::Legend(col_fun = col_fun, title = "Correlation")


    lgd.p <- ComplexHeatmap::Legend(pch = c("*", "**", "***", "****"),
                                    type = "points",
                                    labels = c("<= 0.05", "< 0.01", "< 0.001", "< 0.0001"),
                                    title = "Spearman's P Value")

    heatmap.legend <- ComplexHeatmap::packLegend(lgd.method, lgd.cor, lgd.p,
                                                 direction = "vertical",
                                                 column_gap = grid::unit(1, "cm"))
    heatmap.legend <- grid::grid.grabExpr(ComplexHeatmap::draw(heatmap.legend)) %>%
      ggplotify::as.ggplot()

    p <- cowplot::plot_grid(p, heatmap.legend, rel_widths = c(1,0.1))

    return(p)

  })

  result <- list()
  result[["hub_result"]] <- method.data
  result[["hub_plot"]] <- method.data2
  return(result )

}


