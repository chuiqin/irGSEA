


#' Integrate differential gene set calculated by all enrichment score matrixes
#'
#' Wlicox test is perform to all enrichment score matrixes and gene sets with
#' adjusted p value < 0.05 are used to integrated through RRA.
#'
#' @param object Seurat object including enrichment score matrixes.
#' @param group.by Default ident when it set to NULL. You can specify other
#' columns of metadata.
#' @param metadata A vector add to this Seurat object, and then it is set as
#' the ident of Seurat object to perform differential gene sets analysis.
#' @param col.name A name for metadata.
#' @param method A vector to select enrichment score matrixes. Default all
#' enrichment score matrixes
#' @param p.value p_val_adj or p_val. The gene sets with statistical significance
#' in different scoring methods were filtered based on the p_val_adj or p_val
#' and then RRA analysis was performed. Default p_val_adj.
#'
#' @return A list including the differential gene sets calculated by enrichment
#' score matrixes through wlicox test. Gene sets with adjusted p value < 0.05
#' are regarded as statistically significant. Meanwhile, significant differential
#' gene sets are integrated and the results are saved in a data frame named RRA.
#'  Among them, Gene sets with p value < 0.05 are statistically significant and
#'  common differential in all gene sets enrichment analysis methods.
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
#' result.dge2 <- irGSEA.integrate(object = pbmc3k.final, group.by = NULL,
#' metadata = pbmc3k.final$seurat_annotations, col.name = "ident",
#' method = c("AUCell","UCell","singscore","ssgsea"))
#'
#' }
#'

irGSEA.integrate <- function(object = NULL, group.by = NULL,
                             metadata = NULL, col.name = NULL,
                             method = c("AUCell","UCell","singscore","ssgsea"),
                             p.value = "p_val_adj"){
  # method
  if (all(method %in% Seurat::Assays(object)) == F) {
    stop("Please imput correct `method`.")
  }
  if (! all(method %in% Seurat::Assays(object)) == T) {
    method <- method[method %in% Seurat::Assays(object)]
    message(paste0("object only includes : ",stringr::str_c(method, collapse = ", ")))
    message(paste0("We only calculate : ",stringr::str_c(method, collapse = ", ")))
  }
  # group
  if (purrr::is_null(metadata)) {
    if (purrr::is_null(group.by)) {
      anno.ident <- SeuratObject::Idents(object)
    }else{
      object <- SeuratObject::SetIdent(object, value = group.by)
      anno.ident <- SeuratObject::Idents(object)

    }
  }else{
    if (purrr::is_null(colnames)) {
      stop("You need to input both `metadata` and `col.name`.")
    }
    object <- SeuratObject::AddMetaData(object, metadata = metadata, col.name = col.name)
    object <- SeuratObject::SetIdent(object, value = col.name)
    anno.ident <- SeuratObject::Idents(object)
  }

  anno.ident <- as.factor(as.character(anno.ident))
  SeuratObject::Idents(object) <- anno.ident

  # calculate
  deg.geneset <- list()
  target.gene.geneset <- list()
  avg_diff <- NULL
  unloadNamespace("VISION")
  for (i in seq_along(method)) {
    message(paste0("Calculate differential gene set", " : ", method[i]))
    result.wilcox <- tryCatch({
      marker.geneset <- lapply(levels(anno.ident), function(x){
        #  check Seurat version
        if (utils::packageVersion("Seurat") >= "5.0.2") {
          a <- Seurat::FindMarkers(object = object,
                                   assay = method[i],
                                   slot = "scale.data",
                                   fc.slot = "scale.data",
                                   ident.1 = x,
                                   ident.2 = NULL,
                                   test.use = "wilcox",
                                   min.pct = -Inf,
                                   logfc.threshold = 0,
                                   min.cells.group = 0,
                                   min.diff.pct = -Inf,
                                   verbose = F,
                                   min.cells.feature = 0)
        }else{
          a <- Seurat::FindMarkers(object = object,
                                   assay = method[i],
                                   slot = "scale.data",
                                   ident.1 = x,
                                   ident.2 = NULL,
                                   test.use = "wilcox",
                                   min.pct = -Inf,
                                   logfc.threshold = 0,
                                   min.cells.group = 0,
                                   min.diff.pct = -Inf,
                                   verbose = F,
                                   min.cells.feature = 0)
        }

        a <- a %>% tibble::rownames_to_column(var = "gene") %>%
          dplyr::mutate(cluster = x, direction = dplyr::if_else(avg_diff >0, "up", "down")) %>%
          dplyr::select(-c("pct.1", "pct.2"))
      })
      marker.geneset <- do.call(rbind, marker.geneset)
      deg.geneset[[method[i]]] <- cbind(marker.geneset, methods = method[i])

      # the target gene of gene set
      gene <- NULL
      geneset <- NULL
      if (class(object[[method[i]]])[1] == "Assay5") {
        deg.geneset[[method[i]]] <- deg.geneset[[method[i]]] %>%
          dplyr::mutate(target.gene = plyr::mapvalues(gene,
                                                      from = rownames(object[[method[i]]]),
                                                      to = object[[method[i]]]@meta.data$target.gene))
        target.gene <- NULL
        target.gene.geneset[[i]] <- object[[method[i]]]@meta.data %>%
          dplyr::mutate(geneset = rownames(object[[method[i]]])) %>%
          dplyr::relocate(geneset, .before = target.gene)
      }else{
        deg.geneset[[method[i]]] <- deg.geneset[[method[i]]] %>%
          dplyr::mutate(target.gene = plyr::mapvalues(gene,
                                                      from = rownames(object[[method[i]]]@meta.features),
                                                      to = object[[method[i]]]@meta.features$target.gene))
        target.gene.geneset[[i]] <- object[[method[i]]]@meta.features %>%
          tibble::rownames_to_column(var = "geneset")
      }

      message("Finish!")

    }, error = identity)
    if (methods::is(result.wilcox, "error")) {next}

  }

  deg.geneset.list <- deg.geneset %>% purrr::map( ~.x %>% dplyr::rename(Name = gene))
  deg.geneset <- do.call(rbind, deg.geneset)

  target.gene.geneset <- do.call(rbind, target.gene.geneset)
  target.gene.geneset <- target.gene.geneset %>%
    dplyr::distinct(geneset, .keep_all = T)

  #### RRA  ####
  p_val_adj <- NULL
  p_val <- NULL
  cluster <- NULL
  methods <- NULL

  if (p.value == "p_val_adj") {
    deg.cluster <- deg.geneset %>%
      dplyr::filter(p_val_adj <= 0.05) %>%
      dplyr::select(c("avg_diff", "cluster", "gene","methods"))
  }else{
    deg.cluster <- deg.geneset %>%
      dplyr::filter(p_val <= 0.05) %>%
      dplyr::select(c("avg_diff", "cluster", "gene","methods"))
  }

  # if deg.cluster is empty
  if (nrow(deg.cluster) == 0) {
    stop(cat("None of the scoring methods: (", method,
             ") produced statistically significant gene sets. ",
             "Change the parameter `p.value` from `p_val_adj` to `p_val`, ",
             "or choose other scoring methods.",
             collapse  = ""))
  }else{
    if (any(!method %in% unique(deg.cluster$methods))) {
      cat("None of the scoring methods: (",
          method[!method %in% unique(deg.cluster$methods)],
          ") produced statistically significant gene sets.",collapse  = "")
      cat("\n")
      cat("These methods will not be included in the subsequent RRA analysis.")
    }
  }


  deg.cluster$cluster <- as.factor(as.character(deg.cluster$cluster))
  deg.cluster <- deg.cluster %>%
    dplyr::group_split(cluster) %>%
    purrr::set_names(levels(deg.cluster$cluster))

  deg.cluster.postive <- lapply(deg.cluster, function(x){
    a <- x %>%
      dplyr::filter(avg_diff > 0) %>%
      dplyr::arrange(methods, dplyr::desc(avg_diff))
    a$methods <- as.factor(as.character(a$methods))
    b <- a %>%
      dplyr::group_split(methods) %>%
      purrr::set_names(levels(a$methods))
  })

  deg.cluster.negative <- lapply(deg.cluster, function(x){
    a <- x %>%
      dplyr::filter(avg_diff < 0) %>%
      dplyr::arrange(methods, avg_diff)
    a$methods <- as.factor(as.character(a$methods))
    b <- a %>%
      dplyr::group_split(methods) %>%
      purrr::set_names(levels(a$methods))

  })

  if (! identical(names(deg.cluster.postive), names(deg.cluster.postive %>% purrr::compact()))) {
    a <- setdiff(names(deg.cluster.postive), names(deg.cluster.postive %>% purrr::compact()))
    a <- stringr::str_c(a, collapse = ", ")
    message(paste0("No sigficant genesets in cluster : ", a, " after wilicox test"))
    deg.cluster.postive <- deg.cluster.postive %>% purrr::compact()
  }
  if (! identical(names(deg.cluster.negative), names(deg.cluster.negative %>% purrr::compact()))) {
    a <- setdiff(names(deg.cluster.negative), names(deg.cluster.negative %>% purrr::compact()))
    a <- stringr::str_c(a, collapse = ", ")
    message(paste0("No sigficant genesets in cluster : ", a, " after wilicox test"))
    deg.cluster.negative <- deg.cluster.negative %>% purrr::compact()
  }

  # perform RRA analysis
  sig.genesets.postive <- lapply(deg.cluster.postive, function(x){
    x <- lapply(x,function(y){y$gene})
    x <- RobustRankAggreg::aggregateRanks(glist = x, N = length(unique(unlist(x))))
  })
  sig.genesets.negative <- lapply(deg.cluster.negative,function(x){
    x <- lapply(x,function(y){y$gene})
    x <- RobustRankAggreg::aggregateRanks(glist = x, N = length(unique(unlist(x))))
  })
  Name <- NULL
  sig.genesets.postive <- sig.genesets.postive %>%
    reshape2::melt(id.vars = "Name") %>%
    dplyr::mutate(direction = "up")
  sig.genesets.negative <- sig.genesets.negative %>%
    reshape2::melt(id.vars = "Name") %>%
    dplyr::mutate(direction = "down")
  sig.genesets <- rbind(sig.genesets.postive, sig.genesets.negative)
  colnames(sig.genesets)[4] <- "cluster"
  value <- NULL
  sig.genesets <- sig.genesets %>%
    dplyr::select(c("Name", "value","cluster","direction")) %>%
    dplyr::rename(pvalue = value) %>%
    dplyr::mutate(method = "RRA")

  deg.geneset.list[["RRA"]] <- sig.genesets
  deg.geneset.list[["RRA"]] <- deg.geneset.list[["RRA"]] %>%
    dplyr::mutate(target.gene = plyr::mapvalues(Name,
                                                from = target.gene.geneset$geneset,
                                                to = target.gene.geneset$target.gene))

  return(deg.geneset.list)

}



