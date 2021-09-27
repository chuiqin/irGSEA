#' Calculate enrichment scores from scRNA-seq data
#'
#' Input a Seurat object or scRNA-seq matrix, calculate the enrichment scores of
#' AUCell, UCell, singscore and ssgsea. Then, return a Seurat object including these
#' score matrix.
#'
#' @param object Seurat object or scRNA-seq matrix
#' @param assay Name of assay to use, defaults to the active assay. The parameter
#' works when a seurat object is input.
#' @param slot Default data. The parameter works if a seurat object is input.
#' @param min.cells The minimum detected cells per gene, default 3. The parameter
#' worsk if a scRNA-seq matrix is input.
#' @param min.feature The minimum genes per cell, default 0. The parameter works
#' if a scRNA-seq matrix is input.
#' @param seeds Default 123
#' @param ncores Default 4
#' @param custom Default False. Set it to true when input own genesets.
#' @param geneset Default NULL. Input own genesets as a list. Each element in
#' the list is a gene set. You can also specify positive and negative genes by
#' adding a + or - sign in special gene set.
#' @param msigdb Default True. You can acquire the collection gene sets from
#' msigdb database. It will be ignored when custom is set to true.
#' @param species Default Homo sapiens. Use msigdbr::msigdbr_show_species() to view
#'  all available species. The parameter works if msigdb is True.
#' @param category Default H. You can acquire the Hallmarker gene sets. Use
#' msigdbr::msigdbr_collections to view all available collections gene sets. The
#' parameter works if msigdb is True.
#' @param subcategory Default NULL. The parameter works if msigdb is True.
#' @param geneid Default symbol. Other options are "entrez" and "ensembl". The
#' parameter works if msigdb is True.
#' @param method A vector. Default all method.
#' @param aucell.MaxRank The threshold to calculate the AUC. Default only the top 5%
#'  of the expressed genes are used to checks whether the gene set are within
#'  the top 5% when it set to NULL. You can inputc special number, such as 1500,
#'   to change the threshold. The parameter works if "AUCell" is selected in "method".
#' @param ucell.MaxRank Maximum number of genes to rank per cell. Above this
#' rank, a given gene is considered as not expressed. Default 1500 when it set
#' to NULL. You can input special number, such as 2000, to increase the rank
#' range. The parameter works if "UCell" is selected in "method".
#' @param kcdf Default Gaussian if input expression values are continuous in
#' logarithmic scale. Other option is "Poisson" if input expression values are
#' integer counts. The parameter works if "ssgsea" is selected in "method".
#'
#' @return a Seurat object including score matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' # load PBMC dataset by R package SeuratData
#  # devtools::install_github('satijalab/seurat-data')
#' library(Seurat)
#' library(SeuratData)
#' # download 3k PBMCs from 10X Genomics
#' InstallData("pbmc3k")
#' data("pbmc3k.final")
#' pbmc3k.final <- SeuratObject::UpdateSeuratObject(pbmc3k.final)
#'
#' # Seurat object
#' pbmc3k.final <- irGSEA.score(object = pbmc3k.final, assay = "RNA", slot = "data",
#' msigdb = T, species = "Homo sapiens", category = "H", geneid = "symbol",
#' method = c("AUCell", "UCell", "singscore", "ssgsea"), kcdf = 'Gaussian')
#'
#' # data.fram, matrix, or dgmatrix
#' pbmc3k.final2 <- irGSEA.score(object = pbmc3k.final@@assays$RNA@@counts,
#' assay = "RNA", slot = "data", min.cells = 3, min.feature = 0,
#' method = c("AUCell", "UCell", "singscore", "ssgsea"), kcdf = 'Poisson')
#'
#' # custom geneset
#' markers <- list()
#' markers$Tcell_gd <- c("TRDC+", "TRGC1+", "TRGC2+", "TRDV1+","TRAC-","TRBC1-","TRBC2-")
#' markers$Tcell_NK <- c("FGFBP2+", "SPON2+", "KLRF1+", "FCGR3A+", "CD3E-","CD3G-")
#' markers$Tcell_CD4 <- c("CD4","CD40LG")
#' markers$Tcell_CD8 <- c("CD8A","CD8B")
#' markers$Tcell_Treg <- c("FOXP3","IL2RA")
#' pbmc3k.final3 <- irGSEA.score(object = pbmc3k.final, assay = "RNA", slot = "data",
#' custom = T, geneset = markers, method = c("AUCell", "UCell", "singscore", "ssgsea"),
#' kcdf = 'Gaussian')
#' }
#'
irGSEA.score <- function(object = NULL, assay = NULL, slot = "data",
                         min.cells = 3, min.feature = 0,
                         seeds = 123, ncores = 4,
                         custom = F, geneset = NULL,
                         msigdb = T, species = "Homo sapiens",
                         category = "H", subcategory = NULL, geneid = "symbol",
                         method = c("AUCell", "UCell", "singscore", "ssgsea"),
                         aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                         kcdf = 'Gaussian'){
  set.seed(seeds)
  #### prepare data
  if (all(methods::is(object)=="Seurat")){
    if (purrr::is_null(assay)){assay <- Seurat::DefaultAssay(object)}
    object <- SeuratObject::UpdateSeuratObject(object)
    my.matrix <- SeuratObject::GetAssayData(object, assay = assay, slot = slot)
  }else{
    object <- SeuratObject::CreateSeuratObject(counts = object,
                                               project = "irGSEA",
                                               assay = "RNA",
                                               min.cells = min.cells,
                                               min.feature = min.feature)
    my.matrix <- SeuratObject::GetAssayData(object, assay = "RNA", slot = "counts")}

  #### prepare geneset
  if(msigdb == T){
    h.human <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
    colnames(h.human) <- stringr::str_replace(colnames(h.human), "gene_symbol","symbol")
    colnames(h.human) <- stringr::str_replace(colnames(h.human), "entrez_gene", "entrez")
    colnames(h.human) <- stringr::str_replace(colnames(h.human), "ensembl_gene", "ensembl")
    h.human$gs_name <- as.factor(as.character(h.human$gs_name))
    gs_name <- NULL
    h.sets <- h.human %>%
      dplyr::select(c(gs_name, geneid)) %>%
      dplyr::group_split(gs_name, .keep = F) %>%
      purrr::set_names(levels(h.human$gs_name))

    h.sets <- lapply(h.sets, function(x){unique(unlist(x))})

    # concvert list to GeneSetCollection
    h.gsets <- GSEABase::GeneSetCollection(mapply(function(geneIds, keggId){
      GSEABase::GeneSet(geneIds, geneIdType = GSEABase::EntrezIdentifier(),
                        collectionType = GSEABase::KEGGCollection(keggId),
                        setName = keggId)
    }, h.sets, names(h.sets)))

    # filiter the gene set based on object
    h.gsets <- AUCell::subsetGeneSets(h.gsets, rownames(object))
    h.gsets.list <- GSEABase::geneIds(h.gsets)

    # message: genesets with zero genes after subset
    h.gsets.list.setdiff <- setdiff(names(h.gsets.list),
                                    names(h.gsets.list %>% purrr::compact()))
    if(! purrr::is_empty(h.gsets.list.setdiff)){
      message(paste0("No genes remaining in following genesets: ",
                     stringr::str_c(h.gsets.list.setdiff, collapse = ", ")))
      h.gsets.list <- h.gsets.list %>% purrr::compact()}

    # create new GeneSetCollection
    h.gsets <- GSEABase::GeneSetCollection(mapply(function(geneIds, keggId){
      GSEABase::GeneSet(geneIds, geneIdType = GSEABase::EntrezIdentifier(),
                        collectionType = GSEABase::KEGGCollection(keggId),
                        setName = keggId)
    }, h.gsets.list, names(h.gsets.list)))


  }
  if (custom == T) {
    # list names
    if(! purrr::is_list(geneset)){stop("Custom geneset should be a list.")}
    if(purrr::is_null(names(geneset))){names(geneset) <- paste0("geneset",seq_along(geneset))}

    # filiter the gene sets based on object
    for(i in seq_along(geneset)){
      if(all(stringr::str_detect(geneset[[i]], pattern = "\\+$|-$")==T)){
        geneset[[i]] <- geneset[[i]][stringr::str_remove(geneset[[i]],pattern = "\\+$|-$") %in% rownames(object)]
        if(length(geneset[[i]])==0){
          message(paste0("No genes remaining in following genesets: ",names(geneset[i])))}
        if(all(stringr::str_detect(geneset[[i]], pattern = "\\+$")==F)){
          message(paste0("No positive genes remaining in following genesets: ",names(geneset[i])))}
        if(all(stringr::str_detect(geneset[[i]], pattern = "-$")==F)){
          message(paste0("No negative genes remaining in following genesets: ",names(geneset[i])))}
      }else{
        if(any(stringr::str_detect(geneset[[i]], pattern = "\\+$|-$"))){
          message(paste0("All genes need direction in following genesets: ",names(geneset[i])))}
        geneset[[i]] <- geneset[[i]][geneset[[i]] %in% rownames(object)]
        if(length(geneset[[i]])==0){
          message(paste0("No genes remaining in following genesets: ",names(geneset[i])))}
      }
    }

    h.gsets.list <- geneset %>% purrr::compact()


  }

  #### calculate scores
  if ("AUCell" %in% method) {
    message("Calculate AUCell scores")
    aucell.rank <- AUCell::AUCell_buildRankings(my.matrix,
                                                nCores = ncores,
                                                plotStats = F,
                                                verbose = F)
    h.gsets.list.aucell <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})
    if (purrr::is_null(aucell.MaxRank)){aucell.MaxRank = ceiling(0.05 * nrow(aucell.rank))}
    aucell.scores <- AUCell::AUCell_calcAUC(h.gsets.list.aucell,
                                            aucell.rank,
                                            nCores = ncores,
                                            aucMaxRank = aucell.MaxRank,
                                            verbose = F)
    aucell.scores <- SummarizedExperiment::assay(aucell.scores)
    object[["AUCell"]] <- SeuratObject::CreateAssayObject(counts = aucell.scores)
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = aucell.scores, assay = "AUCell")
    message("Finish calculate AUCell scores")
  }
  if ("UCell" %in% method) {
    message("Calculate UCell scores")
    if (purrr::is_null(ucell.MaxRank)){ucell.MaxRank = 1500}
    ucell.scores <- UCell::ScoreSignatures_UCell(matrix = my.matrix,
                                                 features = h.gsets.list,
                                                 maxRank = ucell.MaxRank,
                                                 w_neg = 1,
                                                 ncores = ncores,
                                                 force.gc = T,
                                                 seed = seeds)
    colnames(ucell.scores) <- stringr::str_remove(colnames(ucell.scores), pattern = "_UCell")
    object[["UCell"]] <- SeuratObject::CreateAssayObject(counts = t(ucell.scores))
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = t(ucell.scores), assay = "UCell")
    message("Finish calculate UCell scores")
  }
  if ("singscore" %in% method) {
    message("Calculate singscore scores")
    # calculate the rank matrix
    attemptsLeft <- 20
    while(attemptsLeft > 0){
      singscore.rank <- tryCatch(singscore::rankGenes(as.data.frame(my.matrix)), warning = identity)
      if(methods::is(singscore.rank, "warning")){
        attemptsLeft <- attemptsLeft - 1
        Sys.sleep(5)
      }else{attemptsLeft <- 0}}

    # calculate separately
    singscore.scores <- list()
    for (i in seq_along(h.gsets.list)){
      if (any(stringr::str_detect(h.gsets.list[[i]], pattern = "\\+$|-$"))) {
        h.gsets.list.positive <- stringr::str_match(h.gsets.list[[i]],pattern = "(.+)\\+")[,2] %>% purrr::discard(is.na)
        h.gsets.list.negative <- stringr::str_match(h.gsets.list[[i]],pattern = "(.+)-")[,2] %>% purrr::discard(is.na)
        if (length(h.gsets.list.positive)==0) {
          singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,
                                                          upSet = h.gsets.list.negative)
        }
        if (length(h.gsets.list.negative)==0) {
          singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,
                                                          upSet = h.gsets.list.positive)
        }
        if ((length(h.gsets.list.positive)!=0)&(length(h.gsets.list.negative)!=0)) {
          singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,
                                                          upSet = h.gsets.list.positive,
                                                          downSet = h.gsets.list.negative)
        }

      }else{
        singscore.scores[[i]] <- singscore::simpleScore(singscore.rank, upSet = h.gsets.list[[i]])
      }
      TotalScore <- NULL
      singscore.scores[[i]] <- singscore.scores[[i]] %>%
        dplyr::select(TotalScore) %>%
        magrittr::set_colnames(names(h.gsets.list)[i])}
    names(singscore.scores) <- names(h.gsets.list)
    singscore.scores <- do.call(cbind, singscore.scores)
    object[["singscore"]] <- SeuratObject::CreateAssayObject(counts = t(singscore.scores))
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = t(singscore.scores), assay = "singscore")
    message("Finish calculate singscore scores")

  }
  if ("ssgsea" %in% method) {
    message("Calculate ssgsea scores")
    h.gsets.list.ssgsea <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})
    ssgsea.scores <- GSVA::gsva(as.matrix(my.matrix),
                                h.gsets.list.ssgsea,
                                method = "ssgsea",
                                kcdf = kcdf,
                                ssgsea.norm = F,
                                parallel.sz = ncores,
                                verbose = F)
    object[["ssgsea"]] <- SeuratObject::CreateAssayObject(counts = ssgsea.scores)
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = ssgsea.scores, assay = "ssgsea")
    message("Finish calculate ssgsea scores")
  }
  return(object)

}




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
                             method = c("AUCell","UCell","singscore","ssgsea")){
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
  avg_diff <- NULL
  for (i in seq_along(method)) {
    message(paste0("Calculate differential gene set", " : ", method[i]))

    marker.geneset <- lapply(levels(anno.ident), function(x){
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
      a <- a %>% tibble::rownames_to_column(var = "gene") %>%
        dplyr::mutate(cluster = x, direction = dplyr::if_else(avg_diff >0, "up", "down")) %>%
        dplyr::select(-c("pct.1", "pct.2"))
    })
    marker.geneset <- do.call(rbind, marker.geneset)
    deg.geneset[[i]] <- cbind(marker.geneset, methods = method[i])

  }
  names(deg.geneset) <- method
  deg.geneset.list <- deg.geneset %>% purrr::map( ~.x %>% dplyr::rename(Name = gene))
  deg.geneset <- do.call(rbind, deg.geneset)
  p_val_adj <- NULL
  cluster <- NULL
  methods <- NULL
  deg.cluster <- deg.geneset %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::select(c("avg_diff", "cluster", "gene","methods"))

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

  sig.genesets.postive <- sig.genesets.postive %>%
    reshape2::melt() %>%
    dplyr::mutate(direction = "up")
  sig.genesets.negative <- sig.genesets.negative %>%
    reshape2::melt() %>%
    dplyr::mutate(direction = "down")
  sig.genesets <- rbind(sig.genesets.postive, sig.genesets.negative)
  colnames(sig.genesets)[4] <- "cluster"
  value <- NULL
  sig.genesets <- sig.genesets %>%
    dplyr::select(c("Name", "value","cluster","direction")) %>%
    dplyr::rename(pvalue = value) %>%
    dplyr::mutate(method = "RRA")
  deg.geneset.list[["RRA"]] <- sig.genesets
  return(deg.geneset.list)

}



