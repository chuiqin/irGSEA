#' Calculate enrichment scores from scRNA-seq data
#'
#' Input a Seurat object or scRNA-seq matrix, calculate the enrichment scores of
#' AUCell, UCell, singscore, ssgsea, JASMINE and viper. Then, return a Seurat object
#' including these score matrix.
#'
#' @param object Seurat object or scRNA-seq matrix
#' @param assay Name of assay to use, defaults to the active assay. The parameter
#' works when a seurat object is input.
#' @param slot Default data. The parameter works if a seurat object is input.
#' VISION or AddModuleScore should be data, and VAM or pagoda2 should be counts.
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
#' @param geneset.weight Default NULL. Input is a list and the length of list should
#' be equal to geneset. Each element in the list represent a gene set. The element
#' is a vector. And the name of vector are gene names. And you can also specify
#' the weight of each gene by a specific number. The parameter works if parameter
#' "geneset" is not null and "wsum", "wmean", "mdt", "viper" are selected in "method".
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
#' @param minGSSize Minimum number of genes in one gene set.
#' @param maxGSSize Maximal number of genes in one gene set.
#' @param method A vector. Default  c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper").
#'               `AUCell (https://doi.org/10.1038/nmeth.4463)`:
#'               AUCell uses the area-under-the-curve (AUC)  to calculate whether
#'               a gene set is enriched within the molecular readouts of each cell.
#'               AUC can be calculated using by default thetop 5% molecular features
#'               in the ranking.
#'
#'               `UCell (https://doi.org/10.1016/j.csbj.2021.06.043)`:
#'               UCell calculates gene signature scores based on the Mann-Whitney
#'               U statistic. And the U statistic is closely related to the
#'               area-under-the-curve (AUC) statistic for ROC curves
#'
#'               `singscore (https://doi.org/10.1093/nar/gkaa802)`:
#'               Singscore uses rank-based statistics to analyze each sample’s
#'               gene expression profile and scores the expression activities of
#'               gene sets at a single-sample level.
#'
#'               `ssgsea (https://doi.org/10.1038/nature08460)`:
#'               ssGSEA ranks gene expression within each cell separately, then
#'               the enrichment score of the gene set of each cell is calculated
#'               by K-S like random walk statistic.
#'
#'               `JASMINE (https://doi.org/10.7554/eLife.71994)`:
#'               JASMINE calculates the approximate mean using gene ranks among
#'               expressed genes and the enrichment of the signature in expressed
#'               genes. The two are then scaled to 0–1 and averaged to result in
#'               the final JASMINE score. JASMINE considers the enrichment of
#'               signature genes in expressed genes to counter dropout effects,
#'               and meanwhile, evaluates the average expression level of
#'               the expressed signature genes.
#'
#'               `VAM (https://doi.org/10.1093/nar/gkaa582)`:
#'               VAM generates scores from scRNA-seq data using a variation of
#'               the classic Mahalanobis multivariate distance measure.
#'
#'               `scSE (https://doi.org/10.1093/nar/gkz601)`:
#'               scSE, single cell signature explorer, measures a signature using
#'               normalized total expression of the signature genes.
#'
#'               `VISION (https://doi.org/10.1038/s41467-019-12235-0)`:
#'               Scores in Vision are calculated by averaging expressed genes for
#'               each gene set. To account for the influence of sample-level
#'               metrics (the number of UMIs/reads per cells), scores are then
#'               corrected by their means and standard deviations.
#'
#'               `wsum (https://doi.org/10.1093/bioadv/vbac016)`:
#'               First, multiply each gene by its associated weight which then
#'               are summed to an enrichment score wsum. Furthermore, permutations
#'               of random gene can be performed to obtain a null distribution that
#'               can be used to compute a z-score norm_wsum, or a corrected
#'               estimate corr_wsum by multiplying wsum by the minus log10 of
#'               the obtained empirical p-value. We use corr_wsum as the enrichment
#'               score for the gene set.
#'
#'               `wmean (https://doi.org/10.1093/bioadv/vbac016)`:
#'               Weighted Mean (WMEAN) is similar to WSUM but it divides the obtained
#'               enrichment score by the sum of the absolute value of weights
#'
#'               `mdt (https://doi.org/10.1093/bioadv/vbac016)`:
#'               MDT fits a multivariate regression random forest for each sample.
#'
#'               `viper (https://doi.org/10.1038/ng.3593)`:
#'               VIPER estimates biological activities by performing a three-tailed
#'               enrichment score calculation. First by a one-tail approach, based
#'               on the absolute value of the gene expression signature (i.e.,
#'               genes are rank-sorted from the less invariant between groups to
#'               the most differentially expressed, regardless of the direction
#'               of change); and then by a two-tail approach, where the positions
#'               of the genes whose expression is repressed by the regulator
#'               are inverted in the gene expression signature before computing
#'               the enrichment score. The one-tail and two-tail enrichment score
#'               estimates are integrated while weighting their contribution based
#'               on the estimated mode of regulation through three-tail approach.
#'               The contribution of each target gene from a given regulon to the
#'               enrichment score is also weighted based on the regulator-target
#'               gene interaction confidence.
#'
#'               `GSVApy (https://doi.org/10.1038/ng.3593)`:
#'               The python version of GSVA is wrapped by the decoupleR package.
#'
#'               `gficf (https://doi.org/10.1093/nargab/lqad024)`:
#'               gficf takes advantage of the informative biological signals spreading
#'               across the latent factors of gene expression values obtained from
#'               non-negative matrix factorization. It uses NMF and FGSEA to estimate
#'               enrichment scores.
#'
#'               `GSVA (https://doi.org/10.1186/1471-2105-14-7)`:
#'               GSVA, Gene Set Variation Analysis, starts by transforming the input
#'               molecular readouts matrix to a readout-level statistic using Gaussian
#'               kernel estimation of the cumulative density function. Then, readout-level
#'               statistics are ranked per sample and normalized to up-weight the two tails
#'               of the rank distribution. Afterwards, an enrichment score (ES)
#'               is calculated as in GSEA, using the running sum statistic. Finally,
#'               the ES can be normalized by subtracting the largest negative ES
#'               from the largest positive ES.
#'
#'               `zscore (https://doi.org/10.1371/journal.pcbi.1000217)`:
#'               zscore aggregates the expression of several genes of the gene set. The
#'               gene expression is scaled by mean and standard deviation over cells.
#'               Then, the enrichment scores for each cell are calculated by averaging
#'               scaled gene expression of all genes within each gene set.
#'
#'               `plage (https://doi.org/10.1186/1471-2105-6-225)`:
#'               plag captures enrichment scores from singular value decompositions (SVD)
#'               strategy. PLAGE first standardizes gene expression matrix across cells.
#'               For a submatrix which genes in a particular gene set, the first
#'               coefficient of right-singular vector in SVD of this matrix is
#'               extracted as enrichment scores.
#'
#'               `AddModuleScore (https://doi.org/10.1126/science.aad0501)`:
#'               Calculate the average expression levels of each program (cluster)
#'               on single cell level, subtracted by the aggregated expression of
#'               control feature sets. All analyzed features are binned based
#'               on averaged expression, and the control features are randomly
#'               selected from each bin.
#'
#'               `pagoda2 (https://doi.org/10.1038/nbt.4038)`:
#'               pagoda2 fits an error model for each cell to depict its properties,
#'               and residual variance of each gene in the cell is re-normalized
#'               subsequently. Then, the enrichment scores of each gene set is
#'               quantified by its first weighted principal component.
#'
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
#' @param JASMINE.method Default oddsratio. You can choose oddsratio or likelihood.
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
                         custom = F, geneset = NULL, geneset.weight = NULL,
                         msigdb = T, species = "Homo sapiens",
                         category = "H", subcategory = NULL,
                         geneid = "symbol", minGSSize = 1, maxGSSize = 500,
                         method = c("AUCell", "UCell", "singscore", "ssgsea"),
                         aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                         kcdf = 'Gaussian', JASMINE.method = "oddsratio"){

  #### Set random seeds
  set.seed(seeds)

  #### 00.prepare data ####
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

  #### 01.prepare geneset ####

  if(msigdb == T & custom == F){

    if ( utils::available.packages()["msigdbr", "Version"] > utils::packageVersion("msigdbr")) {
      message("There is a newer version of the msigdbr package. Watch out for updates!")
      message("You can update the msigdbr package via `install.packages(`msigdbr`)`.")
      message(paste0("Your current version of the msigdbr package: ", utils::packageVersion("msigdbr")))
    }

    h.human <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
    colnames(h.human) <- stringr::str_replace(colnames(h.human), "gene_symbol","symbol")
    colnames(h.human) <- stringr::str_replace(colnames(h.human), "entrez_gene", "entrez")
    colnames(h.human) <- stringr::str_replace(colnames(h.human), "ensembl_gene", "ensembl")
    h.human$gs_name <- as.factor(as.character(h.human$gs_name))

    gs_name <- NULL
    h.sets <- h.human %>%
      dplyr::select(c(gs_name, tidyselect::all_of(geneid))) %>%
      dplyr::group_split(gs_name, .keep = F) %>%
      purrr::set_names(levels(h.human$gs_name))

    # remove duplication
    h.sets <- lapply(h.sets, function(x){unique(unlist(x))})

    # concvert list to GeneSetCollection
    h.gsets <- GSEABase::GeneSetCollection(mapply(function(geneIds, keggId){
      GSEABase::GeneSet(geneIds, geneIdType = GSEABase::EntrezIdentifier(),
                        collectionType = GSEABase::KEGGCollection(keggId),
                        setName = keggId)
    }, h.sets, names(h.sets)))

    # filiter the gene set based on object
    Seurat::DefaultAssay(object) <- assay
    h.gsets <- AUCell::subsetGeneSets(h.gsets, rownames(object))
    h.gsets.list <- GSEABase::geneIds(h.gsets)

    # message: genesets with zero genes after subset
    h.gsets.list.setdiff <- setdiff(names(h.gsets.list),
                                    names(h.gsets.list %>% purrr::compact()))
    if(! purrr::is_empty(h.gsets.list.setdiff)){
      message(paste0("No genes remaining in following genesets: ",
                     stringr::str_c(h.gsets.list.setdiff, collapse = ", ")))
      h.gsets.list <- h.gsets.list %>% purrr::compact()}

    # Filter gene sets by minGSSize and maxGSSize
    h.gsets.list2 <- sapply(h.gsets.list, function(x){
      length(x) >= minGSSize & length(x) <= maxGSSize
    })
    if (! all(h.gsets.list2)) {
      message(paste0("The number of genes is not between ", minGSSize, " and ", maxGSSize, " in following genesets: "))
      message(stringr::str_c(names(h.gsets.list)[!h.gsets.list2], collapse = ", "))
      message("These gene sets will be filtered and will not proceed with the analysis.")
    }
    h.gsets.list <- h.gsets.list[h.gsets.list2]

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
    # Filter gene sets by minGSSize and maxGSSize
    h.gsets.list2 <- sapply(h.gsets.list, function(x){
      length(x) >= minGSSize & length(x) <= maxGSSize
    })
    if (! all(h.gsets.list2)) {
      message(paste0("The number of genes is not between ", minGSSize, " and ", maxGSSize, " in following genesets: "))
      message(stringr::str_c(names(h.gsets.list)[!h.gsets.list2], collapse = ", "))
      message("These gene sets will be filtered and will not proceed with the analysis.")
    }
    h.gsets.list <- h.gsets.list[h.gsets.list2]

  }

  #### 02.check geneset ####
  if (purrr::is_empty(h.gsets.list)) {
    stop("All gene sets are empty after filtering, there may be the following 3 reasons:
1. The species represented by the genes in the gene set do not match the species represented by the genes in the single-cell matrix;
2. The gene of the gene set does not exist in the row name of the single-cell matrix;
3. The values of minGSSize and maxGSSize are inappropriate.")
  }

  if (!is.null(geneset.weight)) {
    geneset.weight <- geneset.weight[names(h.gsets.list)]
    for (i in seq_along(geneset.weight)) {
      geneset.weight[[i]] <- geneset.weight[[i]][h.gsets.list[[i]]]
    }
  }

  #### 03.check matrix and package ####

  # split the matrix if the matrix is too large
  if (ncol(my.matrix) >= 50000) {
    cut.times <- ceiling(ncol(my.matrix)/5000)
    my.matrix.list <- split(seq_along(colnames(my.matrix)),
                            cut(seq_along(colnames(my.matrix)), cut.times))
  }else{
    my.matrix.list <- list()
    my.matrix.list[[1]] <- seq_along(colnames(my.matrix))
  }

  ## Check if the R package is installed and the assay is suitable

  ## VAM
  if ("VAM" %in% method) {
    if (!(assay %in% c("RNA", "SCT")) & (slot %in% c("data"))) {
      stop("VAM only support assay (RNA or SCT) and slot (data).")
    }
    # install package from CRAN
    if (!requireNamespace("VAM", quietly = TRUE)) {
      message("install VAM package from CRAN")
      utils::install.packages("VAM", ask = F, update = F)
    }
  }

  ## VISION
  if ("VISION" %in% method) {
    if (! slot %in% c("counts")) {
      stop("VISION only support slot (counts).")
    }

    # install package from Github
    if (!requireNamespace("VISION", quietly = TRUE)) {
      message("install VISION package from Github")
      devtools::install_github("YosefLab/VISION", force =T)
    }
  }

  # GFICF
  if ("gficf" %in% method) {
    # install.package from Bioconductor
    if (!requireNamespace(c("sva","edgeR", "fgsea"), quietly = TRUE)) {
      BiocManager::install(c("sva","edgeR", "fgsea"), ask = F, update = F)
    }

    if (!requireNamespace("gficf", quietly = TRUE)) {
      message("install gficf package from Github")
      devtools::install_github("gambalab/gficf", force =T)
    }
  }

  # pagoda2
  if ("pagoda2" %in% method) {

    if (! slot %in% c("counts")) {
      stop("pagoda2 only support slot (counts).")
    }

    # install.package from CRAN
    if (!requireNamespace("pagoda2", quietly = TRUE)) {
      utils::install.packages("pagoda2", ask = F, update = F)
    }
    # install.package from Bioconductor
    if (!requireNamespace("scde", quietly = TRUE)) {
      BiocManager::install("scde", ask = F, update = F, force = T)
    }
  }

  # viper
  if (any(c("viper") %in% method)) {

    # install.package from Bioconductor
    if (!requireNamespace("decoupleR", quietly = TRUE)) {
      BiocManager::install("decoupleR", ask = F, update = F, force = T)
    }
    # install.package from Bioconductor
    if (!requireNamespace("viper", quietly = TRUE)) {
      BiocManager::install("viper", ask = F, update = F, force = T)
    }

  }

  # GSVApy
  if (any(c("GSVApy") %in% method)) {

    # install.package from CRAN
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      utils::install.packages("reticulate", ask = F, update = F)
    }

    # install.package from Github
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      devtools::install_github("mojaveazure/seurat-disk", force =T)
    }

    # reticulate::py_config()
    if (! "irGSEA" %in% reticulate::conda_list()$name) {
      reticulate::conda_create("irGSEA")
    }
  }


  #### 04.calculate AUCell scores ####
  # ties.method: "random"
  tryCatch({if ("AUCell" %in% method) {
    message("Calculate AUCell scores")
    h.gsets.list.aucell <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})
    aucell.scores.list <- list()
    for (k in seq_along(my.matrix.list)) {
      # calculate the rank matrix
      aucell.rank <- AUCell::AUCell_buildRankings(my.matrix[, my.matrix.list[[k]]],
                                                  plotStats = F,
                                                  verbose = F,
                                                  splitByBlocks = TRUE)
      if (purrr::is_null(aucell.MaxRank)){aucell.MaxRank = ceiling(0.05 * nrow(aucell.rank))}
      aucell.scores <- AUCell::AUCell_calcAUC(h.gsets.list.aucell,
                                              aucell.rank,
                                              nCores = ncores,
                                              aucMaxRank = aucell.MaxRank,
                                              verbose = F)
      aucell.scores <- SummarizedExperiment::assay(aucell.scores)
      aucell.scores.list[[k]] <- aucell.scores
    }
    aucell.scores.list <- do.call(cbind, aucell.scores.list)
    object[["AUCell"]] <- SeuratObject::CreateAssayObject(counts = aucell.scores.list)
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = aucell.scores.list, assay = "AUCell")
    message("Finish calculate AUCell scores")
    rm(aucell.rank)
    rm(aucell.scores.list)
    gc()

  }}, error = identity)



  #### 05.calculate UCell scores ####
  # ties.method: "average"
  tryCatch({if ("UCell" %in% method) {
    message("Calculate UCell scores")
    if (purrr::is_null(ucell.MaxRank)){ucell.MaxRank = 1500}
    ucell.scores <- UCell::ScoreSignatures_UCell(matrix = my.matrix,
                                                 features = h.gsets.list,
                                                 maxRank = ucell.MaxRank,
                                                 w_neg = 1,
                                                 ncores = ncores,
                                                 force.gc = T)
    colnames(ucell.scores) <- stringr::str_remove(colnames(ucell.scores), pattern = "_UCell")
    object[["UCell"]] <- SeuratObject::CreateAssayObject(counts = t(ucell.scores))
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = t(ucell.scores), assay = "UCell")
    message("Finish calculate UCell scores")
    rm(ucell.scores)
    gc()
  }}, error = identity)


  #### 06.calculate singscore scores ####
  # ties.method: "min"
  tryCatch({if ("singscore" %in% method) {
    message("Calculate singscore scores")
    singscore.scores.list <- list()
    for (k in seq_along(my.matrix.list)) {
      # calculate the rank matrix
      singscore.rank <- singscore::rankGenes(as.matrix(my.matrix[, my.matrix.list[[k]]]))
      # calculate separately
      singscore.scores <- list()

      for (i in seq_along(h.gsets.list)){
        if (any(stringr::str_detect(h.gsets.list[[i]], pattern = "\\+$|-$"))) {
          h.gsets.list.positive <- stringr::str_match(h.gsets.list[[i]],pattern = "(.+)\\+")[,2] %>% purrr::discard(is.na)
          h.gsets.list.negative <- stringr::str_match(h.gsets.list[[i]],pattern = "(.+)-")[,2] %>% purrr::discard(is.na)
          if (length(h.gsets.list.positive)==0) {
            singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,
                                                            upSet = h.gsets.list.negative,
                                                            centerScore = F)
          }
          if (length(h.gsets.list.negative)==0) {
            singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,
                                                            upSet = h.gsets.list.positive,
                                                            centerScore = F)
          }
          if ((length(h.gsets.list.positive)!=0)&(length(h.gsets.list.negative)!=0)) {
            singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,
                                                            upSet = h.gsets.list.positive,
                                                            downSet = h.gsets.list.negative,
                                                            centerScore = F)
          }

        }else{
          singscore.scores[[i]] <- singscore::simpleScore(singscore.rank,
                                                          upSet = h.gsets.list[[i]],
                                                          centerScore = F)
        }
        TotalScore <- NULL
        singscore.scores[[i]] <- singscore.scores[[i]] %>%
          dplyr::select(TotalScore) %>%
          magrittr::set_colnames(names(h.gsets.list)[i])
      }

      names(singscore.scores) <- names(h.gsets.list)
      singscore.scores <- do.call(cbind, singscore.scores)
      singscore.scores.list[[k]] <- singscore.scores
    }
    singscore.scores.list <- do.call(rbind, singscore.scores.list)
    object[["singscore"]] <- SeuratObject::CreateAssayObject(counts = t(singscore.scores.list))
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = t(singscore.scores.list), assay = "singscore")
    message("Finish calculate singscore scores")
    rm(singscore.rank)
    rm(singscore.scores)
    rm(singscore.scores.list)
    gc()

  }}, error = identity)


  #### 07.calculate ssgsea scores ####
  tryCatch({if ("ssgsea" %in% method) {
    message("Calculate ssgsea scores")
    h.gsets.list.ssgsea <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})

    ssgsea.scores.list <- list()
    for (k in seq_along(my.matrix.list)) {
      ssgsea.scores.list[[k]] <- GSVA::gsva(my.matrix[, my.matrix.list[[k]]],
                                            h.gsets.list.ssgsea,
                                            method = "ssgsea",
                                            kcdf = kcdf,
                                            ssgsea.norm = F,
                                            parallel.sz = ncores,
                                            verbose = F)

    }
    ssgsea.scores.list <- do.call(cbind, ssgsea.scores.list)
    object[["ssgsea"]] <- SeuratObject::CreateAssayObject(counts = ssgsea.scores.list)
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = as.matrix(ssgsea.scores.list), assay = "ssgsea")
    message("Finish calculate ssgsea scores")
    rm(ssgsea.scores.list)
    gc()

  }}, error = identity)


  #### 08.calculate JASMINE scores ####
  # ties.method: "average"
  tryCatch({if ("JASMINE" %in% method) {
    message("Calculate JASMINE scores")
    h.gsets.list.jasmine <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})

    # Calculating Mean Ranks for signature genes across each cell
    RankCalculation <- function(x,genes){
      subdata = x[x!=0]                                                                      ### Removing Dropouts from single cell
      DataRanksUpdated=rank(subdata)                                                         ### Calculating ranks of each signature gene per cell
      DataRanksSigGenes = DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)]        ### Shortling rank vector for signature genes
      CumSum = ifelse(length(DataRanksSigGenes),mean(DataRanksSigGenes,na.rm = TRUE),0 )     ### Calculating Mean of ranks for signature genes
      FinalRawRank = CumSum/length(subdata)                                                  ### Normalizing Means by total coverage
      return(FinalRawRank)
    }

    # Calculating enrichment of signature genes across each cell 	(using odds ratio)
    ORCalculation <- function(data, genes){
      GE = data[which(rownames(data) %in% genes),]                                          ### Subsetting data for signature genes
      NGE = data[-which(rownames(data) %in% genes),]                                        ### Subsetting data for non-signature genes
      SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))                                 ### Calculating Number of expressed Signature Genes per cell
      NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))                               ### Calculating Number of expressed Non-Signature Genes per cell
      SigGenesNE = nrow(GE) - SigGenesExp                                                   ### Calculating Number of Not expressed Signature Genes per cell
      SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)									  ### Replacing Zero's with 1
      NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)                                ### Replacing Zero's with 1
      NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)                               ### Calculating Number of Not expressed Non-Signature Genes per cell
      NSigGenesNE = NSigGenesNE - SigGenesNE
      OR = (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp)                         ### Calculating Enrichment (Odds Ratio)
      return(OR)
    }

    # Calculating enrichment of signature genes across each cell (using Likelihood ratio)
    LikelihoodCalculation <- function(data,genes){
      GE = data[which(rownames(data) %in% genes),]
      NGE = data[-which(rownames(data) %in% genes),]
      SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))
      NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))
      SigGenesNE = nrow(GE) - SigGenesExp
      SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)
      NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)
      NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)
      NSigGenesNE = NSigGenesNE - SigGenesNE
      LR1 = SigGenesExp*(NSigGenesExp + NSigGenesNE)
      LR2 = NSigGenesExp * (SigGenesExp + SigGenesNE)
      LR = LR1/LR2
      return(LR)
    }

    # Scalar [0,1] Normalization of Means and Enrichment across set of cells
    NormalizationJAS <- function(JAS_Scores){
      JAS_Scores = (JAS_Scores - min(JAS_Scores))/(max(JAS_Scores)- min(JAS_Scores))
      return(JAS_Scores)
    }

    # Signature Scoring via JASMINE mergining Means and Enrichment
    JASMINE <- function(data,genes,method){
      idx = match(genes,rownames(data))
      idx = idx[!is.na(idx)]
      if(length(idx)> 1){
        RM = apply(data,2,function(x) RankCalculation(x,genes))                              ### Mean RankCalculation for single cell data matrix
        RM = NormalizationJAS(RM)                                                            ### Normalizing Mean Ranks

        if(method == "oddsratio"){
          OR = ORCalculation(data,genes)			                                             ### Signature Enrichment Calculation for single cell data matrix (OR)
          OR = NormalizationJAS(OR)															 ### Normalizing Enrichment Scores (OR)
          JAS_Scores = (RM + OR)/2
        }else if(method == "likelihood"){

          LR = LikelihoodCalculation(data,genes)			                                     ### Signature Enrichment Calculation for single cell data matrix  (LR)
          LR = NormalizationJAS(LR)															 ### Normalizing Enrichment Scores (LR)
          JAS_Scores = (RM + LR)/2
        }
        FinalScores = data.frame(JAS_Scores)
        # FinalScores = data.frame(names(RM),JAS_Scores)                                       ### JASMINE scores
        # colnames(FinalScores)[1]='SampleID'
        return(FinalScores)
      }
    }



    # Either take a BPPARAM object, or make one on the spot using 'ncores'
    BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)

    # calculate
    jasmine.scores.list <- list()
    for (k in seq_along(my.matrix.list)) {
      jasmine.scores <- BiocParallel::bplapply(
        X = seq_along(h.gsets.list.jasmine),
        BPPARAM =  BPPARAM,
        FUN = function(x) {
          data.jasmine <- JASMINE(data = my.matrix[, my.matrix.list[[k]]],
                                  genes = h.gsets.list.jasmine[[x]],
                                  method = JASMINE.method)
          colnames(data.jasmine)[1] <- names(h.gsets.list.jasmine)[[x]]
          return(data.jasmine)
        })

      jasmine.scores <- do.call(cbind, jasmine.scores)
      jasmine.scores.list[[k]] <- jasmine.scores
    }

    jasmine.scores.list <- do.call(rbind, jasmine.scores.list)
    object[["JASMINE"]] <- SeuratObject::CreateAssayObject(counts = t(jasmine.scores.list))
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = t(jasmine.scores.list), assay = "JASMINE")
    message("Finish calculate jasmine scores")
    rm(jasmine.scores)
    rm(jasmine.scores.list)
    gc()
  }}, error = identity)



  #### 09.calculate VAM scores ####
  tryCatch({if ("VAM" %in% method) {
    message("Calculate VAM scores")
    h.gsets.list.vam <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})

    if ((assay %in% c("RNA", "SCT")) & (slot %in% c("data"))) {

      # install package from CRAN
      if (!requireNamespace("VAM", quietly = TRUE)) {
        message("install VAM package from CRAN")
        utils::install.packages("VAM", ask = F, update = F)
      }

      gene.set.collection <- VAM::createGeneSetCollection(gene.ids = rownames(object),
                                                          gene.set.collection = h.gsets.list.vam)
      SeuratObject::DefaultAssay(object) <- assay
      object2 <- VAM::vamForSeurat(seurat.data = object,
                                   gene.set.collection = gene.set.collection,
                                   center = F,
                                   gamma = T,
                                   sample.cov = F,
                                   return.dist = F)
      object[["VAM"]] <- SeuratObject::CreateAssayObject(counts = SeuratObject::GetAssayData(object2, assay = "VAMcdf", slot = "counts"))
      object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                           new.data = as.matrix(SeuratObject::GetAssayData(object2, assay = "VAMcdf", slot = "counts")),
                                           assay = "VAM")
      message("Finish calculate VAM scores")
      rm(object2)
      gc()



    }else{
      message("VAM needs normalized counts. And assay should be RNA or SCT, slot should be data.")
    }

  }}, error = identity)


  #### 10.single-cell-signature-explorer ####
  tryCatch({if ("scSE" %in% method) {
    message("Calculate scSE scores")
    h.gsets.list.scSE <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})


    # Either take a BPPARAM object, or make one on the spot using 'ncores'
    BPPARAM <- BiocParallel::MulticoreParam(workers=ncores)

    # calculate
    scSE.scores.list <- list()
    for (k in seq_along(my.matrix.list)) {
      my.matrix.subset <- my.matrix[, my.matrix.list[[k]]]
      umi.sum <- apply(my.matrix.subset, 2, sum)
      scSE.scores <- BiocParallel::bplapply(
        X = seq_along(h.gsets.list.scSE),
        BPPARAM =  BPPARAM,
        FUN = function(x) {

          umi.geneset <- apply(my.matrix.subset[intersect(h.gsets.list.scSE[[x]], rownames(my.matrix.subset)),],
                               2, sum)
          scse.result <- data.frame(umi.geneset/umi.sum*100,
                                    row.names = colnames(my.matrix.subset))
          colnames(scse.result) <- names(h.gsets.list.scSE)[x]
          return(scse.result)
        })

      scSE.scores <- do.call(cbind, scSE.scores)
      scSE.scores.list[[k]] <- scSE.scores
    }

    scSE.scores.list <- do.call(rbind, scSE.scores.list)

    object[["scSE"]] <- SeuratObject::CreateAssayObject(counts = t(scSE.scores.list))

    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = t(scSE.scores.list),
                                         assay = "scSE")
    message("Finish calculate scSE scores")
    rm(scSE.scores)
    rm(scSE.scores.list)
    gc()

  }}, error = identity)



  #### 11.calculate VISION scores ####
  tryCatch({if (("VISION" %in% method) & (slot %in% c("counts"))) {
    message("Calculate VISION scores")

    # install package from Github
    if (!requireNamespace("VISION", quietly = TRUE)) {
      message("install VISION package from Github")
      devtools::install_github("YosefLab/VISION", force =T)
    }

    # create gene set for VISION
    sigs <- list()
    for (i in seq_along(h.gsets.list)){
      if (any(stringr::str_detect(h.gsets.list[[i]], pattern = "\\+$|-$"))){

        sigData <- sapply(X = h.gsets.list[[i]],
                          FUN = function(x){
                            if (stringr::str_detect(x, pattern = "\\+")==T){
                              x=1
                            }else if (stringr::str_detect(x, pattern = "-")==T){
                              x=-1
                            }else{
                              x=1
                            }}, simplify = TRUE)
        names(sigData) <- stringr::str_remove(names(sigData), "\\+$|-$")

      }else{

        sigData <- structure(rep(1, length(h.gsets.list[[i]])), names = h.gsets.list[[i]])

      }
      sigs.name <- names(h.gsets.list)[i]
      sigs[[sigs.name]] <- VISION::createGeneSignature(name = sigs.name, sigData = sigData)

    }

    vision.obj <- VISION::Vision(object,
                                 signatures = sigs,
                                 projection_methods = NULL,
                                 assay = assay)
    vision.obj <- VISION::analyze(vision.obj)
    sigScores <- VISION::getSignatureScores(vision.obj)

    object[["VISION"]] <- SeuratObject::CreateAssayObject(counts = t(sigScores))
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = t(sigScores),
                                         assay = "VISION")
    message("Finish calculate VISION scores")
    rm(vision.obj)
    rm(sigScores)
    gc()


  }}, error = identity)




  #### 13.calculate decoupler scores ####
  tryCatch({if (any(c("wmean", "wsum", "mdt", "viper",  "GSVApy") %in% method)) {


    # install.package from Bioconductor
    if (!requireNamespace("decoupleR", quietly = TRUE)) {
      BiocManager::install("decoupleR", ask = F, update = F, force = T)
    }


    # create gene set for decoupler
    h.gsets.list.decoupler <- list()
    for (i in seq_along(h.gsets.list)) {

      if (any(stringr::str_detect(h.gsets.list[[i]], pattern = "\\+$|-$"))){

        if (is.null(geneset.weight)) {
          weightData <- sapply(X = h.gsets.list[[i]],
                               FUN = function(x){
                                 if (stringr::str_detect(x, pattern = "\\+")==T){
                                   x=1
                                 }else if (stringr::str_detect(x, pattern = "-")==T){
                                   x=-1
                                 }else{
                                   x=1
                                 }}, simplify = TRUE)
          h.gsets.list.decoupler[[i]] <- data.frame(source = names(h.gsets.list)[i],
                                                    target = h.gsets.list[[i]],
                                                    weight = weightData)
        }else{
          h.gsets.list.decoupler[[i]] <- data.frame(source = names(h.gsets.list)[i],
                                                    target = h.gsets.list[[i]],
                                                    weight = geneset.weight[h.gsets.list[[i]]],
                                                    row.names = NULL)
        }


      }else{


        if (is.null(geneset.weight)) {
          h.gsets.list.decoupler[[i]] <- data.frame(source = names(h.gsets.list)[i],
                                                    target = h.gsets.list[[i]],
                                                    weight = 1)
        }else{
          h.gsets.list.decoupler[[i]] <- data.frame(source = names(h.gsets.list)[i],
                                                    target = h.gsets.list[[i]],
                                                    weight = geneset.weight[h.gsets.list[[i]]],
                                                    row.names = NULL)
        }

      }


    }

    net <- do.call(rbind, h.gsets.list.decoupler)

    if ("wmean" %in% method) {
      message("Calculate wmean scores")
      # Run wmean
      acts <- decoupleR::run_wmean(mat = my.matrix,
                                   net = net,
                                   .source='source',
                                   .target='target',
                                   .mor='weight',
                                   times = 100,
                                   minsize = minGSSize,
                                   seed = seeds,
                                   sparse = T)
      source <- NULL
      condition <- NULL
      score <- NULL
      statistic <- NULL
      acts <- acts %>%
        dplyr::filter(statistic == "corr_wmean") %>%
        dplyr::select(c("source", "condition", "score")) %>%
        tidyr::pivot_wider(names_from = source, values_from = score) %>%
        tibble::column_to_rownames(var = "condition")

      object[["wmean"]] <- SeuratObject::CreateAssayObject(counts = t(acts))
      object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                           new.data = t(acts),
                                           assay = "wmean")
      message("Finish calculate wmean scores")
    }

    if ("wsum" %in% method) {
      message("Calculate wsum scores")
      # Run wsum
      acts <- decoupleR::run_wsum(mat = my.matrix,
                                  net = net,
                                  .source='source',
                                  .target='target',
                                  .mor='weight',
                                  times = 100,
                                  minsize = minGSSize,
                                  seed = seeds,
                                  sparse = T)
      source <- NULL
      condition <- NULL
      score <- NULL
      statistic <- NULL
      acts <- acts %>%
        dplyr::filter(statistic == "corr_wsum") %>%
        dplyr::select(c("source", "condition", "score")) %>%
        tidyr::pivot_wider(names_from = source, values_from = score) %>%
        tibble::column_to_rownames(var = "condition")

      object[["wsum"]] <- SeuratObject::CreateAssayObject(counts = t(acts))
      object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                           new.data = t(acts),
                                           assay = "wsum")
      message("Finish calculate wsum scores")
    }


    if ("mdt" %in% method) {
      message("Calculate mdt scores")
      # Run mdt
      acts <- decoupleR::run_mdt(mat = my.matrix,
                                 net = net,
                                 .source='source',
                                 .target='target',
                                 .mor='weight',
                                 minsize = minGSSize,
                                 nproc = ncores,
                                 seed = seeds,
                                 sparse = T)
      source <- NULL
      condition <- NULL
      score <- NULL
      acts <- acts %>%
        dplyr::select(c("source", "condition", "score")) %>%
        tidyr::pivot_wider(names_from = source, values_from = score) %>%
        tibble::column_to_rownames(var = "condition")

      object[["mdt"]] <- SeuratObject::CreateAssayObject(counts = t(acts))
      object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                           new.data = t(acts),
                                           assay = "mdt")
      message("Finish calculate mdt scores")

    }


    if ("viper" %in% method) {
      message("Calculate viper scores")
      # Run viper
      acts <- decoupleR::run_viper(mat = as.matrix(my.matrix),
                                   net = net,
                                   .source='source',
                                   .target='target',
                                   .mor='weight',
                                   minsize = minGSSize,
                                   cores = ncores)
      source <- NULL
      condition <- NULL
      score <- NULL
      acts <- acts %>%
        dplyr::select(c("source", "condition", "score")) %>%
        tidyr::pivot_wider(names_from = source, values_from = score) %>%
        tibble::column_to_rownames(var = "condition")

      object[["viper"]] <- SeuratObject::CreateAssayObject(counts = t(acts))
      object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                           new.data = t(acts),
                                           assay = "viper")
      message("Finish calculate viper scores")

    }

    tryCatch({if (any(c("GSVApy") %in% method)) {

      message("Calculate GSVApy scores")

      # install.package from CRAN
      if (!requireNamespace("reticulate", quietly = TRUE)) {
        utils::install.packages("reticulate", ask = F, update = F)
      }

      # install.package from Github
      if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
        devtools::install_github("mojaveazure/seurat-disk", force =T)
      }

      # reticulate::py_config()
      if (! "irGSEA" %in% reticulate::conda_list()$name) {
        reticulate::conda_create("irGSEA")
      }

      if (!reticulate::py_module_available(c("anndata", "scanpy", "decoupler", "argparse"))) {
        reticulate::conda_install(envname = "irGSEA",
                                  packages = c("anndata", "scanpy", "decoupler",
                                               "argparse"),
                                  pip = T)
      }

      # convert seurat to h5ad
      seurat2scanpy <- function(x){
        temp <- SeuratObject::CreateSeuratObject(counts = x, slot = "counts")
        SeuratDisk::SaveH5Seurat(temp, filename = "./temp.h5Seurat", overwrite = T)
        SeuratDisk::Convert("./temp.h5Seurat", dest = "h5ad", overwrite = T)

      }

      seurat2scanpy(my.matrix)

      readr::write_csv(net, "./geneset.csv")

      # prepare
      name <- NULL
      python <- NULL
      python = reticulate::conda_list() %>%
        dplyr::filter(name == "irGSEA") %>%
        dplyr::pull(python)

      # GSVApy
      if ("GSVApy" %in% method) {
        message("Calculate GSVApy scores")

        python.file = paste0(system.file(package = 'irGSEA'), '/python/gsva.py')
        if(file.exists(python.file)){
          python.file = python.file
        }else{
          python.file = paste0(system.file(package = 'irGSEA'), '/inst/python/gsva.py')
        }

        if (kcdf == 'Gaussian') {
          kcdf.py = stringr::str_c(c("--kcdf", "True"), collapse = " ")
        }else{
          kcdf.py = stringr::str_c(c("--kcdf", "False"), collapse = " ")
        }

        min_n.py = stringr::str_c(c("--min_n", minGSSize), collapse = " ")
        seed.py = stringr::str_c(c("--seed", seeds), collapse = " ")

        command = stringr::str_c(c(python, python.file, kcdf.py, min_n.py, seed.py), collapse = " ")
        message(command)
        # work
        system(command)

        acts <- readr::read_csv("./matrix.py.result.csv")
        colnames(acts)[1] <- "cell"
        acts <- acts %>% tibble::column_to_rownames(var = "cell")
        object[["GSVApy"]] <- SeuratObject::CreateAssayObject(counts = t(acts))
        object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                             new.data = t(acts),
                                             assay = "GSVApy")
        message("Finish calculate GSVApy scores")
      }
      file.remove("./temp.h5Seurat")
      file.remove("./temp.h5ad")
      file.remove("./geneset.csv")
      file.remove("./matrix.py.result.csv")

    }}, error = identity)


    rm(acts)
    gc()

  }}, error = identity)






  #### 14.calculate gficf scores ####
  # gficf = NMF+ssGSEA
  tryCatch({if ("gficf" %in% method) {
    message("Calculate gficf scores")
    h.gsets.list.gficf <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})

    # install.package from Bioconductor
    if (!requireNamespace(c("sva","edgeR", "fgsea"), quietly = TRUE)) {
      BiocManager::install(c("sva","edgeR", "fgsea"), ask = F, update = F)
    }

    if (!requireNamespace("gficf", quietly = TRUE)) {
      message("install gficf package from Github")
      devtools::install_github("gambalab/gficf", force =T)
    }

    # prepare data
    data <- gficf::gficf(M = my.matrix)
    data <- gficf::runPCA(data = data,dim = 10,use.odgenes = T)

    # custom function
    myrunScGSEA <- function (data, geneID, species, category, subcategory = NULL,
                             pathway.list = NULL, nsim = 10000, nt = 0, minSize = 15,
                             maxSize = Inf, verbose = TRUE, seed = 180582, nmf.k = 100,
                             fdr.th = 0.05, gp = 0, rescale = "none", normalization = "gficf")
    {
      if (nt == 0) {
        nt = parallel::detectCores()
      }
      geneID = base::match.arg(arg = geneID, choices = c("ensamble",
                                                         "symbol"), several.ok = F)
      rescale = base::match.arg(arg = rescale, choices = c("none",
                                                           "byGS", "byCell"), several.ok = F)
      normalization = base::match.arg(arg = normalization, choices = c("gficf",
                                                                       "cpm"), several.ok = F)
      options(RcppML.threads = nt)
      set.seed(seed)
      # def function
      tsmessage <-utils::getFromNamespace("tsmessage", "gficf")

      if (is.null(data$scgsea)) {
        data$scgsea = list()
        if (normalization == "gficf") {
          if (!is.null(data$pca) && data$pca$type == "NMF") {
            if (data$dimPCA < nmf.k || data$pca$use.odgenes) {
              tsmessage("... Performing NMF", verbose = verbose)
              tmp = RcppML::nmf(A = data$gficf, k = nmf.k)
              data$scgsea$nmf.w <- Matrix::Matrix(data = tmp$w,
                                                  sparse = T,
                                                  dimnames = list(rownames(data$gficf),NULL))
              data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp$h,
                                                    sparse = T,
                                                    dimnames = list(NULL,colnames(data$gficf))))
              rm(tmp)
              gc()
            }else {
              tsmessage(paste0("Found NMF reduction with k greaten or equal to ",
                                       nmf.k), verbose = T)
              pointr::ptr("tmp", "data$pca$genes")
              data$scgsea$nmf.w = tmp
              pointr::ptr("tmp2", "data$pca$cells")
              tmp2 <- NULL
              data$scgsea$nmf.h = tmp2
              rm(tmp, tmp2)
              gc()
            }
          }else {
            tsmessage("... Performing NMF", verbose = verbose)
            tmp = RcppML::nmf(data$gficf, k = nmf.k)
            data$scgsea$nmf.w <- Matrix::Matrix(data = tmp$w,
                                                sparse = T,
                                                dimnames = list(rownames(data$gficf),NULL))
            data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp$h,
                                                  sparse = T,
                                                  dimnames = list(NULL,colnames(data$gficf))))
            rm(tmp)
            gc()
          }
        }else {
          tsmessage("... Performing NMF", verbose = verbose)
          tmp = RcppML::nmf(log1p(data$rawCounts), k = nmf.k)
          data$scgsea$nmf.w <- Matrix::Matrix(data = tmp$w,
                                              sparse = T)
          data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp$h,
                                                sparse = T))
          rm(tmp)
          gc()
        }
      }else {
        tsmessage("Found a previous scGSEA, thus the already computed NMF will be used",
                          verbose = T)
        tsmessage("If you want to recompute NMF, please call resetScGSEA first",
                          verbose = T)
        data$scgsea$es <- NULL
        data$scgsea$nes <- NULL
        data$scgsea$pval <- NULL
        data$scgsea$fdr <- NULL
        data$scgsea$pathways <- NULL
        data$scgsea$x <- NULL
        data$scgsea$stat <- NULL
      }
      tsmessage("Loading pathways...", verbose = verbose)
      if (!is.list(pathway.list)) {
        gs = msigdbr::msigdbr(species = species, category = category,
                              subcategory = subcategory)
        if (geneID == "symbol") {
          data$scgsea$pathways = split(x = gs$gene_symbol,
                                       f = gs$gs_name)
        }else {
          data$scgsea$pathways = split(x = gs$ensembl_gene,
                                       f = gs$gs_name)
        }
      }else {
        data$scgsea$pathways = pathway.list
        rm(pathway.list)
      }
      data$scgsea$es = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways),
                                      ncol = ncol(data$scgsea$nmf.w))
      data$scgsea$nes = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways),
                                       ncol = ncol(data$scgsea$nmf.w))
      data$scgsea$pval = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways),
                                        ncol = ncol(data$scgsea$nmf.w))
      data$scgsea$fdr = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways),
                                       ncol = ncol(data$scgsea$nmf.w))
      rownames(data$scgsea$es) = rownames(data$scgsea$nes) = rownames(data$scgsea$pval) = rownames(data$scgsea$fdr) = names(data$scgsea$pathways)
      tsmessage("Performing GSEA...", verbose = verbose)
      oldw <- getOption("warn")
      options(warn = -1)
      pb = utils::txtProgressBar(min = 0, max = ncol(data$scgsea$nmf.w), initial = 0, style = 3)
      nt_fgsea <- ceiling(length(data$scgsea$pathways)/100)
      nt_fgsea <- ifelse(nt_fgsea > nt, nt, nt_fgsea)
      bpparameters <- BiocParallel::SnowParam(nt_fgsea)
      for (i in 1:ncol(data$scgsea$nmf.w)) {
        df = as.data.frame(fgsea::fgseaMultilevel(pathways = data$scgsea$pathways,
                                                  stats = data$scgsea$nmf.w[, i], nPermSimple = nsim,
                                                  gseaParam = gp, BPPARAM = bpparameters, minSize = minSize,
                                                  maxSize = maxSize))[, 1:7]
        data$scgsea$es[df$pathway, i] = df$ES
        data$scgsea$nes[df$pathway, i] = df$NES
        data$scgsea$pval[df$pathway, i] = df$pval
        data$scgsea$fdr[df$pathway, i] = df$padj
        utils::setTxtProgressBar(pb, i)
      }
      base::close(pb)
      on.exit(options(warn = oldw))
      ix = is.na(data$scgsea$nes)
      if (sum(ix) > 0) {
        data$scgsea$nes[ix] = 0
        data$scgsea$pval[ix] = 1
        data$scgsea$fdr[ix] = 1
      }
      data$scgsea$x = data$scgsea$nes
      data$scgsea$x[data$scgsea$x < 0 | data$scgsea$fdr >= fdr.th] = 0
      data$scgsea$x = Matrix::Matrix(data = data$scgsea$nmf.h %*% t(data$scgsea$x), sparse = T)
      data$scgsea$stat = df[, c("pathway", "size")]

      # def function
      armaColSum <-utils::getFromNamespace("armaColSum", "gficf")

      data$scgsea$x = data$scgsea$x[, armaColSum(data$scgsea$x) > 0]
      if (rescale != "none") {
        if (rescale == "byGS") {
          data$scgsea$x = t(data$scgsea$x)
          data$scgsea$x = t((data$scgsea$x - rowMeans(data$scgsea$x))/apply(data$scgsea$x, 1, stats::sd))
        }
        if (rescale == "byCell") {
          data$scgsea$x = (data$scgsea$x - rowMeans(data$scgsea$x))/apply(data$scgsea$x, 1, stats::sd)
        }
      }
      return(data)
    }


    if (ncol(my.matrix) > 10000) {
      nmf.k <- 100
    }else{
      nmf.k <- 50
    }

    # perform
    data <- myrunScGSEA(data = data,
                        geneID = "symbol",
                        minSize = minGSSize,
                        pathway.list = h.gsets.list.gficf,
                        seed = seeds,
                        nmf.k = nmf.k,
                        fdr.th = 1,
                        rescale = "none",
                        verbose = T,
                        nt = ncores)

    object[["gficf"]] <- SeuratObject::CreateAssayObject(counts = t(data$scgsea$x))
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = as.matrix(t(data$scgsea$x)),
                                         assay = "gficf")
    message("Finish calculate gficf scores")
    rm(data)
    gc()


  }}, error = identity)


  #### 15.calculate GSVA scores ####
  tryCatch({if ("GSVA" %in% method) {
    message("Calculate GSVA scores")
    h.gsets.list.GSVA <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})

    GSVA.scores.list <- list()
    GSVA.scores.list[[1]] <- GSVA::gsva(my.matrix,
                                        h.gsets.list.GSVA,
                                        method = "gsva",
                                        kcdf = kcdf,
                                        parallel.sz = ncores,
                                        verbose = F)
    GSVA.scores.list <- do.call(cbind, GSVA.scores.list)
    object[["GSVA"]] <- SeuratObject::CreateAssayObject(counts = GSVA.scores.list)
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = as.matrix(GSVA.scores.list), assay = "GSVA")
    message("Finish calculate GSVA scores")
    rm(GSVA.scores.list)
    gc()

  }}, error = identity)


  #### 16.calculate zscore scores ####
  tryCatch({if ("zscore" %in% method) {
    message("Calculate zscore scores")
    h.gsets.list.zscore <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})

    zscore.scores.list <- list()
    zscore.scores.list[[1]] <- GSVA::gsva(as.matrix(my.matrix),
                                          h.gsets.list.zscore,
                                          method = "zscore",
                                          kcdf = kcdf,
                                          parallel.sz = ncores,
                                          verbose = F)

    zscore.scores.list <- do.call(cbind, zscore.scores.list)
    object[["zscore"]] <- SeuratObject::CreateAssayObject(counts = zscore.scores.list)
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = as.matrix(zscore.scores.list),
                                         assay = "zscore")
    message("Finish calculate zscore scores")
    rm(zscore.scores.list)
    gc()

  }}, error = identity)


  #### 17.calculate plage scores ####
  tryCatch({if ("plage" %in% method) {
    message("Calculate plage scores")
    h.gsets.list.plage <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})

    plage.scores.list <- list()
    plage.scores.list[[1]] <- GSVA::gsva(as.matrix(my.matrix),
                                         h.gsets.list.plage,
                                         method = "plage",
                                         kcdf = kcdf,
                                         parallel.sz = ncores,
                                         verbose = F)

    plage.scores.list <- do.call(cbind, plage.scores.list)
    object[["plage"]] <- SeuratObject::CreateAssayObject(counts = plage.scores.list)
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = as.matrix(plage.scores.list), assay = "plage")
    message("Finish calculate plage scores")
    rm(plage.scores.list)
    gc()

  }}, error = identity)





  #### 18.calculate ssgsea.py scores ####

  # if ("ssGSEApy" %in% method) {
  #   message("Calculate ssGSEApy scores")
  #   h.gsets.list.ssgsea.py <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})
  #
  #   # install.package from CRAN
  #   if (!requireNamespace("reticulate", quietly = TRUE)) {
  #     install.packages("reticulate", ask = F, update = F)
  #   }
  #
  #   # install.package from Github
  #   if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  #     devtools::install_github("mojaveazure/seurat-disk", force =T)
  #   }
  #
  #   # reticulate::py_config()
  #   if (! "irGSEA" %in% reticulate::conda_list()$name) {
  #     reticulate::conda_create("irGSEA")
  #   }
  #   if (! reticulate::py_module_available(c("anndata", "scanpy", "gseapy", "argparse"))) {
  #     reticulate::conda_install(envname = "irGSEA",
  #                               packages = c("anndata", "scanpy", "gseapy",
  #                                            "argparse"),
  #                               pip = T)
  #   }
  #
  #
  #   # convert seurat to h5ad
  #   seurat2scanpy <- function(object, slot, assay){
  #     if (slot == "counts") {
  #       temp <- Seurat::DietSeurat(object, counts = TRUE, data = FALSE,
  #                                  scale.data = FALSE,
  #                                  assays = assay)
  #     }
  #     if (slot == "data") {
  #       temp <- Seurat::DietSeurat(object, counts = FALSE, data = TRUE,
  #                                  scale.data = FALSE,
  #                                  assays = assay)
  #     }
  #     if (slot == "scale.data") {
  #       temp <- Seurat::DietSeurat(object, counts = FALSE, data = FALSE,
  #                                  scale.data = TRUE,
  #                                  assays = assay)
  #     }
  #
  #
  #     SeuratDisk::SaveH5Seurat(temp, filename = "./temp.h5Seurat", overwrite = T)
  #     SeuratDisk::Convert("./temp.h5Seurat", dest = "h5ad", overwrite = T)
  #
  #   }
  #
  #   seurat2scanpy(object, assay, slot)
  #
  #
  #   for (i in seq_along(h.gsets.list.ssgsea.py)) {
  #     h.gsets.list.ssgsea.py[[i]] <- data.frame(source = names(h.gsets.list.ssgsea.py)[i],
  #                                                  target = h.gsets.list.ssgsea.py[[i]])
  #   }
  #   net <- do.call(rbind, h.gsets.list.ssgsea.py)
  #   readr::write_csv(net, "./geneset.csv")
  #
  #   # write gmt file
  #   write.mygmt <- function(geneSet = h.gsets.list.ssgsea.py,  gmt_file ='./geneset.gmt'){
  #     sink(gmt_file)
  #     for (i in 1:length(geneSet)){
  #       cat(names(geneSet)[i])
  #       cat('\ttemp\t')
  #       cat(paste(geneSet[[i]],collapse = '\t'))
  #       cat('\n')
  #     }
  #     sink()
  #   }
  #
  #   write.mygmt(h.gsets.list.ssgsea.py, './geneset.gmt')
  #
  #   python.file = paste0(system.file(package = 'irGSEA'), '/python/ssgsea.py')
  #   if(file.exists(python.file)){
  #     python.file = python.file
  #   }else{
  #     python.file = paste0(system.file(package = 'irGSEA'), '/inst/python/ssgsea.py')
  #   }
  #
  #
  #
  #   min_size.py = stringr::str_c(c("--min_size", minGSSize), collapse = " ")
  #   max_size.py = stringr::str_c(c("--max_size", maxGSSize), collapse = " ")
  #   seed.py = stringr::str_c(c("--seed", seeds), collapse = " ")
  #   threads.py = stringr::str_c(c("--threads", ncores), collapse = " ")
  #
  #   command = stringr::str_c(c(python, python.file, min_size.py, max_size.py, seed.py, threads.py), collapse = " ")
  #   message(command)
  #   # work
  #   system(command)
  #
  #   acts <- readr::read_csv("./matrix.py.result.csv")
  #   colnames(acts)[1] <- "cell"
  #   acts <- acts %>% tibble::column_to_rownames(var = "cell")
  #
  #
  #   object[["ssGSEApy"]] <- SeuratObject::CreateAssayObject(counts = act)
  #   message("Finish calculate ssGSEApy scores")
  #   file.remove("./temp.h5Seurat")
  #   file.remove("./temp.h5ad")
  #   file.remove("./geneset.csv")
  #   file.remove("./matrix.py.result.csv")
  #   rm(acts)
  #   gc()
  #
  #
  # }
  #
  #

  #### 18.calculate AddModuleScore ####
  tryCatch({if (("AddModuleScore" %in% method) & (slot %in% c("data"))) {
    message("Calculate AddModuleScore scores")
    h.gsets.list.AddModuleScore <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})

    object.AddModuleScore <- Seurat::AddModuleScore(object,
                                                    features = h.gsets.list.AddModuleScore,
                                                    name = names(h.gsets.list.AddModuleScore),
                                                    seed = seeds,
                                                    assay = assay)
    object.AddModuleScore <- object.AddModuleScore@meta.data

    path_names = names(h.gsets.list.AddModuleScore)
    score = matrix(NA,nrow=length(path_names),ncol=ncol(my.matrix))
    rownames(score) = path_names
    colnames(score) = colnames(my.matrix)

    for (i in rownames(score)) {
      score[i,] = as.numeric(object.AddModuleScore[, stringr::str_detect(colnames(object.AddModuleScore), pattern = i)])
    }


    object[["AddModuleScore"]] <- SeuratObject::CreateAssayObject(counts = score)
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = score,
                                         assay = "AddModuleScore")
    message("Finish calculate AddModuleScore scores")
    rm(score)
    rm(object.AddModuleScore)
    gc()


  }}, error = identity)





  #### 19.calculate pagoda2 scores ####
  tryCatch({if (("pagoda2" %in% method) & (slot %in% c("counts"))) {
    message("Calculate pagoda2 scores")
    h.gsets.list.pagoda2 <- h.gsets.list %>% purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})
    # install.package from CRAN
    if (!requireNamespace("pagoda2", quietly = TRUE)) {
      utils::install.packages("pagoda2", ask = F, update = F)
    }
    # install.package from Bioconductor
    if (!requireNamespace("scde", quietly = TRUE)) {
      BiocManager::install("scde", ask = F, update = F, force = T)
    }
    # devtools::install_github("cran/flexmix@2.3-13")
    data.pagoda2 <- pagoda2::basicP2proc(my.matrix,
                                         n.cores = ncores,
                                         n.odgenes=2e3,
                                         get.largevis=FALSE,
                                         make.geneknn=FALSE,
                                         get.tsne = F)
    h.gsets.list.pagoda2.env <- list2env(h.gsets.list.pagoda2)
    mytestPathwayOverdispersion=function(setenv, type='counts', max.pathway.size=1e3, min.pathway.size=10,
                                         n.randomizations=5, verbose=FALSE, n.cores=self$n.cores, score.alpha=0.05, plot=FALSE, cells=NULL, adjusted.pvalues=TRUE,
                                         z.score = stats::qnorm(0.05/2, lower.tail = FALSE), use.oe.scale = FALSE, return.table=FALSE, name='pathwayPCA',
                                         correlation.distance.threshold=0.2, loading.distance.threshold=0.01, top.aspects=Inf, recalculate.pca=FALSE, save.pca=TRUE) {

      if (!requireNamespace("scde", quietly=TRUE)){
        stop("You need to install package 'scde' to be able to use testPathwayOverdispersion().")
      }

      nPcs <- 1
      if (type=='counts') {
        x <- self$counts
        # apply scaling if using raw counts
        x@x <- x@x*rep(self$misc[['varinfo']][colnames(x),'gsf'],diff(x@p))
      } else {
        if (!type %in% names(self$reductions)) { stop("Reduction ",type,' not found')}
        x <- self$reductions[[type]]
      }
      if (!is.null(cells)) {
        x <- x[cells,]
      }

      proper.gene.names <- colnames(x)

      if (is.null(self$misc[['pwpca']]) || recalculate.pca) {
        if (verbose) {
          message("determining valid pathways")
        }

        # def function
        sn <- function(x){
          names(x) <- x
          return(x)
        }


        # determine valid pathways
        gsl <- ls(envir = setenv)
        gsl.ng <- unlist(parallel::mclapply(sn(gsl), function(go) sum(unique(get(go, envir = setenv)) %in% proper.gene.names),mc.cores=n.cores,mc.preschedule=TRUE))
        gsl <- gsl[gsl.ng >= min.pathway.size & gsl.ng<= max.pathway.size]
        names(gsl) <- gsl

        if (verbose) {
          message("processing ", length(gsl), " valid pathways")
        }

        cm <- Matrix::colMeans(x)
        # def function
        papply <-utils::getFromNamespace("papply", "pagoda2")
        pwpca <- papply(gsl, function(sn) {
          lab <- proper.gene.names %in% get(sn, envir = setenv)
          if (sum(lab)<1) {
            return(NULL)
          }
          pcs <- irlba::irlba(x[,lab], nv=nPcs, nu=0, center=cm[lab])
          pcs$d <- pcs$d/sqrt(nrow(x))
          pcs$rotation <- pcs$v
          pcs$v <- NULL

          # get standard deviations for the random samples
          ngenes <- sum(lab)
          z <- do.call(rbind,lapply(seq_len(n.randomizations), function(i) {
            si <- sample(ncol(x), ngenes)
            pcs <- irlba::irlba(x[,si], nv=nPcs, nu=0, center=cm[si])$d
          }))
          z <- z/sqrt(nrow(x))

          # local normalization of each component relative to sampled PC1 sd
          avar <- pmax(0, (pcs$d^2-mean(z[, 1]^2))/stats::sd(z[, 1]^2))

          if (T) {
            # flip orientations to roughly correspond with the means
            pcs$scores <- as.matrix(t(x[,lab] %*% pcs$rotation) - as.numeric((cm[lab] %*% pcs$rotation)))
            cs <- unlist(lapply(seq_len(nrow(pcs$scores)), function(i) sign(stats::cor(pcs$scores[i,], colMeans(t(x[, lab, drop = FALSE])*abs(pcs$rotation[, i]))))))
            pcs$scores <- pcs$scores*cs
            pcs$rotation <- pcs$rotation*cs
            rownames(pcs$rotation) <- colnames(x)[lab]
          } # don't bother otherwise - it's not significant
          return(list(xp=pcs,z=z,n=ngenes))
        }, n.cores = n.cores,mc.preschedule=TRUE)
        if (save.pca) {
          self$misc[['pwpca']] <- pwpca
        }
      } else {
        if (verbose) {
          message("reusing previous overdispersion calculations")
          pwpca <- self$misc[['pwpca']]
        }
      }

      if (verbose) {
        message("scoring pathway od signifcance")
      }

      # score overdispersion
      true.n.cells <- nrow(x)

      pagoda.effective.cells <- function(pwpca, start = NULL) {
        n.genes <- unlist(lapply(pwpca, function(x) rep(x$n, nrow(x$z))))
        var <- unlist(lapply(pwpca, function(x) x$z[, 1]))
        if (is.null(start)) { start <- true.n.cells*2 } # start with a high value
        of <- function(p, v, sp) {
          sn <- p[1]
          vfit <- (sn+sp)^2/(sn*sn+1/2) -1.2065335745820*(sn+sp)*((1/sn + 1/sp)^(1/3))/(sn*sn+1/2)
          residuals <- (v-vfit)^2
          return(sum(residuals))
        }
        x <- stats::nlminb(objective = of, start = c(start), v = var, sp = sqrt(n.genes-1/2), lower = c(1), upper = c(true.n.cells))
        return((x$par)^2+1/2)
      }
      n.cells <- pagoda.effective.cells(pwpca)

      vdf <- data.frame(do.call(rbind, lapply(seq_along(pwpca), function(i) {
        vars <- as.numeric((pwpca[[i]]$xp$d))
        cbind(i = i, var = vars, n = pwpca[[i]]$n, npc = seq(1:ncol(pwpca[[i]]$xp$rotation)))
      })))

      # fix p-to-q mistake in qWishartSpike
      qWishartSpikeFixed <- function (q, spike, ndf = NA, pdim = NA, var = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)  {
        params <- RMTstat::WishartSpikePar(spike, ndf, pdim, var, beta)
        stats::qnorm(q, mean = params$centering, sd = params$scaling, lower.tail, log.p)
      }

      # add right tail approximation to ptw, which gives up quite early
      pWishartMaxFixed <- function (q, ndf, pdim, var = 1, beta = 1, lower.tail = TRUE) {
        params <- RMTstat::WishartMaxPar(ndf, pdim, var, beta)
        q.tw <- (q - params$centering)/(params$scaling)
        p <- RMTstat::ptw(q.tw, beta, lower.tail, log.p = TRUE)
        p[p == -Inf] <- stats::pgamma((2/3)*q.tw[p == -Inf]^(3/2), 2/3, lower.tail = FALSE, log.p = TRUE) + lgamma(2/3) + log((2/3)^(1/3))
        p
      }

      vshift <- 0
      ev <- 0

      vdf$var <- vdf$var-(vshift-ev)*vdf$n
      basevar <- 1
      vdf$exp <- RMTstat::qWishartMax(0.5, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
      #vdf$z <- qnorm(pWishartMax(vdf$var, n.cells, vdf$n, log.p = TRUE, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
      vdf$z <- stats::qnorm(pWishartMaxFixed(vdf$var, n.cells, vdf$n, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
      # def function
      bh.adjust <-utils::getFromNamespace("bh.adjust", "pagoda2")
      vdf$cz <- stats::qnorm(bh.adjust(stats::pnorm(as.numeric(vdf$z), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)
      vdf$ub <- RMTstat::qWishartMax(score.alpha/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
      vdf$ub.stringent <- RMTstat::qWishartMax(score.alpha/nrow(vdf)/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)

      if (plot) {
        test_pathway_par <- graphics::par(mfrow = c(1, 1), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
        on.exit(graphics::par(test_pathway_par))
        un <- sort(unique(vdf$n))
        on <- order(vdf$n, decreasing = FALSE)
        pccol <- grDevices::colorRampPalette(c("black", "grey70"), space = "Lab")(max(vdf$npc))
        plot(vdf$n, vdf$var/vdf$n, xlab = "gene set size", ylab = "PC1 var/n", ylim = c(0, max(vdf$var/vdf$n)), col = grDevices::adjustcolor(pccol[vdf$npc],alpha=0.1),pch=19)
        graphics::lines(vdf$n[on], (vdf$exp/vdf$n)[on], col = 2, lty = 1)
        graphics::lines(vdf$n[on], (vdf$ub.stringent/vdf$n)[on], col = 2, lty = 2)
      }

      rs <- (vshift-ev)*vdf$n
      vdf$oe <- (vdf$var+rs)/(vdf$exp+rs)
      vdf$oec <- (vdf$var+rs)/(vdf$ub+rs)

      df <- data.frame(name = names(pwpca)[vdf$i], npc = vdf$npc, n = vdf$n, score = vdf$oe, z = vdf$z, adj.z = vdf$cz, stringsAsFactors = FALSE)
      if (adjusted.pvalues) {
        vdf$valid <- vdf$cz  >=  z.score
      } else {
        vdf$valid <- vdf$z  >=  z.score
      }

      if (!any(vdf$valid)) {
        stop("No significantly overdispersed pathways found at z.score threshold of ",z.score)
      }

      # apply additional filtering based on >0.5 sd above the local random estimate
      vdf$valid <- vdf$valid & unlist(lapply(pwpca,function(x) !is.null(x$xp$scores)))
      vdf$name <- names(pwpca)[vdf$i]

      if (return.table) {
        df <- df[vdf$valid, ]
        df <- df[order(df$score, decreasing = TRUE), ]
        return(df)
      }
      if (verbose) {
        message("compiling pathway reduction")
      }
      # calculate pathway reduction matrix

      # return scaled patterns
      xmv <- do.call(rbind, lapply(pwpca[vdf$valid], function(x) {
        xm <- x$xp$scores
      }))

      if (use.oe.scale) {
        xmv <- (xmv -rowMeans(xmv))* (as.numeric(vdf$oe[vdf$valid])/sqrt(apply(xmv, 1, stats::var)))
        vdf$sd <- as.numeric(vdf$oe)
      } else {
        # chi-squared
        xmv <- (xmv-rowMeans(xmv)) * sqrt((stats::qchisq(stats::pnorm(vdf$z[vdf$valid], lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells)/apply(xmv, 1, stats::var))
        vdf$sd <- sqrt((stats::qchisq(stats::pnorm(vdf$z, lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells))

      }
      rownames(xmv) <- paste("#PC", vdf$npc[vdf$valid], "# ", names(pwpca)[vdf$i[vdf$valid]], sep = "")
      rownames(vdf) <- paste("#PC", vdf$npc, "# ", vdf$name, sep = "")
      self$misc[['pathwayODInfo']] <- vdf

      # collapse gene loading
      if (verbose) {
        message("clustering aspects based on gene loading ... ",appendLF=FALSE)
      }
      tam2 <- pagoda2::pagoda.reduce.loading.redundancy(list(xv=xmv,xvw=matrix(1,ncol=ncol(xmv),nrow=nrow(xmv))),pwpca,NULL,plot=FALSE,distance.threshold=loading.distance.threshold,n.cores=n.cores)
      if (verbose) {
        message(nrow(tam2$xv)," aspects remaining")
      }
      if (verbose) {
        message("clustering aspects based on pattern similarity ... ",appendLF=FALSE)
      }
      tam3 <- pagoda2::pagoda.reduce.redundancy(tam2, distance.threshold=correlation.distance.threshold,top=top.aspects)
      if (verbose) {
        message(nrow(tam3$xv)," aspects remaining\n")
      }
      tam2$xvw <- tam3$xvw <- NULL # to save space
      tam3$env <- setenv

      # clean up aspect names, as GO ids are meaningless
      names(tam3$cnam) <- rownames(tam3$xv) <- paste0('aspect',1:nrow(tam3$xv))

      self$misc[['pathwayOD']] <- tam3
      self$reductions[[name]] <- tam3$xv
      invisible(tam3)
    }

    data.pagoda2$mytestPathwayOverdispersion <- mytestPathwayOverdispersion
    environment(data.pagoda2$mytestPathwayOverdispersion) <- data.pagoda2$.__enclos_env__
    data.pagoda2$mytestPathwayOverdispersion(h.gsets.list.pagoda2.env,
                                             verbose=T,
                                             recalculate.pca=F,
                                             top.aspects=15, score.alpha =0)

    path_names = names(h.gsets.list)
    score = matrix(NA,nrow=length(path_names),ncol=ncol(my.matrix))
    rownames(score) = path_names
    colnames(score) = colnames(my.matrix)

    for(i in 1:length(data.pagoda2$misc$pwpca)){
      if(!is.null(data.pagoda2$misc$pwpca[[i]]$xp$score)){
        score[i,] = as.numeric(data.pagoda2$misc$pwpca[[i]]$xp$scores)
      }
    }

    object[["pagoda2"]] <- SeuratObject::CreateAssayObject(counts = score)
    object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                         new.data = score,
                                         assay = "pagoda2")
    message("Finish calculate pagoda2 scores")
    rm(score)
    rm(data.pagoda2)
    gc()


  }}, error = identity)



  return(object)

}


