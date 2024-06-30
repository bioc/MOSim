#' @import dplyr
#' @import Seurat
#' @import Signac
#' @import stringr
NULL

# Avoid harmless note with R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "n",
                                              "feature", "seurat_annotations"))


#' sc_omicData
#'
#' Checks if the user defined data is in the correct format, or loads the default
#' multiomics pbmc dataset, a subset from SeuratData package
#'
#' @param omics_types A list of strings which can be either "scRNA-seq" or "scATAC-seq"
#' @param data A user input matrix with genes (peaks in case of scATAC-seq) as 
#'    rows and cells as columns. By default, it loads the example data.
#'    If a user input matrix is included, cell columns must be sorted by cell t
#'    ype.
#' @return a named list with omics type as name and the count matrix as value
#' @export
#'
#' @examples
#' # Simulate from PBMC
#' omicsList <- sc_omicData(list("scRNA-seq", "scATAC-seq"))
#'
sc_omicData <- function(omics_types, data = NULL){
  # Check for mandatory parameters
  if (missing(omics_types)){
    stop("You must provide the vector of omics to simulate.")
  }
  
  for (omics in omics_types){
    if (omics != "scRNA-seq" && omics != "scATAC-seq"){
      stop("Omics must be a either 'scRNA-seq' or 'scATAC-seq'")
    }
  }
  
  omics_list <- list()
  # If no data, load default
  if (is.null(data)) {
    cell_types <- list(
      'Treg' = c(1:11),
      'cDC' = c(12:22),
      'CD4_TEM' = c(13:33),
      'Memory_B' = c(14:44)
    )
    message("Loading the default dataset, the cell_types are: 
            list('Treg' = c(1:10),'cDC' = c(11:20),'CD4_TEM' = c(21:30),
      'Memory_B' = c(31:40))")
    for (omics in omics_types){
      if (omics == "scRNA-seq"){
        data_env <- new.env(parent = emptyenv())
        data("scrna", envir = data_env, package = "MOSim")
        counts <- data_env[["scrna"]]
      } else if (omics == "scATAC-seq"){
        data_env <- new.env(parent = emptyenv())
        data("scatac", envir = data_env, package = "MOSim")
        counts <- data_env[["scatac"]]
      }
      omics_list[[omics]] <- counts
    }
    
  } else {
    # If data was inputted by the user, first check
    if (!is.list(data) || length(data) != 1 && length(data) != 2){
      message(paste0("The length of data is ", length(data)))
      stop("Data must be NULL (default) or a list of 1 or 2 elements")
    }
    
    N_data <- length(data)
    for (i in 1:N_data){
      # Then save in a named list
      if (!is.matrix(data[[i]]) && class(data[[i]])[1] != "Assay"){
        stop("Each element of data must be either a matrix or a Seurat object")
      } else if (is.matrix(data[[i]])){
        omics_list[[omics_types[[i]]]] <- data[[i]]
      } else if ("Assay" %in% class(data[[i]]) && omics_types[[i]] == "scRNA-seq"){
        counts <- as.matrix(data[[i]]@counts)
        omics_list[[omics_types[[i]]]] <- counts
      } else if ("Assay" %in% class(data[[i]]) && omics_types[[i]] == "scATAC-seq"){
        counts <- as.matrix(data[[i]]@counts)
        omics_list[[omics_types[[i]]]] <- counts
      }
    }
  }
  return(omics_list)
}
  

#' sc_param_estimation
#'
#' Evaluate the users parameters for single cell simulation and use SPARSim
#' to simulate the main dataset. Internal function
#'
#' @param omics named list containing the omics to simulate as names, which can 
#'    be "scRNA-seq" or "scATAC-seq".
#' @param cellTypes list where the i-th element of the list contains the column 
#'    indices for i-th cell type. List must be a named list.
#' @param diffGenes If number groups > 1, Percentage DE genes to simulate.
#'    List of vectors (one per group to compare to group 1) where the vector
#'    contains absolute number of genes for Up and Down ex: c(250, 500) or a 
#'    percentage for up, down ex: c(0.2, 0.2). The rest will be NE
#' @param minFC Threshold of FC below which are downregulated, by default 0.25
#' @param maxFC Threshold of FC above which are upregulated, by default 4
#' @param numberCells vector of numbers. The numbers correspond to the number 
#'    of cells the user wants to simulate per each cell type. The length of the 
#'       vector must be the same as length of \code{cellTypes}.
#' @param mean vector of numbers of mean depth per each cell type. Must be specified 
#'    just if \code{numberCells} is specified.
#' @param sd vector of numbers of standard deviation per each cell type. Must be 
#'    specified just if \code{numberCells} is specified.
#' @param noiseGroup OPTIONAL. Number indicating the desired standard deviation
#'      between treatment groups
#' @param group Group for which to estimate parameters
#' @param genereggroup List with information of genes, clusters and regulators
#'      that must be related to each other
#' @return a list of Seurat object, one per each omic.
#' @return a named list with simulation parameters for each omics as values.
#' @export
#' @examples
#' omicsList <- sc_omicData(list("scRNA-seq"))
#' cell_types <- list('Treg' = c(1:10),'cDC' = c(11:20),'CD4_TEM' = c(21:30),
#' 'Memory_B' = c(31:40))
#' #estimated_params <- sc_param_estimation(omicsList, cell_types)
#' 
sc_param_estimation <- function(omics, cellTypes, diffGenes = list(c(0.2, 0.2)), 
                                minFC = 0.25, maxFC = 4, numberCells = NULL, 
                                mean = NULL, sd = NULL, noiseGroup = 0.5, 
                                group = 1, genereggroup){
  # Check for mandatory parameters
  if (missing(omics)){
    stop("You must provide the vector of omics to simulate.")
  }
  
  if (missing(cellTypes)) {
    stop("You must provide the correspondence of cells and celltypes")
  }
  
  all_missing <- is.null(numberCells) && is.null(mean) && is.null(sd)
  all_specified <- !is.null(numberCells) && !is.null(mean) && !is.null(sd)
  only_cellnum <- !is.null(numberCells) && is.null(mean) && is.null(sd)
  
  if( !(all_missing || all_specified || only_cellnum)){
    
    stop("The user must either not provide the optional arguments, provide them 
         all or only provide cell numbers")
    
  }
  
  N_omics <- length(omics)
  
  ## Make association dataframe
  if (N_omics > 1){
    
    prov <- MOSim::make_association_dataframe(group, genereggroup, 
                                              dim(omics$`scRNA-seq`)[[1]], 
                                              dim(omics$`scATAC-seq`)[[1]],
                                              minFC, maxFC)
    associationMatrix <- prov$associationMatrix
    
    rownames(omics[[1]]) <- prov$dfGeneNames
    rownames(omics[[2]]) <- prov$dfPeakNames
    
    dfGeneNames <- prov$dfGeneNames
    dfPeakNames <- prov$dfPeakNames
    
    ## Also improve the association gene-regulator
    res <- MOSim::match_gene_regulator_cluster(omics[[1]], omics[[2]], cellTypes, associationMatrix)

    omics[[1]] <- res$rna
    omics[[2]] <- res$atac
    
  } else {
    # If only scRNA define the matrix here with the same format
    dfGeneNames <- rownames(omics[[1]])
    dfPeakNames <- NA
    columns <- c("Gene_ID", "Peak_ID", "RegulatorEffect", "Gene_cluster", "Peak_cluster", 
                 "Gene_DE", "Peak_DE")
    
    associationMatrix <- data.frame(matrix(nrow = length(dfGeneNames), ncol = length(columns)))
    colnames(associationMatrix) <- columns
    
    associationMatrix["Gene_ID"] <- dfGeneNames
    associationMatrix["Peak_ID"] <- rep(NA, length(dfGeneNames))
    associationMatrix["RegulatorEffect"] <- rep("NE", length(dfGeneNames))
    clus <- rep(1:length(genereggroup$`Clusters_scRNA-seq`), 
                each = length(genereggroup$`Clusters_scRNA-seq`[[1]]))
    associationMatrix["Gene_cluster"] <- c(clus, rep(0, length(dfGeneNames) - length(clus)))
    associationMatrix["Peak_cluster"] <- rep(NA, length(dfGeneNames))
    associationMatrix["Gene_DE"] <- c(rep("Up", length(genereggroup[[paste0("GeneExtraUp_G", group)]])),
                                      rep("Down", length(genereggroup[[paste0("GeneExtraDown_G", group)]])),
                                      rep("NE", length(dfGeneNames) - length(genereggroup[[paste0("GeneExtraUp_G", group)]]) - length(genereggroup[[paste0("GeneExtraDown_G", group)]])))
    associationMatrix$Gene_FC <- 1
    associationMatrix$Peak_FC <- NA
  }
  
  # Normalize using scran method, we had to suppress warnings because
  # ATAC has too many zeroes for the normalization to be super comfortable
  norm <- function(om) {
    o <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(om)))
    o <- suppressWarnings(scran::computeSumFactors(o, sizes = seq(20, 100, 5),
                                                   positive = FALSE))
    # Apply normalization factors
    o <- scater::normalizeCounts(o, log = FALSE)
    o[is.na(o)] <- 0
    o[!is.finite(o)] <- 0
    o <- abs(o)
    
    return(o)
  }
  
  norm_list <- lapply(omics, norm)
  param_est_list <- list()
  FC_used_list <- list()
  
  if (group > 1){
    # if its from group 2 upwards, we make the differences by multiplying
    # Times a fold change vector, thus we generate it
    FClist <- list()
    VARlist <- list()
    
    #We cant move this out of the function because it takes the variability
    # of the group into account
    for(i in 1:N_omics){
      message(paste0("Estimating distribution from original data type: ", i))
      param_est <- MOSim::sparsim_estimate_parameter_from_data(raw_data = omics[[i]],
                                                               norm_data = norm_list[[i]],
                                                               conditions = cellTypes)
      param_est_list[[paste0("param_est_", names(omics)[i])]] <- param_est
      
      
      #### Make the variability vector
      VARlist[[paste0("Var_", names(omics)[i])]] <- data.frame(matrix(0, 
             ncol = length(colnames(omics[[i]])), nrow = length(rownames(omics[[i]]))))
      var_group <- stats::rnorm(nrow(as.data.frame(omics[[i]])), 0, noiseGroup)
      # Add variability to each omic compared to group 1
      # transform the matrix into TRUE when > 0 false when 0
      sim_trueFalse <- (omics[[i]] > 0)
      # Multiply the variability vector by a 1/0 to keep the zeros.
      for (c in 1:length(colnames(omics[[i]]))){
        sim_trueFalse[, c] <- as.integer(as.logical(sim_trueFalse[,c]))
        omics[[i]][,c] <- omics[[i]][,c] + (sim_trueFalse[, c] * var_group)
        VARlist[[paste0("Var_", names(omics)[i])]][, c] <- (sim_trueFalse[, c] * var_group)
          
      }
      omics[[i]] <- abs(omics[[i]])
      
    }
  } else if (group == 1){
    # If its the first group, we dont need to add FC, so we multiply by one
    # Instead of a fold change vector
    
    for(i in 1:N_omics){
      message(paste0("Estimating distribution from original data type: ", i))
      param_est <- MOSim::sparsim_estimate_parameter_from_data(raw_data = omics[[i]],
                                                               norm_data = norm_list[[i]],
                                                               conditions = cellTypes)
      param_est_list[[paste0("param_est_", names(omics)[i])]] <- param_est
      
    }
    
    FClist <- list()
    VARlist <- list()
    
    for(i in 1:N_omics){
      # Dont add variability for first group
      VARvec <- rep(1, dim(omics[[i]])[[1]])
      VARlist[[paste0("Var_", names(omics)[i])]] <- VARvec
    }
  }

  N_param_est_list<- length(param_est_list)
  N_cellTypes <- length(cellTypes)
  param_est_list_mod <- list()
  
  for(i in 1:N_param_est_list){
    cell_type_list <- list()
    message(paste0("Creating parameters for omic: ", i))
    
    if (identical(names(omics)[i], "scRNA-seq")){
      FClist[[paste0("FC_est_", names(omics)[i])]] <- as.numeric(stats::na.omit(associationMatrix$Gene_FC))
    } else {
      FClist[[paste0("FC_est_", names(omics)[i])]] <- as.numeric(stats::na.omit(associationMatrix$Peak_FC))
    }
    
    for(j in 1:N_cellTypes){
      message(paste0("Creating parameters for cell type: ", j))
      # Estimate library size  
      if(all_missing){
        libs_param <- param_est_list[[i]][[j]][["lib_size"]]
      } else if (all_specified){
        libs_param <- round(stats::rnorm(n = numberCells[j], mean = mean[j], sd = sd[j]))
      } else if (only_cellnum){
        libs_param <- sample(param_est_list[[i]][[j]][["lib_size"]], 
                             size = numberCells[j], replace = TRUE)
      }
      
      # Since we add the foldchange as a scalar, it's multiplied on top of all
      # celltypes, not only one.
      cond_param <- MOSim::sparsim_create_simulation_parameter(
        intensity = param_est_list[[i]][[j]][["intensity"]] * as.numeric(FClist[[i]]),
        variability = param_est_list[[i]][[j]][["variability"]],
        library_size = libs_param,
        condition_name = param_est_list[[i]][[j]][["name"]],
        feature_names = names(param_est_list[[i]][[j]][["intensity"]]))
      
      cell_type_list[[names(cellTypes)[j]]] <- cond_param
      
    }
    param_est_list_mod[[paste0("param_est_", names(omics)[i])]] <- cell_type_list
  }
  
  return(list(param_list = param_est_list_mod, FClist = FClist,
              VARlist = VARlist, 
              associationMatrix = associationMatrix, dfPeakNames = dfPeakNames, 
              dfGeneNames = dfGeneNames))
}



#' scMOSim
#'
#' Performs multiomic simulation of single cell datasets
#'
#' @param omics named list containing the omic to simulate as names, which can 
#'     be "scRNA-seq" or "scATAC-seq".
#' @param cellTypes list where the i-th element of the list contains the column 
#'     indices for i-th experimental conditions. List must be a named list.
#' @param numberReps OPTIONAL. Number of replicates per group
#' @param numberGroups OPTIONAL. number of different groups
#' @param diffGenes OPTIONAL. If number groups > 1, Percentage DE genes to simulate.
#'    List of vectors (one per group to compare to group 1) where the vector
#'    contains absolute number of genes for Up and Down ex: c(250, 500) or a 
#'    percentage for up, down ex: c(0.2, 0.2). The rest will be NE
#' @param minFC OPTIONAL. Threshold of FC below which are downregulated, by 
#'    default 0.25
#' @param maxFC OPTIONAL. Threshold of FC above which are upregulated, by default 4
#' @param numberCells OPTIONAL. Vector of numbers. The numbers correspond to the number of 
#'     cells the user wants to simulate per each cell type. The length of the 
#'         vector must be the same as length of \code{cellTypes}.
#' @param mean OPTIONAL. Vector of numbers of mean depth per each cell type. Must be specified 
#'     just if \code{numberCells} is specified.The length of the vector must be 
#'         the same as length of \code{cellTypes}.
#' @param sd OPTIONAL. Vector of numbers of standard deviation per each cell type. Must be 
#'     specified just if \code{numberCells} is specified.The length of the vector 
#'         must be the same as length of \code{cellTypes}.
#' @param noiseRep OPTIONAL. Number indicating the desired standard deviation 
#'      between biological replicates.
#' @param noiseGroup OPTIONAL. Number indicating the desired standard deviation
#'      between treatment groups
#' @param regulatorEffect OPTIONAL. To simulate relationship scRNA-scATAC, list 
#'      of vectors (one per group) where the vector contains absolute number of
#'      regulators for Activator and repressor ex: c(150, 200) or a percentage
#'      for Activator and repressor ex: c(0.2, 0.1). The rest will be NE. If not
#'      provided, no table of association between scRNA and scATAC is outputted.
#' @param associationList REQUIRED A 2 columns dataframe reporting peak ids 
#'      related to gene names. If user doesnt have one, load from package
#'      data("associationList")
#' @param feature_no OPTIONAL. If only scRNA-seq to simulate or scRNA and scATAC
#'      but no regulatory constraints, total number of features to be distributed 
#'      between the coexpression clusters.
#' @param clusters OPTIONAL. Number of co-expression patterns the user wants
#'      to simulate
#' @param cluster_size OPTIONAL. It may be inputted by the user. Recommended: 
#'      by default, its the number of features divided by the number of patterns 
#'      to generate.
#' @param TF OPTIONAL default is FALSE, if true, extract TF dataframe
#' @param TFdf OPTIONAL, default is NULL. If an association matrix of TF and 
#'        Target_gene is given the TF expression values are extracted. If no data.frame
#'        is given, using the association of human TF from 
#'        {https://tflink.net/}
#' @return a list of Seurat object, one per each omic.
#' @export
#'
#' @examples
#' omic_list <- sc_omicData(list("scRNA-seq"))
#' cell_types <- list('Treg' = c(1:10),'cDC' = c(11:20),'CD4_TEM' = c(21:30),
#' 'Memory_B' = c(31:40))
#' sim <- scMOSim(omic_list, cell_types, regulatorEffect = list(c(0.1, 0.2)))
#'
scMOSim <- function(omics, cellTypes, numberReps = 1, numberGroups = 1, 
                    diffGenes = NULL, minFC = 0.25, maxFC = 4,
                    numberCells = NULL, mean = NULL, sd = NULL, noiseRep = 0.1 , 
                    noiseGroup = 0.5, regulatorEffect = NULL, associationList = NULL, 
                    feature_no = 8000, clusters = 3, cluster_size = NULL,
                    TF = FALSE, TFdf = NULL){
  
  # Check for mandatory parameters
  if (missing(omics)){
    stop("You must provide the vector of omics to simulate.")
  }
  
  if (missing(cellTypes)) {
    stop("You must provide the correspondence of cells and celltypes")
  }
  
  cellTypesCoex <- cellTypes
  if (is.null(numberCells)){
    numberCells <- lengths(cellTypes)
  }

  ## Check that number of groups and number of differentially expressed
  # probabilities makes sense
  if (numberGroups > 1){
    if (is.null(diffGenes) || length(diffGenes) != (numberGroups - 1)){
      stop(paste0("Number of elements in diffGenes must have a length equal to",
                  " numberGroups -1"))
    }
  }
  
  ## Check that columns of association list are c("Peak_ID", "Gene_ID")
  if (!is.null(associationList) && !identical(colnames(associationList), c("Peak_ID", "Gene_ID"))){
    stop("Column names of the user-inputted association list should be 'Peak_ID' and 'Gene_ID'")
  } else if (is.null(associationList)){
    ## Get the association list loaded in the package
    message("Loading default association list from MOSim package")
    data_env <- new.env(parent = emptyenv())
    data("associationList", envir = data_env, package = "MOSim")
    associationList <- data_env[["associationList"]]
  }
  
  # Check if we have to extract TF
  if (isTRUE(TF)){
    if (is.null(TFdf)){
      # Extract TF from the TF_human vector loaded in the package
      data_env <- new.env(parent = emptyenv())
      data("TF_human", envir = data_env, package = "MOSim")
      TF_human <- data_env[["TF_human"]]
    } else {
      if (!is.data.frame(TFdf) || dim(TFdf)[1] <= 1) {
        stop("Error: variable 'TFdf' is not a dataframe of TF and Target_gene")
      } else{
        TF_human <- TFdf
      }
    }
  }
  
  if (is.null(regulatorEffect) && identical(names(omics[2]), "scATAC-seq")){
    stop("You requested a simulation of scATAC-seq data but did not provide information for variable <regulatorEffect>")
  }
  
  if (numberGroups > 0 && identical(names(omics[2]), "scATAC-seq")){
    if (is.null(regulatorEffect) || length(regulatorEffect) != (numberGroups)){
      stop(paste0("Number of elements in regulatorEffect must have a length equal to",
                  " numberGroups"))
    }
  }
  
  ## Message the experimental dessign
  message(paste0("The experimental design includes: 
                 - ", numberReps, " Biological replicates
                 - ", numberGroups, " Experimental groups
                 - ", paste0(diffGenes), " Differentially expressed genes (Up/Down) per group
                 - ", minFC, " FC below which a gene is downregulated
                 - ", maxFC, " FC above which a gene is upregulated
                 - ", paste0(numberCells), " Number of cells per celltype
                 - ", paste0(regulatorEffect), " Regulators (activator/repressor) per group
                 - ", clusters, " Gene co-expression patterns")[1])
  
  ## format number regulators, if its relative, make absolute
  numfeat <- length(associationList$Gene_ID)
  genereggroup <- list()
  
  if (length(omics) > 1){
    for (i in 1:numberGroups){
      if (regulatorEffect[[i]][1] < 1) {
        # If relative, make absolute numbers
        numActivator <- round(regulatorEffect[[i]][1]*numfeat, digits = 0)
        numRepressor <- round(regulatorEffect[[i]][2]*numfeat, digits = 0)
      } else {
        numActivator <- regulatorEffect[[i]][1]
        numRepressor <- regulatorEffect[[i]][2]
        if (numActivator + numRepressor > numfeat){
          stop(paste0("Number of requested Activators and Repressors is bigger
                    than the possible regulators according to the association
                    dataframe. Activators plus repressors must be below: ",
                      numfeat))
        }
      }
      # Get names of regulators per group
      numActivator <- sample(associationList$Peak_ID, numActivator)
      numRepressor <- sample(setdiff(associationList$Peak_ID, numActivator), numRepressor)
      
      # If group 2 onwards
      if (i > 1){
        ## Get activated, repressed and other diffexp for genes
        if (diffGenes[[i -1]][1] < 1) {
          numup <- round(diffGenes[[i -1]][1]*nrow(omics[[1]]), digits = 0)
          numdown <- round(diffGenes[[i -1]][2]*nrow(omics[[1]]), digits = 0)
        } else {
          numup <- diffGenes[[i -1]][1]
          numdown <- diffGenes[[i -1]][2]
          if (numup + numdown > nrow(omics[[1]])){
            stop(paste0("Number of requested Upregulated and Downregulated genes
                    is bigger than the number of total genes: ", nrow(omics[[1]])))
          }
        }
        
        
        # Get the genes corresponding to the activator regulators and repressors
        genesActivated <- associationList[associationList$Peak_ID %in% numActivator, ]
        genesRepressed <- associationList[associationList$Peak_ID %in% numRepressor, ]
        
        # And add other random genes (not in the association list)
        u <- abs(numup - length(genesActivated[[2]]))
        d <- abs(numdown - length(genesRepressed[[2]]))
        
        availGenes <- setdiff(rownames(omics[[1]]), as.vector(associationList[[2]]))
        genesUp <- sample(availGenes, u)
        availGenes <- setdiff(availGenes, genesUp)
        genesDown <- sample(availGenes, d)
        
        remaining <- setdiff(rownames(omics[[1]]), genesActivated[[2]])
        remaining <- setdiff(remaining, genesRepressed[[2]])
        remaining <- setdiff(remaining, genesUp)
        remaining <- setdiff(remaining, genesDown)
        genereggroup[[paste0("GeneActivated_G", i)]] <- genesActivated
        genereggroup[[paste0("GeneRepressed_G", i)]] <- genesRepressed
        genereggroup[[paste0("GeneExtraUp_G", i)]] <- genesUp
        genereggroup[[paste0("GeneExtraDown_G", i)]] <- genesDown
        genereggroup[[paste0("GeneRemaining_G", i)]] <- remaining
        
        # Here the up and down, have to be as many as possible from the association
        # list, so take rows where its activator and make them up, when its repressor
        # make them down
        if (numup < length(numActivator) || numdown < length(numRepressor)){
          stop("You have asked for many regulators, but there aren't enough 
           differentially expressed genes to be regulated")
        }
        # Get the info for the ATAC
        u <- numup - length(unique(genesActivated[[1]]))
        d <- numdown - length(unique(genesRepressed[[1]]))
        availFeat <- setdiff(rownames(omics[[2]]), unique(associationList[[1]]))
        regAct <- sample(availFeat, u)
        availFeat <- setdiff(availFeat, regAct)
        regRep <- sample(availFeat, d)
        
        remaining <- setdiff(rownames(omics[[2]]), as.vector(genesActivated[[1]]))
        remaining <- setdiff(remaining, as.vector(genesRepressed[[1]]))
        remaining <- setdiff(remaining, regAct)
        remaining <- setdiff(remaining, regRep)
        
        genereggroup[[paste0("FeatExtraUp_G", i)]] <- regAct
        genereggroup[[paste0("FeatExtraDown_G", i)]] <- regRep
        genereggroup[[paste0("FeatRemaining_G", i)]] <- remaining
      } else {
        ## Here say what to put in the genereggroup if we dont have two groups
        activated_genes <- associationList[associationList$Peak_ID %in% numActivator, ]
        repressed_genes <- associationList[associationList$Peak_ID %in% numRepressor, ]
        genereggroup[["FeatExtraUp_G1"]] <- numActivator
        genereggroup[["FeatExtraDown_G1"]] <- numRepressor
        genereggroup[["FeatRemaining_G1"]] <- setdiff(rownames(omics[[2]]), c(numActivator, numRepressor))
        genereggroup[["GeneActivated_G1"]] <- activated_genes
        genereggroup[["GeneRepressed_G1"]] <- repressed_genes
        genereggroup[["GeneExtraUp_G1"]] <- NA
        genereggroup[["GeneExtraDown_G1"]] <- NA
        genereggroup[["GeneRemaining_G1"]] <- setdiff(rownames(omics[[1]]), c(as.character(activated_genes$Gene_ID), as.character(repressed_genes$Gene_ID)))
      }
    }
  } else {
    if (numberGroups > 1){
      ## Here say what to do if we don't have ATAC-seq
      for (i in 2:numberGroups){
        if (diffGenes[[i -1]][1] < 1) {
          numup <- round(diffGenes[[i -1]][1]*nrow(omics[[1]]), digits = 0)
          numdown <- round(diffGenes[[i -1]][2]*nrow(omics[[1]]), digits = 0)
        } else {
          numup <- diffGenes[[i -1]][1]
          numdown <- diffGenes[[i -1]][2]
          if (numup + numdown > nrow(omics[[1]])){
            stop(paste0("Number of requested Upregulated and Downregulated genes
                      is bigger than the number of total genes: ", nrow(omics[[1]])))
          }
        }
        u <- sample(rownames(omics[[1]]), numup)
        availGenes <- setdiff(rownames(omics[[1]]), u)
        d <- sample(availGenes, numdown)
        genereggroup[[paste0("GeneExtraUp_G", i)]] <- u
        genereggroup[[paste0("GeneExtraDown_G", i)]] <- d
        genereggroup[[paste0("GeneRemaining_G", i)]] <- setdiff(availGenes, d)
      }
    }
  }

  
  ### Start working on the data
  
  N_omics <- length(omics)
  
  ## Subset only columns of interest (our celltypes)
  for (om in 1:N_omics){
    omics[[om]] <- omics[[om]][, unname(unlist(cellTypes))]
  }
  # Reorganize our cellTypes variable according to the subset
  createRangeList <- function(vector, nam) {
    rangeList <- list()
    totalValues <- sum(vector)
    start <- 1
    
    for (i in 1:length(vector)) {
      end <- start + vector[i] - 1
      rangeList[[i]] <- seq(start, end)
      start <- end + 1
    }
    names(rangeList) <- nam
    return(rangeList)
  }
  
  cellTypes <- createRangeList(numberCells, names(cellTypes))
  
  # Make the patterns to simulate coexpression
  lpatterns <- make_cluster_patterns(length(cellTypes), clusters = clusters)
  # Get also the indices of cluster patterns that are opposite and will be 
  # important for the regulators
  genereggroup[["opposite_indices"]] <- lpatterns$opposite_indices
  patterns <- lpatterns$patterns
  
  ## Check if there is any group that is asking for more repressors
  # than what is possible according to the number of opposing genes or peaks 
  # we have available.
  
  
  
  
  
  
  # Simulate coexpression
  for (i in 1:N_omics){
    coexpr_results <- MOSim::simulate_coexpression(omics[[i]],
                                            feature_no = feature_no, cellTypes = cellTypesCoex,
                                            patterns = patterns, cluster_size = cluster_size)
    
    # Get the coexpressed matrix out
    omics[[i]] <- as.data.frame(coexpr_results$sim_matrix)
    rownames(omics[[i]]) <- omics[[i]]$feature
    omics[[i]]$feature <- NULL
    # Get the clusters out
    genereggroup[[paste0("Clusters_", names(omics)[[i]])]] <- coexpr_results$sim_clusters
  }

  # Start the lists we will need to include in the output
  seu_groups <- list()
  Association_list <- list()
  FC_used_list <- list()
  VAR_used_list <- list()
  param_list <- list()
  
  for (g in 1:numberGroups){
    seu_replicates <- list()

    message(paste0("Estimating parameters for experimental group ", g))
    
    param_l <- MOSim::sc_param_estimation(omics, cellTypes, diffGenes, minFC, maxFC, 
                                          numberCells, mean, sd, noiseGroup, g, 
                                          genereggroup)
    

    Association_list[[paste0("AssociationMatrix_Group_", g)]] <- param_l$associationMatrix
    FC_used_list[[paste0("FC_Group_", g)]] <- param_l$FClist
    VAR_used_list[[paste0("VAR_Group_", g)]] <- param_l$VARlist
    param_list <- param_l$param_list
    
    for (r in 1:numberReps){
      message(paste0("Simulating parameters for replicate ", r))
      
      sim_list <- list()
      
      
      for(i in 1:N_omics){
        # Simulate the replicate
        sim <- MOSim::sparsim_simulation(dataset_parameter = param_list[[i]])
        sim <- sim[["count_matrix"]]
        # Generate a standard deviation to add to the matrix
        var_rep <- stats::rnorm(nrow(as.data.frame(omics[[i]])), 0, noiseRep)
        # Add the standard deviation of the replicate
        # transform the matrix into TRUE when > 0 false when 0
        sim_trueFalse <- (sim > 0)
        # Multiply the variability vector by a 1/0 to keep the zeros.
        for (e in 1:length(colnames(sim))){
          sim_trueFalse[, e] <- as.integer(as.logical(sim_trueFalse[,e]))
          sim[,e] <- sim[,e] + (sim_trueFalse[,e] * var_rep)
        }
        # Make sure there are no negative numbers
        sim <- abs(sim)
        # Now pass this modified matrix back
        sim_list[[paste0("sim_", names(omics)[i])]] <- sim
        
      }
      
      seu_obj <- list()
      N_sim <- length(sim_list)
      
      for(i in 1:N_sim){
        
        assay_name <- str_split(names(sim_list)[i], "-")[[1]][1]
        assay_name <- sub("sim_sc","",assay_name)
        
        options(Seurat.object.assay.version = "v3")
        seu <- Seurat::CreateAssayObject(counts = sim_list[[i]],
                                          assay = assay_name,
                                          rownames = rownames(sim_list[[i]]), #explicitly specifying rownames
                                          colnames = colnames(sim_list[[i]])) #and colnames for Seurat obj
        seu_obj[[names(sim_list)[i]]] <- seu
        
      }
      seu_replicates[[paste0("Rep_", r)]] <- seu_obj
    }
    
    seu_groups[[paste0("Group_", g)]] <- seu_replicates
    
  } 
  
  # Bring back final celltypes
  seu_groups[["cellTypes"]] <- cellTypes
  seu_groups[["patterns"]] <- patterns
  # Bring back FC vector
  seu_groups[["FC"]] <- FC_used_list
  seu_groups[["AssociationMatrices"]] <- Association_list
  seu_groups[["Variability"]] <- VAR_used_list
  
  
  # If TF, extract their count matrix
  if (isTRUE(TF)){
    # Extract the TFs according to the users list
    for (group_name in names(seu_groups)) {
      if (grepl( "Group", group_name, fixed = TRUE)){
        group <- seu_groups[[group_name]]
        for (replicate_name in names(group)) {
          replicate <- group[[replicate_name]]
          seu_groups[[group_name]][[replicate_name]]$sim_TF <- subset(replicate$`sim_scRNA-seq`, features = TF_human$TF)
        }
      }
    }
    seu_groups[["TFtoGene"]] <- TF_human
  }
  
  return(seu_groups)
  
}



#' scOmicSettings
#'
#' @param sim a simulated object from scMOSim function
#' @param TF OPTIONAL default is FALSE, if true, extract TF association matrix
#'
#' @return list of Association matrices explaining the effects of each
#' regulator to each gene
#' @export
#'
#' @examples
#'
#' cell_types <- list('Treg' = c(1:10),'cDC' = c(11:20),'CD4_TEM' = c(21:30),
#' 'Memory_B' = c(31:40))
#' omicsList <- sc_omicData(list("scRNA-seq"))
#' sim <- scMOSim(omicsList, cell_types)
#' res <- scOmicSettings(sim)
scOmicSettings <- function(sim, TF = FALSE){
  asma <- sim$AssociationMatrices
  FC <- sim$FC
  
  # Add columns for FC
  for (i in c(1:length(asma))){
    # I have NAs in Gene_ID and Peak_ID, so make the vectors accordingly
    asma[[i]]$Gene_FC <- ifelse(!is.na(asma[[i]]$Gene_ID), FC[[i]][[1]], NA)
    if (length(FC) > 1){
      
      asma[[i]]$Peak_FC <- ifelse(!is.na(asma[[i]]$Peak_ID), FC[[i]][[2]], NA)
    }
  }
  
  TFtoGene_Association <- list()
  
  ## Filter TFtoGene so only genes in our clusters are there
  if (isTRUE(TF)){
    TFtoGene <- sim$TFtoGene
    colnames(TFtoGene) <- c("Gene_ID", "Target_ID")
    # Extract only genes with a cluster from association matrix
    clu <- asma$AssociationMatrix_Group_1[asma$AssociationMatrix_Group_1$Gene_cluster >= 1,]
    TFtoGene <- TFtoGene[TFtoGene$Gene_ID %in% clu$Gene_ID,]
    TFtoGene <- TFtoGene[TFtoGene$Target_ID %in% clu$Gene_ID, ]
    
    # Get opposing patterns
    opposite_indices <- MOSim::check_patterns(sim$patterns)
    
    for (group_name in names(sim)){
      if (grepl( "Group", group_name, fixed = TRUE)){
        TFtoGene_matrix <- data.frame(Gene_ID = character(),
                                      Target_ID = character(),
                                      Cluster_ID = character(),
                                      Cluster_Target = character(),
                                      stringsAsFactors = FALSE)
        
        # Get association name
        association_name <- grep(paste0(group_name, "$"), names(asma), value = TRUE)
        mat <- asma[[association_name]]
        
        # Loop through each row of dataframe B
        for (i in 1:nrow(TFtoGene)) {
          # Extract the Gene_ID and Target_ID from dataframe B
          gene_id <- TFtoGene$Gene_ID[i]
          target_id <- TFtoGene$Target_ID[i]
          
          gene_cluster_id <- stats::na.omit(mat$Gene_cluster[mat$Gene_ID == gene_id])
          target_cluster_id <- stats::na.omit(mat$Gene_cluster[mat$Gene_ID == target_id])
          
          # If a matching Cluster_ID is found, add it to dataframe C
          if (gene_cluster_id[1] > 0 && target_cluster_id[1] > 0) {
            # Determine regulator_effect
            if (gene_cluster_id == target_cluster_id) {
              regulator_effect <- "activator"
            } else if (any(sapply(opposite_indices, function(x) all(x %in% c(gene_cluster_id, target_cluster_id))))) {
              regulator_effect <- "repressor"
            } else {
              regulator_effect <- "NE"
            }
            
            TFtoGene_matrix <- rbind(TFtoGene_matrix, 
                                     data.frame(Gene_ID = gene_id,
                                                Target_ID = target_id,
                                                Cluster_ID = gene_cluster_id,
                                                Cluster_Target = target_cluster_id,
                                                regulator_effect = regulator_effect,
                                                stringsAsFactors = FALSE))
            
          }
        }
        TFtoGene_Association[[group_name]] <- TFtoGene_matrix
      }
    }
    asma[["TFtoGene_Association"]] <- TFtoGene_Association
  }
  
  return(asma)
}

#' scOmicResults
#'
#' @param sim a simulated object from scMOSim function
#' @return list of seurat objects with simulated data
#' @export
#' @examples
#'
#' cell_types <- list('Treg' = c(1:10),'cDC' = c(11:20),'CD4_TEM' = c(21:30),
#' 'Memory_B' = c(31:40))
#' omicsList <- sc_omicData(list("scRNA-seq"))
#' sim <- scMOSim(omicsList, cell_types)
#' res <- scOmicResults(sim)
scOmicResults <- function(sim){
  df <- sim[grepl("Group_", names(sim))]
  return(df)
  }