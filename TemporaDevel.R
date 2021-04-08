# Tempora related (Bader lab)


IdentifyVaryingPWsBJ <- function(object, pval_threshold=0.05){
  
  if (class(object)[1] != "Tempora"){
    stop("Not a valid Tempora object")
  }
  if (is.null(object@n.pcs)){
    stop("BuildTrajectory has not been run. See ?Tempora::BuildTrajectory for details")
  }
  if (is.null(object@cluster.pathways)){
    stop("CalculatePWProfiles has not been run. See ?Tempora::CalculatePWProfiles for details")
  }
  gsva_bycluster <- object@cluster.pathways
  
  significant_pathways <- c()
  for (i in 1:object@n.pcs){
    genes_scaled <- scale(object@cluster.pathways.dr$rotation[,i])
    significant_pathways <- c(names(which(genes_scaled[,1] > 1.5 | genes_scaled[,1] < -1.5)), significant_pathways)
  }
  
  pca_pathways <- sub("%.*", "", significant_pathways)
  pca_pathways <- gsub("\\s*\\([^\\)]+\\)","",pca_pathways)
  pca_pathways_cleaned <- gsub("[[:punct:]]", "", pca_pathways)
  themes <- pca_pathways_cleaned
  
  cat("Fitting GAM models...")
  
  # create all possible paths from beginning state to final state (time point) (white list)
  library(tidyr)
  clTime = object@cluster.metadata[,c("Id", "Cluster_time_score")]
  clTime$Timescore = round(clTime$Cluster_time_score)
  clTime$Clusters = clTime$Id
  combList = list()
  addNext = function (x, tidx, combStr) {
    # cat(file = stderr(), combStr)
    if (!tidx %in% x$Timescore) {combList[[length(combList)+1]] <<- combStr; return(combStr)}
    sepStr = "_"
    if(combStr == "") sepStr = ""
    for (pId in x[x$Timescore == tidx, "Clusters"]) {
      addNext(clTime, tidx + 1, sprintf("%s%s%03d", combStr, sepStr, as.numeric(pId)))
    }
  }
  addNext(clTime, 1, "") 
  goodComb = data.frame(combList = combList %>% unlist)
  goodComb$found = 0
  object@trajectory[apply(object@trajectory[,c("from", "to")], 1,FUN = function(x) paste(x,collapse  = "_")) %>%
    lapply(FUN= function(x) any(grepl(x, goodComb))) %>% unlist,]
  
  apply(object@trajectory[,c("from", "to")], 1,FUN = function(x) sprintf("%03d_%03d", x[1], x[2])) %>%
    lapply(FUN= function(x) {
      mIdx = grepl(x, goodComb$combList);
      goodComb[mIdx, "found"] <<- goodComb[mIdx, "found"] +1
      } ) %>% unlist

  fullPathLength = max(clTime$Timescore) - 1 #how many connections are expected
  goodComb = goodComb[which(goodComb$found == fullPathLength),]
  
  if (nrow(goodComb)<1) {
    cat("No trajectory goes from beginnig to end")
    return(object)
  }
  

  p_vals <- gams <- list()
  for (i in 1:length(themes)){ # loop over pathways (pathway themes)
    print(i)
    if(length(grep(themes[i], rownames(gsva_bycluster))) == 0) {
      p_vals[[i]] <- 1
      gams[[i]] <- NA
      next
    }
    if (length(grep(themes[i], rownames(gsva_bycluster))) > 1){
      plot_df <- data.frame(
        cluster=colnames(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), 
        value=colMeans(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ], na.rm=T))
    } else if (length(grep(themes[i], rownames(gsva_bycluster))) == 1){
      plot_df <- data.frame(
          cluster=names(gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ]), 
          value=gsva_bycluster[grep(themes[i], rownames(gsva_bycluster)), ])
      }
    plot_df$time <- object@cluster.metadata$Cluster_time_score
    
    # here we only want to look at the clusters that are connected in the right manner (maybe even the raw data!!!) 
    round(object@cluster.metadata$Cluster_time_score) # position on time axis
    object@meta.data$Clusters # Clusterzugehoerigkeit
    object@meta.data$barcode # cell names
    object@meta.data$Timepoints # position on time axis as factor in correct order
    as.integer(object@meta.data$Timepoints) # 
    object@trajectory$from # has to be from list of parents (as per time axis)
    object@trajectory$to # has to be child (timepoint -1)
    
    
    
    gams[[i]] <- mgcv::gam(value ~ s(time, k=3, bs='cr'), data=plot_df)
    temp_anova <- mgcv::anova.gam(gams[[i]])
    p_vals[[i]] <- temp_anova$s.pv
  }
  
  names(p_vals) <- names(gams) <- themes
  
  pval_threshold = pval_threshold
  p_vals_adj <- p.adjust(unlist(p_vals[which(unlist(p_vals) > 0)]), method = "BH")
  varying_pathways <- p_vals_adj[which(p_vals_adj < pval_threshold)]
  varying_pathways <- varying_pathways[!duplicated(names(varying_pathways))]
  
  if (length(varying_pathways)==0){
    cat("No temporally varying pathways detected. Please try running IdentifyVaryingPWs with a more relaxed p-value cutoff.")
    #eventhough the function was not successful return the object because in the vignette
    # this function call sets the original object to what is returned and if it is null
    # you loose all the processing you have done until now. 
    return(object)
  } else {
    object@varying.pws <- varying_pathways
    object@gams <- gams
    return(object)
  }
}


