# Title: heat scree plot 
# Date: September 14, 2021
# Author: Sam

# NB: code originally from Meaghan

# initialize the heatplot function 
# outputs heatmap of PCs along with scree plot of variance attributed to each pc
heat_scree_plot <- function(Loadings, Importance, Num, Order, Categorical, Continuous){
  
  
  # adjust according to importance of first PC 
  adjust <- 1-Importance[1]
  pca_adjusted <- Importance[2:length(Importance)]/adjust
  pca_df <- data.frame(adjusted_variance = pca_adjusted, 
                       PC = seq(1:length(pca_adjusted)))
  
  ##############
  # scree plot #
  ##############
  
  # plot variance that each adjusted PC accounts for 
  scree <- ggplot(pca_df[which(pca_df$PC<=Num),],aes(PC,adjusted_variance)) + 
    geom_bar(stat = "identity",color="black",fill="grey") +
    theme_bw()+
    theme(axis.text = element_text(size =12),
          axis.title = element_text(size =15),
          plot.margin=unit(c(1,1.5,0.2,2.25),"cm"))+ylab("Variance") +
    scale_x_continuous(breaks = seq(1,Num,1))
  
  
  ##########################################################
  # Heat map of variance in each variable explained by PCs #
  ##########################################################
  
  # correlate metadata (variables) with PCS
  # Run anova of each PC on each meta data variable
  
  aov_PC_meta <- NULL
  # categorical variables
  if(!is.null(Categorical)) {
    # run ANOVA on each PC
    aov_PC_meta <- lapply(1:ncol(Categorical), 
                          function(covar) sapply(1:ncol(Loadings), 
                                                 function(PC) summary(aov(Loadings[,PC]~
                                                                            Categorical[,covar]))[[1]]$"Pr(>F)"[1]))
    
    # set names according to names of categorical variables
    names(aov_PC_meta) <- colnames(Categorical)
    # create matrix from list
    aov_PC_meta <- do.call(rbind, aov_PC_meta)
  }
  
  cor_PC_meta <- NULL
  # continuous variables 
  if(!is.null(Continuous)) {
    # conduct spearman correlation for continuous vars
    cor_PC_meta <- lapply(1:ncol(Continuous),
                          function(covar) sapply(1:ncol(Loadings),
                          function(PC) (cor.test(Loadings[,PC],
                                                 as.numeric(Continuous[,covar]),
                                                 alternative = "two.sided", method="spearman",
                                                 na.action=na.omit, exact=FALSE)$p.value)))
    
    # rename according to names of continuous variables
    names(cor_PC_meta)<-colnames(Continuous)
    # create matrix from list
    cor_PC_meta <- do.call(rbind, cor_PC_meta)
  }
  

  #######################################################
  # Prepare continous and categorical data for heat map #
  #######################################################
  
  # first, remove first pc since all others have been adjusted in relation
  # do sep for whether they are cont, cat, or both var types
  
  allvar_PC_adjust <- NULL
  # both continuous and categorical vars
  if(!is.null(Categorical) & !is.null(Continuous)){
    allvar_PC <- as.data.frame(rbind(aov_PC_meta, cor_PC_meta))
    allvar_PC_adjust <- allvar_PC[,2:ncol(allvar_PC)]
  }
  
  # only categorical vars
  if(!is.null(Categorical) & is.null(Continuous)){
    allvar_PC_adjust <- aov_PC_meta[,2:ncol(aov_PC_meta)]
  }
  
  # only continuous vars
  if(is.null(Categorical) & !is.null(Continuous)){
    allvar_PC_adjust <- cor_PC_meta[,2:ncol(cor_PC_meta)]
  }
  
  
  #######################################
  # Clean all variable PCs for plotting #
  #######################################
  
  # plot number of PCs specified by user; reduces number of columns to be equal to "num"
  plotting_PCs <- allvar_PC_adjust[,1:Num]
  
  # convert plotting PCs to numeric
  plotting_PCs_num <- apply(plotting_PCs,2, as.numeric)
  
  # convert plotting PCs to dataframe
  plotting_PCs_num <- as.data.frame(plotting_PCs_num)
  
  # Rename columns to PC 1, 2, 3, etc...
  colnames(plotting_PCs_num) <- sapply(1:Num, function(x) paste("PC",x, sep=""))
  
  # Create column to name rows accordiong to variable names (metadata)
  plotting_PCs_num$meta <- rownames(plotting_PCs[1:nrow(plotting_PCs),])
  
  # Melt dataframe to long format for plotting
  plotting_PCs_melt <- reshape2::melt(plotting_PCs_num, id.vars="meta")
  
   #Cluster metadata according to order specified by user
  ord <- Order
  plotting_PCs_order <- unique(plotting_PCs_melt$meta)[rev(ord)]
  plotting_PCs_melt$meta <- factor(plotting_PCs_melt$meta, levels = plotting_PCs_order)
  
  # hard code colours for heat map into dataframe according to PC significance
  plotting_PCs_melt$Pvalue<-sapply(1:nrow(plotting_PCs_melt), function(x)
    if(plotting_PCs_melt$value[x]<=0.001){"<=0.001"}else{
      if(plotting_PCs_melt$value[x]<=0.01){"<=0.01"}else{
        if(plotting_PCs_melt$value[x]<=0.05){"<=0.05"}else{">0.05"}}})
  
  
  ######################
  # heat map in ggplot #
  ######################
  
  heat <- ggplot(plotting_PCs_melt, aes(variable, meta, fill = Pvalue)) +
    geom_tile(color = "black",size=0.5) +
    theme_gray(8) + 
    scale_fill_manual(breaks = c("<=0.001", "<=0.01", "<=0.05", ">0.05"),
                      values=c("#084594","#4292c6","#9ecae1","#ffffff")) + 
    theme(axis.text = element_text(size =10, color="black"),
          axis.text.x = element_text(),
          axis.title = element_text(size =15),
          legend.text = element_text(size =14),
          legend.title = element_text(size =12),
          legend.position = "bottom",
          plot.margin=unit(c(0,2.25,1,1),"cm"))+
    xlab("Principal Component")+ylab(NULL)
  
  
  ###################
  # heat scree plot #
  ###################
  
  # use cow plot function to make it look nice
  cowplot::plot_grid(scree, heat, ncol=1)
}

