# Title: Correlation matrix calculations
# Date: August 30, 2021
# Author: Sam Lee


# This script contains the functions to calculate correlations and pvalues need to make a heatmap.
# The corr.mat function calculates correlations. A polychoric calculation is used for factor vs factor (new 2021; as per Lisa Lix).
# tryCatch statements added so variables with few observations do not need to be omitted beforehand.


# load required packages
require(stats) # correlations 
require(polycor) # polychoric correlation for factored vars
require(lsr) # for getting correlation from aov
require(reshape2) # for melting correlation matrix 
require(ggplot2) # for plotting the matrix


############################
# correlation calculations #
############################

# function to examine correlations in data frame
corr.mat = function(df){
  
  print("checking dataframe clases")
  stopifnot(inherits(df, "data.frame")) # must be a data frame
  stopifnot(sapply(df, class) %in% c("integer", "numeric", "factor", "character")) # columns must be one of these classes
  print("calculating correlation matrix")
  
  cor_fun <- function(pos_1, pos_2){
    
    #################### 
    # both are numeric #
    ####################
    
    if(class(df[[pos_1]]) %in% c("integer", "numeric") &&
       class(df[[pos_2]]) %in% c("integer", "numeric")){
      
      print("both are numeric") # print statement for debugging
      print(paste(colnames(df)[[pos_1]], colnames(df)[[pos_2]])) # print statement for debugging
      
      r <- tryCatch({cor(df[[pos_1]], df[[pos_2]],  use = "pairwise.complete.obs")},  # correlation from stats package
                    error = function(e) {NA})
    }
    
    ####################################### 
    # first is numeric,  second is factor #
    #######################################  
    
    if(class(df[[pos_1]]) %in% c("integer", "numeric") &&
       class(df[[pos_2]]) %in% c("factor", "character")){
      
      print("first is numeric, second is factor") # print statements for debugging
      print(paste(colnames(df)[[pos_1]], colnames(df)[[pos_2]])) # print colnames for debugging
      print("checking levels of pos_2 factor") # print statement for debugging
      
      # only use complete cases 
      tmp <- data.frame(df[[pos_1]], df[[pos_2]], row.names = NULL)
      tmp2 <- tmp[complete.cases(tmp),]
      tmp2[,2] <- as.factor(tmp2[,2])
      
      if(is.factor(tmp2[,2]) & (nlevels(droplevels(tmp2[,2]))<2)){
        r <- NA # if less than two levels, cant compare, set to NA 
      } else {
        if(is.factor(tmp2[,2]) & (nlevels(droplevels(tmp2[,2]))>=2)){
          r <- tryCatch({stats::aov(df[[pos_1]] ~ as.factor(df[[pos_2]])) %>% # aov from stats package
              lsr::etaSquared(type=3) %>% `[`(1) %>% sqrt()},  # eta from lsr package
              error = function(e) {NA})
        } else {
          r <- NA
        }
      }
      rm(tmp,tmp2)
    }
    
    ####################################### 
    # first is factor,  second is numeric #
    ####################################### 
    
    if(class(df[[pos_2]]) %in% c("integer", "numeric") &&
       class(df[[pos_1]]) %in% c("factor", "character")){
      
      print("first is factor, second is numeric") # print statements for debugging
      print(paste(colnames(df)[[pos_1]], colnames(df)[[pos_2]])) # print colnames for debugging
      print("checking levels of first factor variable") # print statement for debugging
      
      # use complete cases in calculations 
      tmp <- data.frame(df[[pos_1]], df[[pos_2]], row.names = NULL)
      tmp2 <- tmp[complete.cases(tmp),]
      tmp2[,1] <- as.factor(tmp2[,1])
      
      if(is.factor(tmp2[,1]) & (nlevels(droplevels(tmp2[,1]))<2)){
        r <- NA # if less than two levels, cant compare, set to NA 
      } else {
        if(is.factor(tmp2[,1]) & (nlevels(droplevels(tmp2[,1]))>=2)){
          r <-  tryCatch({ stats::aov(df[[pos_2]] ~ as.factor(df[[pos_1]])) %>% # aov from stats package
              lsr::etaSquared(type=3) %>% `[`(1) %>% sqrt()},  # eta from lsr package
              error = function(e) {NA})
        } else {
          r <- NA
        }
      }
      rm(tmp,tmp2)
    }
    
    ###########################  
    #both are factor/character#
    ###########################
    
    if(class(df[[pos_1]]) %in% c("factor", "character") &&
       class(df[[pos_2]]) %in% c("factor", "character")){
      
      print("both are characters") # print statement for debugging
      print(paste(colnames(df)[[pos_1]], colnames(df)[[pos_2]])) # print colnames for debugging 
      
      print("checking levels of of both positions") # print statement for debugging
      
      # only use complete cases in calculations 
      tmp <- data.frame(df[[pos_1]], df[[pos_2]], row.names = NULL)
      tmp2 <- tmp[complete.cases(tmp),]
      tmp2[,1] <- as.factor(tmp2[,1])
      tmp2[,2] <- as.factor(tmp2[,2])
      
      
      if(colnames(df)[[pos_1]] == colnames(df)[[pos_2]]){
        r <- 1 # when comparing var to itself set pvalue to 1
      } else {
        if(is.factor(tmp2[,1]) && (nlevels(droplevels(tmp2[,1]))<2)){
          r <- NA # if fewer than 2 levels set to NA
        } else {
          if(is.factor(tmp2[,2]) && (nlevels(droplevels(tmp2[,2]))<2)){
            r <- NA # if fewer than 2 levels set to NA
          } else {
           r <- tryCatch({polychor(tmp2[,2], tmp2[,1])},  # eta from lsr package
                          error = function(e) {NA})
          }
        }
      }
      rm(tmp,tmp2)
    }
    return(r)
  } 
  
  # create correlation matrix
  cor_fun <- Vectorize(cor_fun)
  corrmat <- outer(1:ncol(df), 1:ncol(df), function(x, y) cor_fun(x, y))
  
  # set row and col names of matrix to colnames of df
  rownames(corrmat) <- colnames(df)
  colnames(corrmat) <- colnames(df)
  return(corrmat)
}


########################
# p-value calculations #
########################

# function to examine correlation PVALUES in data frame
pval.mat <- function(df){
  
  #check that classes are appropriate
  print("checking dataframe classes...")
  stopifnot(inherits(df, "data.frame")) # must be data frame
  stopifnot(sapply(df, class) %in% c("integer", "numeric", "factor", "character")) # columns must be one of these classes
  print("calculating correlation matrix p-values...")
  
  ##################################
  #function for calculating pvalues#
  ##################################
  
  pval_fun <- function(pos_1, pos_2){
    
    ####################
    # both are numeric #
    ####################
    if(class(df[[pos_1]]) %in% c("integer", "numeric") &&
       class(df[[pos_2]]) %in% c("integer", "numeric")){
      
      print("both are numeric") # print statement for debugging
      
      # case where comparing var against itself
      if(colnames(df)[[pos_1]] == colnames(df)[[pos_2]]){
        p <- 1
      } else {
        print(paste(colnames(df)[[pos_1]], colnames(df)[[pos_2]])) # print col names for debugging
        p <- tryCatch({unname((summary(lm(df[[pos_1]] ~ df[[pos_2]]))$coefficients[,4])[2])}, # lm function for pvalue 
                      error = function(e) {NA})
      }
    }
    
    ##########################################
    # first is integer and second is numeric #
    ##########################################
    
    if(class(df[[pos_1]]) %in% c("integer", "numeric") &&
       class(df[[pos_2]]) %in% c("factor", "character")){
      
      print("first is numeric, second is factor") # print statement for debugging
      print(paste(colnames(df)[[pos_1]], colnames(df)[[pos_2]])) # print names for debugging 
      print("Checking levels of second factor variable...") # print statement for debugging
      
      # use only complete cases
      tmp <- data.frame(df[[pos_1]], df[[pos_2]], row.names = NULL)
      tmp2 <- tmp[complete.cases(tmp),]
      tmp2[,2] <- as.factor(tmp2[,2])
      
      if(is.factor(tmp2[,2]) & (nlevels(droplevels(tmp2[,2]))<2)){
        p <- NA # cases where dropped values result in a factor with one or fewer levels and set to NA
      } else {
        if(nlevels(droplevels(tmp2[,2]))>2){
          # ANOVA for factor with 3 or more levels
          lm_p <- lm(df[[pos_1]] ~ as.factor(df[[pos_2]]))
          p <-  tryCatch({anova(lm_p)$`Pr(>F)`[1]}, 
                         error = function(e) {NA})
        } else{
          # t-test for 2 factors
          p <- tryCatch({t.test(df[[pos_1]] ~ as.factor(df[[pos_2]]), alternative = "two.sided")$p.value}, 
                        error = function(e) {NA})
        }
      }
      rm(tmp,tmp2)
    }
    
    ##########################################
    # second is integer and first is numeric #
    ##########################################
    
    if(class(df[[pos_2]]) %in% c("integer", "numeric") &&
       class(df[[pos_1]]) %in% c("factor", "character")){
      
      print("first is factor, second is numeric") # print statement for debugging
      print(paste(colnames(df)[[pos_1]], colnames(df)[[pos_2]])) # print names for debugging
      print("checking levels of first factor variable...") # print statemetn for debugging
      
      # use only complete cases in calculations 
      tmp <- data.frame(df[[pos_1]], df[[pos_2]], row.names = NULL)
      tmp2 <- tmp[complete.cases(tmp),]
      tmp2[,1] <- as.factor(tmp2[,1])
      
      
      if(is.factor(tmp2[,1]) && (nlevels(droplevels(tmp2[,1]))<2)){
        p <- NA # cases where dropped values result in a factor with one or fewer levels and set to NA
      } else {
        if(nlevels(droplevels(tmp2[,1]))>2){
          # ANOVA for factor with 3 or more levels
          lm_p <- lm(df[[pos_2]] ~ as.factor(df[[pos_1]]))
          p <-  tryCatch({anova(lm_p)$`Pr(>F)`[1]}, 
                         error = function(e) {NA})
        } else{
          # t-test for 2 factors
          p <-  tryCatch({t.test(df[[pos_2]] ~ as.factor(df[[pos_1]]), alternative = "two.sided")$p.value}, 
                         error = function(e) {NA})
        }
      }
      rm(tmp,tmp2)
    }
    
    #############################
    # both are factor/character #
    #############################
    
    if(class(df[[pos_1]]) %in% c("factor", "character") &&
       class(df[[pos_2]]) %in% c("factor", "character")){
      
      print("both are characters") # print statement for debugging
      print(paste(colnames(df)[[pos_1]], colnames(df)[[pos_2]])) # print colnames for debugging 
      print("checking levels of both positions") # print statement for debugging
      
      # use complete cases for calculations 
      tmp <- data.frame(df[[pos_1]], df[[pos_2]], row.names = NULL)
      tmp2 <- tmp[complete.cases(tmp),]
      tmp2[,1] <- as.factor(tmp2[,1])
      tmp2[,2] <- as.factor(tmp2[,2])
      
      if(colnames(df)[[pos_1]] == colnames(df)[[pos_2]]){
        p <- 1 # when comparing var to itself set pvalue to 1
      } else {
        if(is.factor(tmp2[,1]) && (nlevels(droplevels(tmp2[,1]))<2)){
          p <- NA # cases where dropped values result in a factor with one or fewer levels and set to NA
        } else {
          if(is.factor(tmp2[,2]) && (nlevels(droplevels(tmp2[,2]))<2)){
            p <- NA # cases where dropped values result in a factor with one or fewer levels and set to NA
          } else {
            #get p value
            fish_table <- table(df[[pos_1]], df[[pos_2]])
            p <- tryCatch({fisher.test(fish_table, simulate.p.value = T)$p.value}, 
                          error = function(e) {NA})
          }
        }
      }
      rm(tmp,tmp2)
    }
    return(p)
  }
  
  # create matrix 
  pval_fun <- Vectorize(pval_fun)
  pvalmat <- outer(1:ncol(df), 1:ncol(df), function(x, y) pval_fun(x, y))
  
  # set row names and colnames to colnames of df
  rownames(pvalmat) <- colnames(df)
  colnames(pvalmat) <- colnames(df)
  
  # return pvalue matrix 
  return(pvalmat)
}


####################################################################
# functions to get lower and upper triangles of correlation matrix #
####################################################################

# upper and lower funs modified by kurt to remove diagonal? 

# lower triangle 
get.lower.tri <-  function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[upper.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}


# upper triangle 
get.upper.tri <- function(cormat, show.diag = FALSE) {
  if (is.null(cormat)) {
    return(cormat)
  }
  cormat[lower.tri(cormat)] <- NA
  if (!show.diag) {
    diag(cormat) <- NA
  }
  return(cormat)
}


##############################
# cluster correlation matrix #
##############################

# hc.order correlation matrix
corr.hclust <- function(cormat, hc.method = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")) {
  hc.method <- match.arg(hc.method)
  dd <- stats::as.dist((1 - cormat) / 2)
  hc <- stats::hclust(dd, method = hc.method)
  hc$order
}

###############################
# correlation matrix plotting #
###############################

corr.plot <- function(corr,
                      pval = NULL,
                      method = c("square", "circle"), 
                      type = c( "lower", "full", "upper"), 
                      show.legend = TRUE, 
                      show.diag = FALSE, 
                      colors = c("blue", "white", "red"), 
                      outline.colour = "gray", 
                      hclust = FALSE,
                      hclust.method = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"), 
                      label.sig = FALSE, 
                      label.colour = "black",
                      label.size = 4, 
                      sig.level = 0.05, 
                      label.insig = FALSE, 
                      shape = 4, 
                      shape.colour = "black", 
                      shape.size = 5,
                      axis.text.size = 12, 
                      axis.text.colour = "black", 
                      no.na = TRUE,
                      x.axis.text.angle = 45, 
                      digits = 2) {
  
  # defaults to value at first index
  type <- match.arg(type)
  method <- match.arg(method)

  
  # must be a matrix or dataframe 
  if (!is.matrix(corr) & !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  
  # if a data frame, convert to matrix 
  if (is.data.frame(corr)){
    corr <- as.matrix(corr) 
  }
  
  # round correlation coefficients based on user input
  corr <- base::round(x = corr, digits = digits)
  
  # if the user has specified clustering should occur 
  if (hclust) {
    hc.method <- match.arg(hclust.method)
    ord <- corr.hclust(corr, hc.method)
    corr <- corr[ord, ord]
    if (!is.null(pval)) {
      pval <- pval[ord, ord]
      pval <- base::round(x = pval, digits = digits)
    }
  }
  
  # if the user has specified upper or lower triangle, get these
  if (type == "lower") {
    corr <- get.lower.tri(corr, show.diag)
    pval <- get.lower.tri(pval, show.diag)
  }
  else if (type == "upper") {
    corr <- get.upper.tri(corr, show.diag)
    pval <- get.upper.tri(pval, show.diag)
  }
  
  # melt matrix 
  corr <- reshape2::melt(corr, na.rm = F) # melt the correlation matrix 
  colnames(corr) <- c("Var1", "Var2", "value") # set column names
  corr$pvalue <- rep(NA, nrow(corr)) # create columns for p values
  corr$signif <- rep(NA, nrow(corr)) # create column for specifying significance
  
  # for pvalue matrix, melt and get significance
  if (!is.null(pval)) {
    pval <- reshape2::melt(pval, na.rm = F)
    corr$coef <- corr$value
    corr$pvalue <- pval$value
    corr$signif <- as.numeric(pval$value <= sig.level)
  }
  
  
  # remove rows with missing data
  if(no.na){
    # missing either pvalue or correlation coefficient 
    nas <- unique(c(which(is.na(corr$value)), which(is.na(pval$value))))
    pval <- pval[-nas, ]
    corr <- corr[-nas, ]
  }
  
  # create the base plot using ggplot
  p <- ggplot(data = corr,
              mapping = aes_string(x = "Var1",
                                   y = "Var2",
                                   fill = "value"))
  
  # set whether squares or circles in matrix
  if (method == "square") {
    p <- p + geom_tile(color = outline.colour)
  } 
  else if (method == "circle") {
    p <- p + geom_point(color = outline.colour, 
                        shape = 21,
                        aes_string(size = abs(corr$value)*10)) +
      scale_size(range = c(4,10)) +
      guides(size = FALSE)
  }
  
  
  # set colour scale based on user input
  p <- p + scale_fill_gradient2(low = colors[1],
                                high = colors[3],
                                mid = colors[2],
                                midpoint = 0,
                                limit = c(-1, 1),
                                na.value="white", 
                                space = "Lab",
                                name = "Correlation coefficient")
  
  
  # label squares/circles with pvalue
  label <- round(x = corr[, "signif"], digits = digits)
  if (!is.null(pval) & !label.insig) {
    ns <- corr$pvalue > sig.level
    }
  else if (sum(ns) > 0) {
      label[ns] <- " "
  }
  
  # if significance labels selected 
  if (label.sig) {
    p <- p + geom_point(data = corr %>% subset(signif==1),
                        mapping = aes_string(x = "Var1", y = "Var2"),
                        shape = shape, size = shape.size, color = shape.colour)
  }
  
  # if pval not null 
  if (!is.null(pval) & label.insig) {
    p <- p + geom_point(data = pval,
                        mapping = aes_string(x = "Var1", y = "Var2"),
                        shape = shape, size = shape.size, color = shape.colour)
  }
  
  # final aesthetic changes
  p <- p + 
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(angle = x.axis.text.angle,
                                     colour = axis.text.colour,
                                     size = axis.text.size,
                                     vjust = 1,
                                     hjust = 1),
          axis.text.y = element_text(size = axis.text.size,
                                     colour = axis.text.colour))
  
  # output plot
  p
  
}