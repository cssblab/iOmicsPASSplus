#######################################
#######################################
### Software:      iOmicsPASSv2plus ###
### Developer:     Hiromi K.W.L    ####
### Version:       R-Beta Prog.    ####
### Last modified: 27.02.2020      ####
#######################################
#######################################

#' Network inference module
#'
#' This function helps to create an inferred partial correlation network via two approaches:(1) Supervised and (2) Hybrid.
#' The supervised approach uses glasso to estimate a sparse inversed covariance matrix completely from the data alone.
#' The hybrid approach combines a prior network with an estimated network from the supervised method to create a fused network.
#'
#'
#' @param inputDat a list object containing up to three matrices (X, Y and Z) with features as the row.names and sample IDs on the columns.
#' @param detectOutliers whether to carry out outlier filtering using PCA. Outliers will be removed if outside of 4 SDs from the median of the first two PCs.
#' @param option option=1 for supervised method and option=2 for hybrid method.
#' @param NetworkFile prior network file. File should indicate which data type the feature comes from (i.e. X or Y) using "NodeA_DT" and "NodeB_DT"
#' in case of same gene name. Required if option=2 (hybrid).
#' @param cutoff_fusedmat cut-off value between 0 to 1 used to convert the fused matrix into an adjacency matrix.
#' @param standardization whether to carry out standardization for each data in inputDat. Recommended if integrating multiple datasets before computing cross-covariance matrix.
#' @param log.transform whether to log-transform each matrix in inputDat (default=TRUE).
#' @param toSkipPCA whether to carry out principal component analysis, not recommended if there are too many missing entries in data (default = FALSE).
#' @param useCorrelation whether to compute correlation instead of covariance matrix (default = FALSE). If correlation is selected, standardization of variables is not required.
#' @param Plot.PCA whether to output the plots of PCA of the input data (default=TRUE)
#' @param Plot.modelSelect whether to output the plots to help pick a regularization parameter based on AIC, BIC, eBIC and CV (default=TRUE).
#' @param Plot.Covariance whether to output the heatmap illustrating the cross-covariance matrix. Not recommended if the dimension of the matrix exceeds 2000.
#' @param lambda.vec a vector of regularization parameter to use in the soft-thresholding in glasso.
#' @param numLambda number of values in the regularization parameter.
#' @param numFolds number of folds in the CV model selection
#' @param criterion c("AIC","BIC","eBIC","CV"). Model selection criteria to use to pick the parameter that yield the best penalized log-likelihood.
#' @param optLambda value of the regularization parameter to use to estimate the sparse inverse covariance.
#' @param Calibration whether to carry out the model selection step by fitting glasso on a grid of lambda values. To be turned off after picking an optimal lambda.
#' @param tag a string that will be appended to the end of the output files like a tag.
#' @param verbose whether to output steps and update status of running the function.
#' @return Multiple output files with inferred network and estimated partial correlation matrix.
#' \itemize{
#'   \item PCAplot.pdf - A panel of PCA plots of individual -omics data and combined data highlighting outliers.
#'   \item Boxplots_outliers.pdf - A panel of boxplots of each -omics data across the samples, highlighting outliers
#'   \item Heatmap_CrossCovarianceMatrix.pdf - A heatmap of the calculated cross-covariance matrix of the combined data.
#'   \item CalibrationPlots_glasso.pdf - A multi-panel plot of four model selection criteria against a grid of lambda values to help user pick the optimal value.
#'   \item Plots_glasso.png - a multipanel plot illustrating the derivation of the precision matrix from the cross-covariance matrix.
#'   \item Plots_Hybridmethod.png - a multipanel plot illustrating the derivation of the precision matrix from supervised approach and combining the prior information to form the fused matrix in the hybrid approach.
#'   \item glasso_estimated_icov.txt - a text file containing the estimated sparse inverse of the covariance matrix (also known as precision matrix).
#'   \item PartialCorrelation_icov.txt - a text file containing the corresponding partial correlation calculated from the estimated inverse of covariance.
#'   \item Combined_data.txt - A dataframe with the various -omics data concatenated together after standardization.
#'   \item Estimated_Network_glasso.txt - The corresponding network file from the estimated precision matrix from the supervised approach.
#'   \item Fusednetwork_hybridmethod.txt - The corresponding network file from the estimated fused matrix from the hybrid approach.
#' }
#' @examples
#'
#' data(Tulip_Protein)
#' data(Tulip_microRNA)
#'
#' row.names(Tulip_Protein) = Tulip_Protein$Protein
#' row.names(Tulip_microRNA) = Tulip_microRNA$miRNA
#'
#' Tulip_Protein = Tulip_Protein[,-1]
#' Tulip_microRNA = Tulip_microRNA[,-1]
#'
#' inputDat=list(Tulip_Protein, Tulip_microRNA)
#' names(inputDat) = c("Protein","microRNA")
#'
#' #######################
#' # supervised approach #
#' #######################
#'
#' ## using covariance matrix to estimate network ##
#' NetDeconvolute(inputDat, option=1,log.transform=TRUE, tag="supervised",criterion="eBIC",
#' Calibration=TRUE, verbose=T)
#'
#' ###using correlation matrix to estimate network ###
#' NetDeconvolute(inputDat, option=1,log.transform=TRUE, tag= "correlation",criterion="eBIC",
#' Calibration=TRUE, useCorrelation=T, standardization=F, verbose=T)
#'
#' #' ### Note that after standardization, covariance and correlation is the same ###
#' ### Correlation can be used if data cannot be standardized, for example when looking at changes (post-pre) ###
#'
#' ## continuing if using covariance since we are working with expression data,
#' ## we refine the lambda vector to zoom into a specific window ##
#' lambda_new=exp(seq(log(0.5),log(0.01), length=30))
#'
#' NetDeconvolute(inputDat, option=1,log.transform=TRUE, tag= "supervised2",criterion="eBIC",
#' Calibration=TRUE, lambda.vec=lambda_new, verbose=T)
#'
#' # pick a threshold using the calibration plot with highest CVscore or lowest AIC/BIC/eBIC.
#' # Then turn off calibration to FALSE to re-running the cross-validation for glasso.
#' NetDeconvolute(inputDat, option=1,log.transform=TRUE,tag= "supervised",criterion="eBIC",
#' Calibration=FALSE, optLambda=0.38,verbose=TRUE)
#'
#'
#' ###################
#' # hybrid approach #
#' ###################
#' data(TargetScan_network)
#' data(PPI_network)
#' TargetScan_network$NodeA_DT = "X"
#' TargetScan_network$NodeB_DT = "Y"
#' PPI_network$NodeA_DT = "X"
#' PPI_network$NodeB_DT = "X"
#'
#' colnames(TargetScan_network) = c("NodeA","NodeB","Dir","NodeA_DT","NodeB_DT")
#' colnames(PPI_network) = c("NodeA","NodeB","Dir","NodeA_DT","NodeB_DT")
#' PriorNet = rbind(TargetScan_network,PPI_network)
#'
#' NetDeconvolute(inputDat, option=2, NetworkFile=PriorNet, tag="hybrid",criterion="eBIC",
#' log.transform=TRUE,Calibration=FALSE, optLambda=0.382, verbose=TRUE)
#'
#'
#' @Note The higher the total dimensionality of the data, the longer it will take to run the
#' program. It's suggested to turn Plot.Covariance = FALSE when dimensionality is larger than 2000.
#' @export
#' @export
#' @import huge gplots gridGraphics grid gridExtra pheatmap
#' @importFrom corpcor make.positive.definite
#' @importFrom nnet multinom
#' @importFrom pracma nearest_spd eig
#' @importFrom RColorBrewer brewer.pal
NetDeconvolute <- function(inputDat ,detectOutliers=TRUE, option, Calibration=TRUE, tag= NULL,log.transform =TRUE,cutoff_fusedmat=0.5,
                           standardization=TRUE, NetworkFile= NULL,  toSkipPCA = FALSE, Plot.PCA=TRUE , Plot.modelSelect = TRUE, Plot.Covariance = TRUE,useCorrelation=FALSE,
                           lambda.vec = NULL, numLambda=30, numFolds=5,criterion=c("AIC","BIC","eBIC","CV"),optLambda=NULL, verbose = F){


  if (!file.exists("Output")){
    dir.create("Output")
  }
   ## data cleaning ##
  ntype = length(inputDat)

  if(ntype>3) {
    stop("Only at most 3 types of data is allowed to create network. Please try again.")
  }

  str = names(inputDat)
  commonS = colnames(inputDat[[1]])
  if(ntype>1){
    for(i in 2:ntype) commonS = intersect(commonS,colnames(inputDat[[i]]))
  }


  if(length(commonS)==0){
    stop("No common sample labels in datasets. Please try again.")
  }
  if(length(commonS)<10){
    stop("Number of common sample labels in datasets is less than 10. Recommend to have a sample size of at least 20. Please try again!")
  }
  if(is.null(optLambda) & is.null(criterion)){
    stop("Please select a model selection criterion (AIC, BIC, eBIC or CV) for picking an optimal lambda.")
  }
  if(is.null(optLambda) & !Calibration){
    stop("Please specify a lambda value if calibration is set to FALSE.")
  }

  ndim = sum(unlist(lapply(inputDat,function(x) nrow(x))))
  if(ndim>2000 & (Plot.Covariance == TRUE)){
    stop("Dimension of data is more than 2000. Program may take a while to run...Please be patient! Thank you!\n
    Suggest to set Plot.Covariance = FALSE as it's no longer practical to plot a heatmap with such high dimensionality and to prevent R from crashing.")
  }
  for(i in 1:ntype){
    tmpD = inputDat[[i]]
    tmpDnew = tmpD[,match(commonS, colnames(tmpD))]
    if((any(tmpDnew==0) & log.transform==TRUE)){
      stop("Unable to take log-transform as there are zero entries in the input data. Please check again!\n")
    }
    inputDat[[i]] = tmpDnew
  }
  datY = NULL
  datZ = NULL
  datX =  inputDat[[1]]
  if(ntype>1) datY =  inputDat[[2]]
  if(ntype==3) datZ =  inputDat[[3]]

  if(log.transform){
    for(i in 1:ntype){
      inputDat[[i]] = log(inputDat[[i]])
    }
  }

  Zdat = inputDat
  if(standardization){
    for(i in 1:ntype){
      tmpDD = inputDat[[i]]
      mnX = apply(tmpDD,1,function(x) mean(x,na.rm=T))
      sdX = apply(tmpDD,1, function(x) sd(x,na.rm=T))
      cenX = sweep(tmpDD,1,mnX)
      #mad = apply(tmpDD,1,median)
      Zdat[[i]] = cenX/sdX
    }
  }

  ZdatY = NULL
  ZdatZ = NULL
  ZdatX = Zdat[[1]]
  if(ntype>1) ZdatY = Zdat[[2]]
  if(ntype==3) ZdatZ = Zdat[[3]]

  vec_n = c(nrow(ZdatX))
  listx <- paste0("X__", 1:nrow(ZdatX))
  KEYx = data.frame(listx, row.names(ZdatX))
  colnames(KEYx) = c("Label", "Name")
  KEYx$NewName = paste0(KEYx$Name, "_dataX")
  KEYx$DataType = "X"
  row.names(ZdatX) = listx
  KEY = KEYx
  DAT = ZdatX
  List = c(listx)

  if(ntype>1){
    vec_n = c(vec_n,nrow(ZdatY))
    listy <- paste0("Y__", 1:nrow(ZdatY))
    KEYy = data.frame(listy, row.names(ZdatY))
    colnames(KEYy) = c("Label", "Name")
    KEYy$NewName = paste0(KEYy$Name, "_dataY")
    KEYy$DataType = "Y"
    KEY = rbind(KEY, KEYy)
    row.names(ZdatY) = listy
    DAT = rbind(DAT, ZdatY)
    List = c(List, listy)
  }
  if(ntype==3){
    vec_n = c(vec_n,nrow(ZdatZ))
    listz <- paste0("Z__", 1:nrow(ZdatZ))
    KEYz = data.frame(listz, row.names(ZdatZ))
    colnames(KEYz) = c("Label", "Name")
    KEYz$NewName = paste0(KEYz$Name, "_dataZ")
    KEYz$DataType = "Z"
    KEY = rbind(KEY, KEYz)
    row.names(ZdatZ) = listz
    DAT = rbind(DAT, ZdatZ)
    List = c(List, listz)
  }


  CC = c("royalblue","cyan1","slateblue1")
  COL = NULL
  for(i in 1:ntype){
    COL = c(COL, rep(CC[i],vec_n[i]))
  }

  limits = c(min(DAT, na.rm=T), max(DAT, na.rm=T))

  if(any(is.na(DAT))){
    cid_rm = apply(DAT,1,function(x) sum(is.na(x)))
    if(any(cid_rm> 0.5*ncol(DAT))) stop("There are features with more than 50% missing measurements across all the samples. Please impute or reduce the proportion of missingness before running NetDeconvolute()!\n")
    cid_na = apply(DAT,1,function(x) any(is.na(x)))
    if(sum(cid_na) >0.3*nrow(DAT)){
      cat("There are >30% of the features with missing measurements, we will skip PCA and not perform outlier removal.\n")
      toSkipPCA=TRUE
      detectOutliers=TRUE
    }
    if(sum(cid_na) <=0.3*nrow(DAT)) DAT.pca = DAT[!cid_na,]
  }

  if(all(!is.na(DAT))) DAT.pca = DAT

  if(!toSkipPCA){
    tmpPCA = prcomp(t(DAT.pca), center = T)
    xx = tmpPCA$x
    dev = (tmpPCA$sdev/sum(tmpPCA$sdev))*100

    if(ntype>1){
      if(any(is.na(ZdatX))){
        cid_rm = apply(ZdatX,2,function(x) any(is.na(x)))
        ZdatX.pca = ZdatX[,!cid_rm]
      }
      if(all(!is.na(ZdatX))) ZdatX.pca = ZdatX
      tmpPCA2 = prcomp(t(ZdatX.pca), center = T)
      xx2 = tmpPCA2$x
      dev2 = (tmpPCA2$sdev/sum(tmpPCA2$sdev))*100

      if(any(is.na(ZdatY))){
        cid_rm = apply(ZdatY,2,function(x) any(is.na(x)))
        ZdatY.pca = ZdatY[,!cid_rm]
      }
      if(all(!is.na(ZdatY))) ZdatY.pca = ZdatY
      tmpPCA3 = prcomp(t(ZdatY.pca), center = T)
      xx3 = tmpPCA3$x
      dev3 = (tmpPCA3$sdev/sum(tmpPCA3$sdev))*100

      if(ntype==3){
        if(any(is.na(ZdatZ))){
          cid_rm = apply(ZdatZ,2,function(x) any(is.na(x)))
          ZdatZ.pca = ZdatZ[,!cid_rm]
        }
        if(all(!is.na(ZdatZ))) ZdatZ.pca = ZdatZ
        tmpPCA4 = prcomp(t(ZdatZ.pca), center = T)
        xx4 = tmpPCA4$x
        dev4 = (tmpPCA4$sdev/sum(tmpPCA4$sdev))*100
      }
      ### identifying outliers ###

      if(detectOutliers){
        cat("Carrying out outlier detection step...\n")
        bound1 = 4*sd(xx[,1])
        bound2 = 4*sd(xx[,2])
        mm1 = median(xx[,1])
        mm2 = median(xx[,2])
        cid_out = which((xx[,1] > (mm1+bound1))|(xx[,1]< (mm1 - bound1))|(xx[,2] > (mm2+bound2))|(xx[,2]< (mm2 - bound2)))
        outlier = colnames(DAT)[cid_out]

        ht = 5
        if(ntype>1) ht = ht*(ntype)

        file <- "Output/Boxplot_outliers"
        if(!is.null(tag)) filenew = paste0(file,"_",tag,".pdf")
        else filenew <- "Output/Boxplot_outliers.pdf"

        pdf(filenew,height=ht,width=12,useDingbats = F)
        if(ntype>1) par(mfrow=c(ntype,1),las=2)
        else par(las=2)
        cc = rep("grey70", ncol(DAT))
        if(length(cid_out)>0) cc[cid_out]= "firebrick1"

        boxplot(ZdatX,col = cc, cex=0.7,outline=F, main=paste0(str[1]," data"),xaxt="n")
        if(length(cid_out)>0)axis(1, at=c(1:ncol(ZdatX))[-cid_out],labels=colnames(DAT)[-cid_out], cex.axis=0.7, las=2)
        if(length(cid_out)>0)axis(1, at=cid_out,labels=outlier, cex.axis=0.7, las=2, col.axis="red")
        if(ntype>1){
          boxplot(ZdatY,col = cc, cex=0.7,outline=F, main=paste0(str[2]," data"),xaxt="n")
          if(length(cid_out)>0)axis(1, at=c(1:ncol(ZdatY))[-cid_out],labels=colnames(DAT)[-cid_out], cex.axis=0.7, las=2)
          if(length(cid_out)>0)axis(1, at=cid_out,labels=outlier, cex.axis=0.7, las=2, col.axis="red")
          if(ntype==3){
            boxplot(ZdatZ,col = cc, cex=0.7,outline=F, main=paste0(str[3]," data"),xaxt="n")
            if(length(cid_out)>0)axis(1, at=c(1:ncol(ZdatZ))[-cid_out],labels=colnames(DAT)[-cid_out], cex.axis=0.7, las=2)
            if(length(cid_out)>0)axis(1, at=cid_out,labels=outlier, cex.axis=0.7, las=2, col.axis="red")
          }
        }
        dev.off()
        if(verbose){
          if(length(outlier)==0) cat("No outliers were identified in this step.\n\nProceeding on...\n")
          if(length(outlier)==1) cat(paste0("1 outlier (", outlier  ,") was removed.\n\nProceeding on...\n"))
          if(length(outlier)>1) cat(paste0(length(outlier)," outliers (",paste0(outlier,collapse=","),") were removed.\n\nProceeding on...\n"))
        }
        ## remove outliers identified ##
        if(length(cid_out)>0){
          RM = match(outlier, colnames(DAT))
          DAT = DAT[,-RM]
          ZdatX = ZdatX[,-RM]
          if(ntype>1)ZdatY = ZdatY[,-RM]
          if(ntype==3)ZdatZ = ZdatZ[,-RM]
        }
      }
      if(!detectOutliers) cat("Skipping outlier detection step...\n")

      if(Plot.PCA){
        len = 5.5
        new.main = "all"
        if(ntype>1) len = len*(ntype+1)
        if(ntype==1) new.main = str[1]

        file <- "Output/PCAplot"
        if(!is.null(tag)) filenew = paste0(file,"_",tag,".pdf")
        else filenew <- "Output/PCAplot.pdf"

        pdf(filenew,height=5,width=len,useDingbats = F)
        if(ntype>1) par(mfrow=c(1,ntype+1))
        plot(xx[,1], xx[,2], xlab=paste0("PC1 (",round(dev[1], 1),"%)"), ylab = paste0("PC2 (",round(dev[2], 1),"%)"), main=paste0("PCA on ",new.main," data"),  pch=16)
        if(length(outlier)>0){
          points(xx[cid_out,1], xx[cid_out,2], col=2, pch=16, cex=1.1)
          text(xx[cid_out,1], xx[cid_out,2], col=2,labels =outlier, pos = 4)
        }
        if(ntype>1){
          plot(xx2[,1], xx2[,2], xlab=paste0("PC1 (",round(dev2[1], 1),"%)"), ylab = paste0("PC2 (",round(dev2[2], 1),"%)"), main=paste0("PCA on ", str[1]," data"), pch=16)
          if(length(outlier)>0){
            points(xx2[cid_out,1], xx2[cid_out,2], col=2, pch=16, cex=1.1)
            text(xx2[cid_out,1], xx2[cid_out,2], col=2,labels =outlier, pos = 4)
          }
          plot(xx3[,1], xx3[,2], xlab=paste0("PC1 (",round(dev3[1], 1),"%)"), ylab = paste0("PC2 (",round(dev3[2], 1),"%)"), main=paste0("PCA on ", str[2]," data"), pch=16)
          if(length(outlier)>0){
            points(xx3[cid_out,1], xx3[cid_out,2], col=2, pch=16, cex=1.1)
            text(xx3[cid_out,1], xx3[cid_out,2], col=2,labels =outlier, pos = 4)
          }
          if(ntype==3){
            plot(xx4[,1], xx4[,2], xlab=paste0("PC1 (",round(dev4[1], 1),"%)"), ylab = paste0("PC2 (",round(dev4[2], 1),"%)"), main=paste0("PCA on ", str[3]," data"), pch=16)
            if(length(outlier)>0){
              points(xx4[cid_out,1], xx4[cid_out,2], col=2, pch=16, cex=1.1)
              text(xx4[cid_out,1], xx4[cid_out,2], col=2,labels =outlier, pos = 4)
            }
          }
        }
        dev.off()
      }
    }
  }


  ########################
  ### Analysis options ###
  ########################
  ## Option 1: unsupervised approach (uses glasso to identify sparse precision matrix to derive adjacency matrix)
  ## Option 2: hybrid approach (fusion of networks)
    if(option==1|option==2){

    if(!useCorrelation) COVmat = cov(t(DAT),use = "pairwise.complete.obs")
    if(useCorrelation ) COVmat = cor(t(DAT),use = "pairwise.complete.obs")
    col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(20)
    type =NULL
    for(i in 1:ntype){
      type =c(type, rep(str[i], vec_n[i]))
    }
    COLpal = NULL
    colBAR = NULL

    if(ntype>1){
      CC = c("royalblue","cyan1","slateblue1")
      for(i in 1:ntype){
        #colBAR = c(colBAR,rep(str[i], vec_n[i]))
        COLpal = c(COLpal, rep(CC[i],vec_n[i]))
      }
    }
    rowBAR = colBAR = COLpal
    rr = quantile(as.numeric(COVmat))
    if(abs(rr[2])<0.4|abs(rr[4])<0.4) lk = 0.5
    if(abs(rr[2])>0.4|abs(rr[4])>0.4) lk = ceiling(rr[4])
    if(Plot.Covariance){
      file <- "Output/Heatmap_CrossCovarianceMatrix"

      if(!is.null(tag)) filenew = paste0(file,"_",tag,".png")
      if(is.null(tag)) filenew <- "Output/Heatmap_CrossCovarianceMatrix.png"

      if(!is.null(colBAR)){
        ttt= factor(KEY$DataType[match(row.names(COVmat), KEY$Label)],labels= unique(KEY$DataType))
        if(ntype==2) {
          CCnew = CC[1:2]
          ttt2 = ifelse(ttt=="X", str[1],str[2])
          names(CCnew) = str
          ann_col<-list(DataType=CCnew)
        }
        if(ntype==3){
          CCnew = CC
          ttt2 = rep(NA,length(ttt))
          for(i in 1:length(ttt)){
            if(ttt[i]=="X") ttt2[i] = str[1]
            if(ttt[i]=="Y") ttt2[i] = str[2]
            if(ttt[i]=="Z") ttt2[i] = str[3]
          }

          names(CCnew) = str
          ann_col=list(DataType=CCnew)
        }
        categoriesrow <- data.frame(DataType=factor(ttt2, labels=str, levels=str))
        row.names(categoriesrow) <- row.names(COVmat)

        HM_COVmat= pheatmap(COVmat, cluster_cols=T, cluster_rows=T, col=rev(col), breaks=seq(-lk,lk,by=(lk+lk)/20),filename = filenew, width=10, height=8.5,
                    annotation_names_row=F, show_rownames = F, show_colnames = F,annotation_colors = ann_col,annotation_names_col=F,
                    legend=T,treeheight_row=75,treeheight_col=50,annotation_row = categoriesrow, annotation_col =categoriesrow, main = "Cross-Covariance Matrix" ,fontsize=14)
      }
      if(is.null(colBAR)) HM_COVmat = pheatmap(COVmat, cluster_cols=T, cluster_rows=T, col=rev(col), breaks=seq(-lk,lk,by=(lk+lk)/20),
                                               show_rownames = F, show_colnames = F,legend=T,filename = filenew, width=10, height=8.5,
                                               treeheight_row=75,treeheight_col=50, main = "Cross-Covariance Matrix" ,fontsize=14)


    }
    if(useCorrelation) Snew = COVmat
    if(!useCorrelation){
      ### Given S is your covariance matrix ###
      if(!checkSPD(COVmat)){
        Snew = make.positive.definite(COVmat)
        Fnorm.diff = norm(Snew - COVmat, type = 'F')  ## you can use this to check the difference in the norm of the PSD corrected matrix
        cat("Input cross-covariance matrix is not positive definite.\nFinding nearest positive semi-definite matrix with a frobenius norm difference of", Fnorm.diff,"\n")
      }
      if(checkSPD(COVmat)) Snew = COVmat
    }

    cat("\nCarrying out graphical lasso, this may take a while...\n\n")
    n<- ncol(DAT)

    if(Calibration){
      cat("Starting calibration to pick an optimal lambda value.\n")
      if(is.null(lambda.vec)){
        cat("No lambda values were supplemented by user. Creating a lambda vector of length ",numLambda," using data.\n")
      }

      #lambda.vec=exp(seq(log(0.5),log(0.01), length=30))

      a = Sys.time()

      ModCriterion = glasso_Select(Snew, sample.size = n, rholist=lambda.vec,nlambda=numLambda,verbose=verbose, penalty=penalty, gamma = 0.5)
      b = Sys.time()
      timeDiff = b-a
      timePar = " minutes"
      if(timeDiff>60) {
        timeDiff=timeDiff/60
        timePar = " hour(s)"
      }
      cat(paste0("Completed in ",round(timeDiff,3), timePar, ", moving on to cross validation.\n"))
      CVtest = glasso_cv(t(DAT), nfolds=numFolds,cor=useCorrelation, rholist=lambda.vec,nlambda=numLambda,verbose=verbose)
      c = Sys.time()
      timeDiff = c-b
      timePar = " minutes"
      if(timeDiff>60) {
        timeDiff=timeDiff/60
        timePar = " hour(s)"
      }
      cat(paste0("Completed in ",round(timeDiff,3), timePar, "...\n"))

      if(is.null(lambda.vec)) lambda.vec = ModCriterion$rholist

      if(Plot.modelSelect){
        file <- "Output/CalibrationPlots_glasso"
        if(!is.null(tag)) filenew = paste0(file,"_",tag,".pdf")
        else filenew <- "Output/CalibrationPlots_glasso.pdf"

        pdf(filenew,height=12, width=11.5, useDingbats = F)
        par(mfrow=c(2,2))
        plot(lambda.vec, ModCriterion$AIC,lwd=1.2,pch=16, col="grey70", type="o", xlab="Lambda values", ylab="",main="AIC")
        pick = which.min(ModCriterion$AIC)
        #TP = detectTurningPoints(ModCriterion$AIC, "min")
        points(lambda.vec[pick], ModCriterion$AIC[pick], pch=4,cex=1.6,lwd=2,col="red")
        mtext(text = paste0("Lambda = ", round(lambda.vec[pick],3),"\nAIC = ",round(ModCriterion$AIC[pick],0)),side = 1,at = (lambda.vec[pick]),line = -2.5 ,cex=0.7)

        plot(lambda.vec, ModCriterion$BIC,lwd=1.2,pch=16, col="grey70", type="o", xlab="Lambda values", ylab="",
             main="BIC")
        pick = which.min(ModCriterion$BIC)
        points(lambda.vec[pick], ModCriterion$BIC[pick],  pch=4,cex=1.6,lwd=2,col="red")
        mtext(text = paste0("Lambda = ", round(lambda.vec[pick],3),"\nBIC = ",round(ModCriterion$BIC[pick],0)),side = 1,at = (lambda.vec[pick]),line = -2.5 ,cex=0.7)

        plot(lambda.vec, ModCriterion$eBIC,lwd=1.2,pch=16, col="grey70", type="o", xlab="Lambda values", ylab="",
             main="extended-BIC")
        pick = which.min(ModCriterion$eBIC)
        points(lambda.vec[pick], ModCriterion$eBIC[pick], pch=4,cex=1.6,lwd=2,col="red")
        mtext(text = paste0("Lambda = ", round(lambda.vec[pick],3),"\neBIC = ",round(ModCriterion$eBIC[pick],0)),side = 1,at = (lambda.vec[pick]),line = -2.5 ,cex=0.7)

        plot(CVtest$rholist, CVtest$CVvalue, cex=1.2,pch=16, col="grey70", type="o", xlab="Lambda values", ylab="CV(rho)",
             main=paste0(numFolds,"-fold Cross-validation"))
        pick =which(CVtest$rholist==CVtest$optLambda)
        points(CVtest$optLambda, CVtest$CVvalue[pick], pch=4,cex=1.6,lwd=2,col="red")
        mtext(text = paste0("Lambda = ", round(CVtest$optLambda,3),"\nCV value= ",round(CVtest$CVvalue[pick],0)),side = 3,at = (lambda.vec[pick]),line = -3 ,cex=0.7)
        dev.off()
      }
      if(is.null(optLambda)){

        optLambda_AIC = lambda.vec[which.min(ModCriterion$AIC)]
        optLambda_BIC = lambda.vec[which.min(ModCriterion$BIC)]
        optLambda_eBIC = lambda.vec[which.min(ModCriterion$eBIC)]
        optLambda_CV = lambda.vec[which.max(CVtest$CVvalue)]

        lastopt = lambda.vec[1]
        criterion_new = criterion
        if((criterion=="AIC") & (optLambda_AIC == lastopt)){
          if(optLambda_BIC != lastopt) criterion_new="BIC"
          if(optLambda_eBIC != lastopt) criterion_new="eBIC"
          if(optLambda_CV != lastopt) criterion_new="CV"
        }
        if((criterion=="BIC") & (optLambda_BIC == lastopt)){
          if(optLambda_CV != lastopt) criterion_new="CV"
          if(optLambda_AIC != lastopt) criterion_new="AIC"
          if(optLambda_eBIC != lastopt) criterion_new="eBIC"
        }
        if((criterion=="eBIC") & (optLambda_eBIC == lastopt)){
          if(optLambda_CV != lastopt) criterion_new="CV"
          if(optLambda_AIC != lastopt) criterion_new="AIC"
          if(optLambda_BIC != lastopt) criterion_new="BIC"
        }
        if((criterion=="CV") & (optLambda_CV == lastopt)){
          if(optLambda_BIC != lastopt) criterion_new="BIC"
          if(optLambda_eBIC != lastopt) criterion_new="eBIC"
          if(optLambda_AIC!= lastopt) criterion_new="AIC"
        }
        if(criterion!= criterion_new) cat("Specified criterion resulted in the largest lambda selected, modifying to use new criterion instead: ",criterion_new,"\n")
        criterion=criterion_new
        if(criterion=="AIC") optLambda = optLambda_AIC
        if(criterion=="BIC") optLambda = optLambda_BIC
        if(criterion=="eBIC") optLambda = optLambda_eBIC
        if(criterion=="CV") optLambda = optLambda_CV
      }
      cat(paste0("Using ",criterion," as the model selection criterion, we pick lambda = ", round(optLambda,3), " as the optimal threshold in glasso.\n\n"))

      return("End of Calibration. Please input a new lambda vector and re-calibration or re-run NetDeconvolute() after picking optimal lambda value and set Calibration to FALSE.\n")
    }
    cat("Skipping calibration step and using user-specified lambda value to estimate precision matrix.\n")

    GLP.final<-huge(Snew,method="glasso",verbose=T, scr=F,lambda =optLambda, cov.output=T)
    AdjacencyMat =GLP.final$path[[1]]
    pp = sum(AdjacencyMat[upper.tri(AdjacencyMat)]>0)/length(AdjacencyMat[upper.tri(AdjacencyMat)])

    cat(paste0(round(pp,3)*100,"% of all the possible edges were non-zeros in the sparse precision matrix.\n\n"))

    PRECISION = data.frame(GLP.final$icov)
    row.names(PRECISION) = colnames(PRECISION) = row.names(DAT)
    row.names(AdjacencyMat) = colnames(AdjacencyMat) = row.names(DAT)
    PRECISION_out = data.frame(row.names(PRECISION),row.names(PRECISION), PRECISION, check.names=F)
    colnames(PRECISION_out)[1:2] = c("Key","Feature")

    PRECISION_out$Key = KEY$NewName[match(row.names(PRECISION),KEY$Label)]
    PRECISION_out$Feature = KEY$Name[match(row.names(PRECISION),KEY$Label)]
    colnames(PRECISION_out)[3:ncol(PRECISION_out)] = KEY$NewName[match(row.names(PRECISION),KEY$Label)]

    PartialCor = convertPC(as.matrix(PRECISION))
    PartialCor_out = data.frame(row.names(PRECISION),row.names(PRECISION),PartialCor, check.names=F)
    colnames(PartialCor_out)[1:2] = c("Key","Feature")

    PartialCor_out$Key = KEY$NewName[match(row.names(PRECISION),KEY$Label)]
    PartialCor_out$Feature = KEY$Name[match(row.names(PRECISION),KEY$Label)]
    colnames(PartialCor_out)[3:ncol(PartialCor_out)] = KEY$NewName[match(row.names(PRECISION),KEY$Label)]

    write.table(PRECISION_out,"Output/glasso_estimated_icov.txt",sep="\t",row.names=F,quote=F)
    write.table(PartialCor_out,"Output/PartialCorrelation_icov.txt",sep="\t",row.names=F,quote=F)

    ordrow = HM_COVmat$tree_row$order
    ordcol = HM_COVmat$tree_col$order

    AdjacencyMatN = AdjacencyMat[ordrow,ordcol]
    PRECISION_N=PRECISION[ordrow,ordcol]
    rowBARN =rowBAR[ordrow]
    colBARN =colBAR[ordcol]

    COVmatN =COVmat[ordrow,ordcol]
    SnewN =Snew[ordrow,ordcol]


    if(Plot.Covariance & option!=2){

      col2 = c("black","red")

      arr = list(COVmatN, SnewN ,PRECISION_N, AdjacencyMatN)
      title = c("Sample covariance","Nearest semi-PD covariance",
                "Estimated Precision matrix (GLASSO)", "Adjacency matrix (GLASSO)")

      rr = quantile(Snew)
      tk1 = ceiling(abs(rr[2]))

      rr2 = checkRange(PRECISION_N)
      mm = min(abs(rr2[2]), abs(rr2[4]))
      if(mm>0.01) tk2 = 0.01
      if(mm<0.01) tk2 = 0.005
      tk3 = 1

      bks = list(seq(-tk1,tk1,by=tk1/10), seq(-tk1,tk1,by=tk1/10),
                 seq(-tk2,tk2,by=tk2/10),  seq(0,tk3,by=tk3/10))

      GL= list()
      for(i in 1:4){
        newMat = as.matrix(arr[[i]])
        rr = row.names(newMat)
        cc = colnames(newMat)
        tttrow= factor(KEY$DataType[match(rr, KEY$Label)],labels= unique(KEY$DataType))
        tttcol= factor(KEY$DataType[match(cc, KEY$Label)],labels= unique(KEY$DataType))
        if(ntype==2) {
          CCnew = CC[1:2]
          ttt2row = ifelse(tttrow=="X", str[1],str[2])
          ttt2col = ifelse(tttcol=="X", str[1],str[2])
          names(CCnew) = str
          ann_col<-list(DataType=CCnew)
        }
        if(ntype==3){
          CCnew = CC
          ttt2row=ttt2col = rep(NA,length(tttrow))
          for(j in 1:length(tttrow)){
            if(tttrow[j]=="X") ttt2row[j] = str[1]
            if(tttrow[j]=="Y") ttt2row[j] = str[2]
            if(tttrow[j]=="Z") ttt2row[j] = str[3]
            if(tttcol[j]=="X") ttt2col[j] = str[1]
            if(tttcol[j]=="Y") ttt2col[j] = str[2]
            if(tttcol[j]=="Z") ttt2col[j] = str[3]
          }
          names(CCnew) = str
          ann_col=list(DataType=CCnew)
        }
        categoriesrow <- data.frame(DataType=factor(ttt2row, labels=str, levels=str))
        row.names(categoriesrow) <- rr
        categoriescol <- data.frame(DataType=factor(ttt2col, labels=str, levels=str))
        row.names(categoriescol) <- cc


        if(i!=4){
          if(is.null(rowBARN)){
             h= pheatmap(newMat, cluster_cols=F, cluster_rows=F, col=bluered(20), breaks=bks[[i]],  dendrogram="none",mar=c(3,3),
                        legend=FALSE, show_rownames = F, show_colnames = F, main = title[i] ,fontsize=14,silent=T)
          }
          if(!is.null(rowBARN)){
              h= pheatmap(newMat, cluster_cols=F, cluster_rows=F, col=bluered(20), breaks=bks[[i]],  dendrogram="none",mar=c(3,3),
                      annotation_names_row=F,legend=FALSE, show_rownames = F, show_colnames = F,annotation_colors = ann_col,annotation_names_col=F,
                      annotation_row = categoriesrow, annotation_col =categoriescol, main = title[i] ,fontsize=14,silent=T)
          }
        }

       if(i==4){
          if(is.null(rowBARN)){
            h= pheatmap(newMat, cluster_cols=F, cluster_rows=F, col=rev(gray.colors(10)), breaks=bks[[i]],  dendrogram="none",mar=c(3,3),
                      legend=FALSE, show_rownames = F, show_colnames = F, main = title[i] ,fontsize=14,silent=T)
          }
          if(!is.null(rowBARN)){
            h= pheatmap(newMat, cluster_cols=F, cluster_rows=F, col=rev(gray.colors(10)), breaks=bks[[i]],  dendrogram="none",mar=c(3,3),
                      legend=FALSE, show_rownames = F, show_colnames = F, annotation_row = categoriesrow, annotation_col =categoriescol,annotation_names_col=F,
                      annotation_names_row=F,annotation_colors = ann_col,main = title[i] ,fontsize=14,silent=T)
          }
        }
        GL[[i]] = h[[4]]
      }

      file <- "Output/Plots_glasso"
      if(!is.null(tag)) filenew = paste0(file,"_",tag,".png")
      if(is.null(tag)) filenew <- "Output/Plots_glasso.png"
      png(filenew, width = 14, height =12, units = 'in', res = 300)
      par( mar=c(3,4,5,4))
      grid.arrange(grobs=GL, ncol=2, clip=TRUE)
      dev.off()

    }


    ## output integrated data ##
    Z_COMB = DAT
    Z_COMB = data.frame(KEY$Name[match(row.names(Z_COMB), KEY$Label)], Z_COMB, check.names=F)
    colnames(Z_COMB)[1] = "Feature"
    write.table(Z_COMB,"Output/Combined_data.txt",sep="\t",row.names=F,quote=F)

    NET_est = Network2Adj(PRECISION, dir = 2, ff = colnames(PRECISION),PartialCor = PartialCor,reportCor=T)
    mA = match(NET_est$nodeA, KEY$Label)
    mB = match(NET_est$nodeB, KEY$Label)
    NET_est$nodeA=KEY$Name[mA]
    NET_est$nodeB=KEY$Name[mB]

    NET_est$DatatypeA = str[as.numeric(factor(KEY$DataType[mA]))]
    NET_est$DatatypeB = str[as.numeric(factor(KEY$DataType[mB]))]
    NET_est$EdgeType = "Within"
    NET_est$EdgeType[which(NET_est$DatatypeA!=NET_est$DatatypeB)] = "Between"

    cat("\nIn the sparse covariance matrix, there are a total of",sum(NET_est$EdgeType=="Within"),"within-data types and",
        sum(NET_est$EdgeType=="Between"), "between-data types with conditional dependencies between features.\n")


    write.table(NET_est,"Output/Estimated_Network_glasso.txt",sep="\t",row.names=F,quote=F)
    system("cp Output/Combined_data.txt Output/Estimated_Network_glasso.txt iOmicsPASS/inputFiles/")
  }

  if(option==2){
    ## Option 2:semi-supervised approach combining glasso with known prior network to form a fused network

    Prior = NetworkFile[,c(1:3)]
    if(any(NetworkFile[,4]=="X"|NetworkFile[,5]=="X")){
      cid1_x = which(toupper(as.character(NetworkFile[,4]))=="X")
      cid2_x = which(toupper(as.character(NetworkFile[,5]))=="X")
      if(length(cid1_x)!=0) Prior[cid1_x,1] = KEYx$Label[match(NetworkFile[cid1_x,1],KEYx$Name)]
      if(length(cid2_x)!=0) Prior[cid2_x,2] = KEYx$Label[match(NetworkFile[cid2_x,2],KEYx$Name)]
    }
    if(any(NetworkFile[,4]=="Y"|NetworkFile[,5]=="Y")){
      cid1_y = which(toupper(as.character(NetworkFile[,4]))=="Y")
      cid2_y = which(toupper(as.character(NetworkFile[,5]))=="Y")
      if(length(cid1_y)!=0)Prior[cid1_y,1] = KEYy$Label[match(NetworkFile[cid1_y,1],KEYy$Name)]
      if(length(cid2_y)!=0) Prior[cid2_y,2] = KEYy$Label[match(NetworkFile[cid2_y,2],KEYy$Name)]
    }
    if(any(NetworkFile[,4]=="Z"|NetworkFile[,5]=="Z")){
      cid1_z = which(toupper(as.character(NetworkFile[,4]))=="Z")
      cid2_z = which(toupper(as.character(NetworkFile[,5]))=="Z")
      if(length(cid1_z)!=0)Prior[cid1_z,1] = KEYz$Label[match(NetworkFile[cid1_z,1],KEYz$Name)]
      if(length(cid2_z)!=0) Prior[cid2_z,2] = KEYz$Label[match(NetworkFile[cid2_z,2],KEYz$Name)]
    }
    Prior = Prior[!apply(Prior[,c(1:2)],1,function(x)any(is.na(x))),]

    PriorMat = Network2Adj(Prior, dir = 1, ff = colnames(PRECISION))

    tmpL = createMatrixB(COVmat, PRECISION, PriorMat)
    matrixA = tmpL$matrixA
    matrixB = tmpL$matrixB

    W = 0.5* (matrixA + PriorMat)
    U = unionBinMat(matrixA,PriorMat)
    phi = getsigns(matrixA,PriorMat)
    Smat = computeScore(U, phi,W, matrixB)

    Fusedmat = PriorMat + Smat
    FusedmatN =cutFusedmat(Fusedmat, cutoff=cutoff_fusedmat)
    Fmat = convertBinary(FusedmatN)

    PriorMatN =PriorMat[match(row.names(COVmatN),row.names(PriorMat)),match(colnames(COVmatN),colnames(PriorMat))]
    FusedmatN =Fusedmat[match(row.names(COVmatN),row.names(Fusedmat)),match(colnames(COVmatN),colnames(Fusedmat))]
    FmatN =Fmat[match(row.names(COVmatN),row.names(Fmat)),match(colnames(COVmatN),colnames(Fmat))]
    if(Plot.Covariance){
      arr = list(COVmatN, PRECISION_N ,AdjacencyMatN,
                  PriorMatN, FusedmatN, FmatN)

      title = c("Sample covariance","Estimated Precision matrix (glasso)",
                "Adjacency matrix (glasso)",  "Prior Matrix","Fused matrix (Hybrid)",
                "Adjacency matrix (Hybrid)")

      rr = quantile(COVmatN)
      tk1 = ceiling(abs(rr[2]))

      rr2 = checkRange(PRECISION_N)
      mm = min(abs(rr2[2]), abs(rr2[4]))
      if(mm>0.01) tk2 = 0.01
      if(mm<0.01) tk2 = 0.005
      tk3 = 1

      bks = list(seq(-tk1,tk1,by=tk1/10), seq(-tk2,tk2,by=tk2/10),
                 seq(0,tk3,by=tk3/10), seq(0,tk3,by=tk3/10),
                 seq(-tk3,tk3,by=tk3/10), seq(0,tk3,by=tk3/10))

      GL= list()
      for(i in 1:6){
        newMat = as.matrix(arr[[i]])
        rr = row.names(newMat)
        cc = colnames(newMat)

        tttrow= factor(KEY$DataType[match(rr, KEY$Label)],labels= unique(KEY$DataType))
        tttcol= factor(KEY$DataType[match(cc, KEY$Label)],labels= unique(KEY$DataType))
        if(ntype==2) {
          CCnew = CC[1:2]
          ttt2row = ifelse(tttrow=="X", str[1],str[2])
          ttt2col = ifelse(tttcol=="X", str[1],str[2])
          names(CCnew) = str
          ann_col<-list(DataType=CCnew)
        }
        if(ntype==3){
          CCnew = CC
          ttt2row=ttt2col = rep(NA,length(tttrow))
          for(i in 1:length(tttrow)){
            if(tttrow[j]=="X") ttt2row[j] = str[1]
            if(tttrow[j]=="Y") ttt2row[j] = str[2]
            if(tttrow[j]=="Z") ttt2row[j] = str[3]
            if(tttcol[j]=="X") ttt2col[j] = str[1]
            if(tttcol[j]=="Y") ttt2col[j] = str[2]
            if(tttcol[j]=="Z") ttt2col[j] = str[3]
          }
          names(CCnew) = str
          ann_col=list(DataType=CCnew)
        }
        categoriesrow <- data.frame(DataType=factor(ttt2row, labels=str, levels=str))
        row.names(categoriesrow) <- rr
        categoriescol <- data.frame(DataType=factor(ttt2col, labels=str, levels=str))
        row.names(categoriescol) <- cc

        if(i%in% c(1,2,5)){
          if(is.null(rowBARN)){
            h= pheatmap(newMat, cluster_cols=F, cluster_rows=F, col=bluered(20), breaks=bks[[i]],  dendrogram="none",mar=c(3,3),
                        legend=FALSE, show_rownames = F, show_colnames = F, main = title[i] ,fontsize=11,silent=T)
          }
          if(!is.null(rowBARN)){
            categoriesrow <- data.frame(DataType = factor(KEY$DataType[match(rr, KEY$Label)], labels = unique(KEY$DataType)))
            categoriescol <- data.frame(DataType = factor(KEY$DataType[match(cc, KEY$Label)], labels = unique(KEY$DataType)))
            row.names(categoriesrow) <- rr
            row.names(categoriescol) <- cc
            if(ntype==2) ann_col=list(DataType=c(X=CC[1], Y=CC[2]))
            if(ntype==3) ann_col=list(DataType=c(X=CC[1], Y=CC[2], Z=CC[3]))
            h= pheatmap(newMat, cluster_cols=F, cluster_rows=F, col=bluered(20), breaks=bks[[i]],  dendrogram="none",mar=c(3,3),
                        annotation_names_row=F,legend=FALSE, show_rownames = F, show_colnames = F,annotation_colors = ann_col,annotation_names_col=F,
                        annotation_row = categoriesrow, annotation_col =categoriescol, main = title[i] ,fontsize=11,silent=T)

          }
        }

        if(i%in% c(3,4,6)){
          if(is.null(rowBARN)){
            h= pheatmap(newMat, cluster_cols=F, cluster_rows=F, col=rev(gray.colors(10)), breaks=bks[[i]],  dendrogram="none",mar=c(3,3),
                        legend=FALSE, show_rownames = F, show_colnames = F, main = title[i] ,fontsize=11,silent=T)
          }
          if(!is.null(rowBARN)){
            h= pheatmap(newMat, cluster_cols=F, cluster_rows=F, col=rev(gray.colors(10)), breaks=bks[[i]],  dendrogram="none",mar=c(3,3),
                        legend=FALSE, show_rownames = F, show_colnames = F, annotation_row = categoriesrow, annotation_col =categoriescol,annotation_names_col=F,
                        annotation_names_row=F,annotation_colors = ann_col,main = title[i] ,fontsize=11,silent=T)
          }
        }
        GL[[i]] = h[[4]]
      }

      file <- "Output/Plots_Hybridmethod"
      if(!is.null(tag)) filenew = paste0(file,"_",tag,".png")
      if(is.null(tag)) filenew <- "Output/Plots_Hybridmethod.png"

      png(filenew, width=16.5, height=10, units="in", res=300)
      par(oma = c(1,1,1,1),mar=c(3,4,5,4))
      grid.arrange(grobs=GL, ncol=3, clip=TRUE)
      dev.off()
    }

    FusedNetwork = Network2Adj(Fmat, dir = 2,ff = colnames(Fmat), PartialCor = COVmat,reportCor=F)

    mA = match(FusedNetwork$nodeA, KEY$Label)
    mB = match(FusedNetwork$nodeB, KEY$Label)
    FusedNetwork$nodeA=KEY$Name[mA]
    FusedNetwork$nodeB=KEY$Name[mB]

    FusedNetwork$DatatypeA = str[as.numeric(factor(KEY$DataType[mA]))]
    FusedNetwork$DatatypeB = str[as.numeric(factor(KEY$DataType[mB]))]
    FusedNetwork$EdgeType = "Within"
    FusedNetwork$EdgeType[which(FusedNetwork$DatatypeA!=FusedNetwork$DatatypeB)] = "Between"

    cat("\nIn the fused matrix, there are a total of",sum(FusedNetwork$EdgeType=="Within"),"within-data types and",
        sum(FusedNetwork$EdgeType=="Between"),"between-data types with conditional dependencies between features.\n")

    checkOverlaps(NET_est,FusedNetwork)

    write.table(FusedNetwork,"Output/Fusednetwork_hybridmethod.txt",sep="\t",row.names=F, quote=F)
    system("cp Output/Combined_data.txt Output/Fusednetwork_hybridmethod.txt iOmicsPASS/inputFiles/")
  }
}

#' Create Input parameter file for iOmicsPASS
#'
#' @param data.X dataframe X
#' @param data.Y dataframe Y
#' @param data.Z dataframe Z
#' @param dir directory pointing to the location of the files
#' @param phenotype file with group sample information
#' @param pathway pathway file for network enrichment
#' @param within.net network file for creating edges within X
#' @param btw.net network file for creating edges between X and Y
#' @param standardize.data whether to standardize each data (default=TRUE)
#' @param log.transform whether to log-transform each data (default=FALSE)
#' @param normalizeBy Either "Y" or "Z". if "Y", coexpression scores in X are normalized by Y and if "Z", coexpression scores in X are normalized by Z.
#' @param usePrior whether to use an input prior to classify samples in discriminant model. If false, equal prior will be used.(default=FALSE)
#' @param priorfile file of the prior probabilities calculated from createPrior().
#' @param Enrichment Whether to run network enrichment (default=TRUE)
#' @param Cross.Validate Whether to run cross-validation (default=TRUE)
#' @param num.Kfold number of folds in CV
#' @param knn.impute whether to perform K-nearest neighbor imputation for missing entries
#' @param knn.k number of folds for K-nearest neighbor imputation
#' @param max.block.impute number of blocks of samples to consider in KNN imputation
#' @param min.obs minimum number of non-missing observations requires across each feature in each phenotypic outcome.
#' @param min.prop minimum proportion of non-missing observations requires across each feature in each phenotypic outcome.
#' @param min.thres threshold to be used to derive the shrunken centroids in iOmicsPASS
#' @param tag a string that will be appended to the end of the output files like a tag.
#' @param bg.prop minimum proportion of features in each pathway that are also present in the background list.
#' @param min.bg.size minimum number of features in each pathway that are also present in the background list.
#' @param min.sig.size minimum number of features in each pathway that are both part of the signature and present in the background list.
#' @return a parameter file for input into iOmicsPASS.R()
#' @examples
#' data(PhenotypeFile)
#' data(bioPathways)
#'
#' ## Running with estimated network from NetDeconvolute() ##
#' createInputParam(data.X="Combined_data.txt", within.net="Estimated_Network_glasso.txt",
#' phenotype =PhenotypeFile, Enrichment=FALSE)
#'
#' ## with Network-based pathway Enrichment ##
#' createInputParam(data.X="Combined_data.txt", within.net="Estimated_Network_glasso.txt",
#' phenotype =PhenotypeFile, pathway =bioPathways, Enrichment=TRUE)
#'
#' data(Tulip_Protein)
#' data(Tulip_microRNA)
#' data(PPI_network)
#' data(TargetScan_network)
#'
#' ## Using known biological networks ##
#' createInputParam(data.X=Tulip_Protein,data.Y=Tulip_microRNA,phenotype =PhenotypeFile,
#' log.transform=TRUE,btw.net=TargetScan_network,within.net=PPI_network, Enrichment=FALSE)
#'
#'
#' @export
createInputParam <- function(data.X, data.Y=NULL, data.Z=NULL,within.net=NULL, btw.net=NULL, dir = "iOmicsPASS/inputFiles/",
                             phenotype =NULL, pathway = NULL,standardize.data = T, log.transform = F, normalizeBy =NULL,
                             min.obs = 5, min.prop = 0.8, knn.impute = F,knn.k = 10, max.block.impute = 1000,
                             Cross.Validate = TRUE, num.Kfold = 5, min.thres = NULL, Enrichment=TRUE, tag= NULL,
                             usePrior = FALSE, priorfile = NULL,bg.prop = 0.2,min.bg.size = 3, min.sig.size = 1){

  str.datX = "dataX.txt"
  str.datY = "dataY.txt"
  str.datZ = "dataZ.txt"
  str.phenotype = "subtypeF.txt"
  str.pathway = "pathwayF.txt"
  str.withinNet = "network_Within.txt"
  str.betweenNet = "network_Between.txt"
  str.priorfile = "PriorProb.txt"

  ## log.transform can be a vector of TRUE or FALSE corresponding to X,Y,Z
  ## standardize.data can be a vector of TRUE or FALSE corresponding to X,Y,Z


  bool.impute = "false"
  if(knn.impute) bool.impute ="true"
  bool.CV = "false"
  if(Cross.Validate) bool.CV ="true"
  bool.analyzeY = "false"
  if(!is.null(data.Y)) bool.analyzeY = "true"
  bool.analyzeZ = "false"
  if(!is.null(data.Z)) bool.analyzeZ = "true"
  bool.enrich = "false"
  if(Enrichment) bool.enrich = "true"
  bool.prior = "false"
  if(usePrior) bool.prior = "true"

  bool.stand = rep("false",3)
  if(length(standardize.data)==1) if(standardize.data) bool.stand = rep("true",3)
  if(length(standardize.data)>1) {
    for(i in 1:length(standardize.data))if(standardize.data[i]) bool.stand[i] ="true"
  }
  bool.logtrans = rep("false",3)
  if(length(log.transform)==1) if(log.transform) bool.logtrans = rep("true",3)
  if(length(log.transform)>1) {
    for(i in 1:length(log.transform))if(log.transform[i]) bool.logtrans[i] ="true"
  }

  if(is.null(data.X))stop("Data X required. Please try again!\n")

  if(is.null(phenotype))stop("Phenotype information is required. Please try again!\n")

  if(is.null(pathway) & Enrichment==TRUE)stop("Pathway information is required for pathway enrichment. Please try again!\n")

  if(!is.null(data.Y) & is.null(btw.net))stop("A between-network file connecting data X and Y is required. Please try again!\n")

  if(!is.null(data.Z) & is.null(normalizeBy))stop("Data Z is provided, please specify which data (X or Y) to normalize data Z by.\n")

  if(is.null(within.net)&is.null(btw.net))stop("At least one network information is required. Please try again!\n")

  if(class(data.X)== "character") str.datX = data.X
  if(class(data.X)== "data.frame"|class(data.X)== "matrix") write.table(data.X,paste0(dir,str.datX), sep="\t",row.names=F,quote=F)
  NAMES = c("DATA_X","ZTRANS_X","LOG_TRANSFORM_X")
  VALUES = c(paste0(dir,str.datX), bool.stand[1], bool.logtrans[1])

  if(class(phenotype)== "character") str.phenotype = phenotype
  if(class(phenotype)== "data.frame"|class(phenotype)== "matrix") write.table(phenotype,paste0(dir,str.phenotype), sep="\t",row.names=F,quote=F)
  NAMES = c(NAMES,"SUBTYPE_FILE")
  VALUES = c(VALUES,paste0(dir,str.phenotype))

  if(Enrichment){
    if(class(pathway)== "character") str.pathway = pathway
    if(class(pathway)== "data.frame"|class(pathway)== "matrix") write.table(pathway,paste0(dir,str.pathway), sep="\t",row.names=F,quote=F)
    NAMES = c(NAMES,"MODULE_FILE")
    VALUES = c(VALUES,paste0(dir,str.pathway))
  }
  if(!is.null(data.Y)){
    if(class(data.Y)== "character") str.datY = data.Y
    if(class(data.Y)== "data.frame"|class(data.Y)== "matrix")write.table(data.Y,paste0(dir,str.datY), sep="\t",row.names=F,quote=F)
    NAMES = c(NAMES,"DATA_Y","ZTRANS_Y","LOG_TRANSFORM_Y")
    VALUES = c(VALUES,paste0(dir,str.datY), bool.stand[2], bool.logtrans[2])
  }
  if(!is.null(data.Z)){
    if(class(data.Z)== "character") str.datZ = data.Z
    if(class(data.Z)== "data.frame"|class(data.Z)== "matrix")write.table(data.Z,paste0(dir,str.datZ), sep="\t",row.names=F,quote=F)
    NAMES = c(NAMES,"DATA_Z","ZTRANS_Z","LOG_TRANSFORM_Z")
    VALUES = c(VALUES,paste0(dir,str.datZ), bool.stand[3], bool.logtrans[3])
  }
  if(!is.null(within.net)){
    if(class(within.net)== "character") str.withinNet = within.net
    if(class(within.net)== "data.frame"|class(within.net)== "matrix")write.table(within.net,paste0(dir,str.withinNet), sep="\t",row.names=F,quote=F)
    NAMES = c(NAMES,"NETWORK_WITHIN")
    VALUES = c(VALUES,paste0(dir,str.withinNet))
  }
  if(!is.null(btw.net)){
    if(class(btw.net)== "character") str.betweenNet = btw.net
    if(class(btw.net)== "data.frame"|class(btw.net)== "matrix")write.table(btw.net,paste0(dir,str.betweenNet), sep="\t",row.names=F,quote=F)
    NAMES = c(NAMES,"NETWORK_BTW")
    VALUES = c(VALUES,paste0(dir,str.betweenNet))
  }

  if(!is.null(priorfile)){
    if(class(priorfile)== "character") str.priorfile = priorfile
    if(class(priorfile)== "data.frame"|class(priorfile)== "matrix") write.table(priorfile,paste0(dir,str.priorfile), sep="\t",row.names=F,quote=F)
    NAMES = c(NAMES,"PRIOR_FILE")
    VALUES = c(VALUES,paste0(dir,str.priorfile))
  }


  NAMES = c(NAMES, "KNN_IMPUTE","KNN_K","MAX_BLOCKSIZE", "CROSS_VALIDATION","CV_FOLD","BACKGROUND_PROP","ENRICHMENT","MINBG_SIZE","MINSIG_SIZE","ANALYZE_Y", "ANALYZE_Z","USE_PRIOR")
  VALUES = c(VALUES, bool.impute,knn.k, max.block.impute, bool.CV,num.Kfold, bg.prop, bool.enrich, min.bg.size, min.sig.size, bool.analyzeY, bool.analyzeZ, bool.prior)
  if(!is.null(normalizeBy)){
    NAMES = c(NAMES,"NORMALIZE_ZBY")
    VALUES = c(VALUES,normalizeBy)
  }
  if(!is.null(min.thres)){
    NAMES = c(NAMES,"MIN_THRES")
    VALUES = c(VALUES,min.thres)
  }
  if(!is.null(min.prop)){
    NAMES = c(NAMES,"MIN_PROP")
    VALUES = c(VALUES,min.prop)
  }
  if(!is.null(min.obs)){
    NAMES = c(NAMES,"MIN_OBS")
    VALUES = c(VALUES,min.obs)
  }
  input_file_dataframe <- data.frame(param_name=NAMES, param_value=VALUES,stringsAsFactors=FALSE)
  fnew = "input_param"
  if(!is.null(tag)) fnew = paste0(fnew,"_",tag)
  write.table(input_file_dataframe, file=paste0("iOmicsPASS/", fnew), sep=" = ", quote=F, col.names = F,row.names=F)

}

#' Carry out Predictive Analysis of Subnetwork Signatures
#'
#' @param ff input parameter
#' @param plotCV Whether or not to plot the performance of the cross-validation (default = TRUE)
#' @param Cross.Validate Whether or not cross-validation was performed in iOmicsPASS (default = TRUE)
#' @param outputDir directory to write the output files (default output to "iOmicsPASS/Output/")
#' @return dataframe with iOmicsPASS parameters used for prediction, multiple text files and plots written to output directory.
#' \itemize{
#'   \item AttributesTable.txt - Attributes table for every node in the network to be used in Cytoscape.
#'   \item BGlist.txt - A list of edges that are present in both the network and the input data.
#'   \item CVerrors.txt - A table of grid of threshold and their corresponding mean cross-validation errors and selected edges.
#'   \item CVplot_Penalty.pdf - A plot of the mean cross-validation error against the grid of threshold used to shrink the centroid.
#'   \item EdgesSelected_minThres.txt - dataframe with selected predictive edges and the dik scores for each phenotype that can be used to visualize directly in Cytoscape.
#'   \item Features_Neighbors.txt - Attributes table with neighboring node information for each node in the network.
#'   \item PredictiveEdges_Parameters.txt - A dataframe with with iOmicsPASS parameters used for prediction.
#'   \item SampleClass_Probabilities.txt - A dataframe with class probabilities assigned to the samples input.
#'   \item Ztransform_XXX.txt - Standardized data X/Y/Z.
#'   \item XXX_Enrichment_up.txt - a dataframe with the results of the enrichment of the selected edges (over-expressed) out of all the edges in the network for each class.
#'   \item XXX_Enrichment_down.txt - a dataframe with the results of the enrichment of the selected edges (under-expressed) out of all the edges in the network for each class.
#' }
#' @examples
#' data(PhenotypeFile)
#' data(bioPathways)
#'
#' ## Running with estimated network from NetDeconvolute() ##
#' createInputParam(data.X="Combined_data.txt", within.net="Estimated_Network_glasso.txt",
#' phenotype =PhenotypeFile,Enrichment=FALSE)
#' iOmicsPASS.output<-iOmicsPASS.R(ff="input_param")
#'
#' ## pick optimal threshold as 2.4 and turn off CV ##
#' ## rerun with network-based pathway enrichment ##
#' createInputParam(data.X="Combined_data.txt", within.net="Estimated_Network_glasso.txt",
#' phenotype=PhenotypeFile,pathway=bioPathways,Enrichment=TRUE,min.thres=2.4,Cross.Validate=FALSE)
#'
#' iOmicsPASS.output<-iOmicsPASS.R(ff="input_param",Cross.Validate=FALSE, plotCV=FALSE)
#'
#' data(Tulip_Protein)
#' data(Tulip_microRNA)
#' data(PPI_network)
#' data(TargetScan_network)
#'
#' ## Using known biological networks ##
#' createInputParam(data.X=Tulip_Protein,data.Y=Tulip_microRNA,phenotype=PhenotypeFile,
#' log.transform=TRUE,within.net=PPI_network,btw.net=TargetScan_network,Cross.Validate=TRUE,
#' Enrichment=FALSE, tag="KnownNetwork")
#' iOmicsPASS.output<-iOmicsPASS.R(ff="input_param_KnownNetwork")
#'
#'
#' @export
iOmicsPASS.R <- function( ff= "input.param", outputDir ="iOmicsPASS/Output/", Cross.Validate=TRUE, plotCV=TRUE){

  ret = NULL
  cmdStr = paste0("bin/iOmicsPASS iOmicsPASS/",ff)
  cmdStr =paste(cmdStr ,outputDir)
  system(cmdStr)

  if(plotCV) Plot_CrossValidation(outputDir)

  ## read in the required model parameters ##
  if(!Cross.Validate) ret = read.delim(paste0(outputDir,"/PredictiveEdges_Parameters.txt"), as.is=T,header=T,check.names=F)
  ret
}

#' Create prior probabilities for discriminant model in iOmicsPASS
#'
#' @param FILE dataframe with sample ID in the row names and clinical characteristics (for adjustment) in the columns.
#' @param y character string of the phenotypic or group variable.
#' @param tag a string that will be appended to the end of the output files like a tag.
#' @param var.cat a string or a vector of string matching the column names of FILE that are categorical variables.
#' @param testFile dataframe for test samples matching the colnames and variables in FILE
#' @param predict Whether or not to predict probabilities on a test dataset
#' @param outputDir directory to write the outputfile
#' @return a dataframe with class probabilities assigned to each sample, used as prior in iOmicsPASS
#' @examples
#' data(PhenotypeFile)
#' PhenotypeFile2 = PhenotypeFile
#' row.names(PhenotypeFile2) = PhenotypeFile$TulipID
#' PhenotypeFile2=PhenotypeFile2[,-1]
#' prior_out = createPrior(PhenotypeFile2, y = "Group",outputDir = "iOmicsPASS/inputFiles/")
#' prior_test = createPrior(PhenotypeFile2, y = "Group",predict=TRUE, testFile = TestData)
#' @export
createPrior <-function(FILE, y = NULL,var.cat=NULL,testFile=NULL,predict=FALSE, outputDir = "iOmicsPASS/inputFiles/", tag=NULL){

  #require(nnet)
  ## if predict is set to TRUE, model is built on input data and prior is created by predicting on test data
  ## if predict is set to FALSE, model is built on input data and prior is created for training data
  index_grp = which(colnames(FILE)==y)
  if(predict)if(is.null(testFile)) stop("Test data should be provided if predict is set to TRUE, please re-specify again!\n")


  if(length(index_grp)==0) stop("Outcome variable y is not found in the file, please re-specify again!\n")
  if(is.null(var.cat)) print("No categorical variables were specified, data will be all analyzed as continuous.\n")

  y.grp= factor(FILE[,index_grp])
  num.cat = length(unique(y.grp))
  if(num.cat==1) stop("Error! There is only 1 level of outcome y found.")

  FILE_new = FILE[,-index_grp]
  cid_cat= which(colnames(FILE_new)%in% var.cat)
  if(!is.null(var.cat)){
    if(length(cid_cat)==0) stop("The categorical variables were not found in the file, please re-specify again!\n")
    for(i in cid_cat) FILE_new[,i] = factor(FILE_new[,i])
  }
  if(predict){
    if(!all(colnames(testFile)%in%colnames(FILE_new))) stop("Not all the variables in input file is found in the test file. Please input same variable names in the test file for prediction.\n")
    testFile_new = testFile[,match(colnames(FILE_new),colnames(testFile))]
    if(!is.null(cid_cat))for(i in cid_cat) testFile_new[,i] = factor(testFile_new[,i])

  }

  if(num.cat==2){
    print("There are 2 level of outcome y found. A binary Logistic regression will be employed.\n")
    mod = glm(y.grp~.,family=binomial(link="logit"),data=FILE_new)
    if(!predict){
      y.pred =predict(mod, newdata = FILE_new,type= "response")
      y.pred =data.frame(1-y.pred, y.pred)
      colnames(y.pred) = unique(y.grp)
    }
    if(predict){
      y.pred =predict(mod, newdata = testFile_new,type= "response")
      y.pred =data.frame(1-y.pred, y.pred)
      colnames(y.pred) = unique(y.grp)
    }
    #predClass = apply(y.pred,1,function(x) colnames(y.pred)[which.max(x)])
  }

  if(num.cat>2){
    print(paste0("There are ",num.cat," level of outcome y found. A multinomial logistic regression will be employed.\n"))

    mod = multinom(y.grp ~ ., data = FILE_new)
    if(!predict){
      y.pred =predict(mod, newdata = FILE_new, "probs")
    }
    if(predict){
      y.pred =predict(mod, newdata = testFile_new, "probs")
    }
  }

  if(!predict){
    OUT = data.frame(row.names(FILE), y.pred, check.names=F)
    colnames(OUT)[1] = "ID"
    if(is.null(tag)) FF= "PriorProb.txt"
    if(!is.null(tag)) FF = paste0("PriorProb_", tag,".txt")
    write.table(OUT,paste0(outputDir,FF), sep="\t",row.names=F,quote=F)
  }
  if(predict){
    OUT = data.frame(row.names(testFile), y.pred, check.names=F)
    colnames(OUT)[1] = "ID"
    if(is.null(tag)) FF= "PriorProb_testData.txt"
    if(!is.null(tag)) FF = paste0("PriorProb_testData", tag,".txt")
    write.table(OUT,paste0(outputDir,FF), sep="\t",row.names=F,quote=F)
  }

  return(OUT)
}

#' Compile and build the C++ program for running iOmicsPASS
#'
#' @param currDir Current directory where iOmicsPASSv2plus is unzipped and where the makefile is.
#' @return NULL
#' @examples
#' # set your working directory to iOmicsPASSv2plus folder #
#' INSTALL.iOmicsPASS()
#'
#' @export
INSTALL.iOmicsPASS<- function(currDir=NULL){
  if(is.null(currDir)) currDir = getwd()
  system(paste0("cd ",currDir))

  if (!file.exists("bin")){
    dir.create("bin")
  }

  if (!file.exists("Input")){
    dir.create("Input")
  }

  if (!file.exists("Output")){
    dir.create("Output")
  }
  if (!file.exists("iOmicsPASS/inputFiles")){
    dir.create("iOmicsPASS/inputFiles")
  }
  if (!file.exists("iOmicsPASS/Output")){
    dir.create("iOmicsPASS/Output")
  }

  system("make")
}

#' Plot Cross-validation performance
#'
#' @param outputdir directory pointing to CVerrors.txt.
#' @return a PDF file with the performance plot
#' @examples
#' Plot_CrossValidation("iOmicsPASS/output/")
#' @export
Plot_CrossValidation<-function(outputdir){
  CVdat = read.delim(paste0(outputdir,"CVerrors.txt"),as.is=T,check.names=F)

  minError = min(CVdat$CVerror)
  minThres = CVdat$Threshold[which(CVdat$CVerror==minError)][1]
  genesurv = CVdat$NumEdgesSelected[which(CVdat$CVerror==minError)][1]

  cid = grep("CVerror_", colnames(CVdat))
  ll=gsub("CVerror_","",colnames(CVdat)[cid])
  cc = 2:(length(cid)+1)
  legend_lab = c("Overall")
  for(i in 1:length(ll)) legend_lab = c(legend_lab, ll[i])

  sd = sd(CVdat$CVerror[-nrow(CVdat)])
  within = CVdat$Threshold[which(abs(CVdat$CVerror)< (minError+sd))]
  xpt = max(within)
  ypt = CVdat$CVerror[which(CVdat$Threshold==xpt)]
  zpt = CVdat$MeanCV_NumEdgesSelected[which(CVdat$Threshold==xpt)]

  pdf(paste0(outputdir,"CVplot_Penalty.pdf"),height=5, width=7, useDingbats = F)
  par(mar=c(4,4,7,2),mai=c(1,1,1,0.5))
  plot(CVdat$Threshold,CVdat$CVerror,ylim=c(0,1),type="l",col=1,lwd=2.5, cex.axis=0.8,ylab="Mean misclassification error", xlab="Threshold")
  mtext("Mean number of edges selected",side=3,line=3)
  axis(side = 3,at=CVdat$Threshold, lab=CVdat$MeanCV_NumEdgesSelected,las=2,srt= 45,cex.axis=0.8)
  for(i in 1:length(cid)) lines(CVdat$Threshold, CVdat[,cid[i]], col=cc[i],lty=1)
  legend("topleft",legend_lab, col=c(1,cc),lty=1 ,lwd=c(4,rep(4,length(ll))),cex=0.8)

  points(minThres,CVdat$CVerror[match(minThres,CVdat$Threshold)], pch=4,lwd=2,col=2)
  points(xpt,ypt, pch=4,lwd=2,col=2)
  text(minThres,CVdat$CVerror[match(minThres,CVdat$Threshold)],pos=1, offset=1,lab=paste0("Minimum Threshold = ",paste(round(minThres,3),collapse="\t"),"\nClassification Error = ",round(minError,3),"\nSelected Edges = ",paste(genesurv,collapse="\t")),cex=0.6,las=1)
  if(xpt!=minThres) text(xpt, ypt,pos=1, offset=1,lab=paste0( "Threshold = ",paste(round(xpt,3),collapse="\t"),"\nClassification Error = ",round(ypt,3),"\nSelected Edges = ",paste(zpt,collapse="\t")),cex=0.6,las=1)
  abline(h=(minError+sd),lty=2,col="grey40")
  text(-0.15,(minError+sd)+0.02, pos =4, font =2,lab="1 SD above minimum threshold",cex=0.7,col="grey40")
  dev.off()
}

#' Classification on external dataset
#'
#' @param file output from running iOmicsPASS.R().
#' @param newData Test data where the various -omics data are concatenated over the same samples. First column should be the list of features to be matched to the predictive signatures.
#' @param standardize whether to perform standardization on the data.
#' @param usePrior boolean.
#' @param prior filename of prior.
#' @param prop proportion of features in test data that are part of the predictive signature. (default=0.8)
#' @return a dataframe with the class probabilities for each sample
#' @examples
#'
#' data(Tulip_Protein)
#' data(Tulip_microRNA)
#' data(PhenotypeFile)
#' data(PPI_network)
#' data(TargetScan_network)
#'
#' row.names(Tulip_Protein) = Tulip_Protein$Protein
#' row.names(Tulip_microRNA) = Tulip_microRNA$miRNA
#'
#'
#' ## create a testData by randomly sampling from original data ##
#' set.seed(22)
#' sample_pick = sample(PhenotypeFile$TulipID, 6)   ## 4 LIS and 2 OIR
#'
#' Tulip_Protein_test = Tulip_Protein[,c(1,match(sample_pick, colnames(Tulip_Protein)))]
#' Tulip_microRNA_test = Tulip_microRNA[,c(1,match(sample_pick, colnames(Tulip_microRNA)))]
#'
#' PhenotypeFile_test = PhenotypeFile[which(PhenotypeFile$TulipID %in% sample_pick),]
#'
#' row.names(PhenotypeFile_test) = PhenotypeFile_test$TulipID
#' row.names(PhenotypeFile) = PhenotypeFile$TulipID
#'
#' prior_train = createPrior(PhenotypeFile[,-1], y = "Group",predict=FALSE)
#'
#' createInputParam(data.X=Tulip_Protein, data.Y=Tulip_microRNA, within.net=PPI_network,
#'  btw.net=TargetScan_network,log.transform=TRUE,phenotype =PhenotypeFile,
#'  Enrichment=FALSE,usePrior = TRUE, priorfile = prior_train, tag="train")
#'
#' iOmicsPASS.train <- iOmicsPASS.R(ff = "input_param_train")
#'
#' createInputParam(data.X=Tulip_Protein, data.Y=Tulip_microRNA, within.net=PPI_network,
#'  btw.net=TargetScan_network,log.transform=TRUE,phenotype =PhenotypeFile,Cross.Validate=FALSE,
#'  min.thres=1.8,Enrichment=FALSE,usePrior = TRUE, priorfile = prior_train, tag="train")
#'  iOmicsPASS.train <- iOmicsPASS.R(ff = "input_param_train", plotCV=FALSE)
#'
#'  ## predict test data using signatures from training data ##
#' testData=rbind(Tulip_Protein_test[,-1],Tulip_microRNA_test[,-1])
#' testData=log2(testData)
#' testData = data.frame(row.names(testData), testData, check.names=FALSE)
#' colnames(testData)[1] = "Feature"
#'
#'  ## remove group information ##
#' PhenotypeFile_test=PhenotypeFile_test[,-2]
#' prior_test = createPrior(PhenotypeFile[,-1], y = "Group",predict=TRUE,
#' testFile = PhenotypeFile_test[,-1])
#' pred.out <- Predict.iOmicsPASS(iOmicsPASS.train , testData, standardize = TRUE,prop = 0.9,
#'  usePrior=T,prior=prior_test )
#'
#' @export
Predict.iOmicsPASS <- function(file, newData ,standardize = TRUE, usePrior=FALSE, prior=NULL,prop = 0.8){

  MAT = as.matrix(newData[,-1])
  row.names(MAT) = newData[,1]
  if(usePrior)if(is.null(prior)) stop("Prior should be provided if usePrior is set to TRUE, please run createPrior() with a test data!\n")
  if(standardize){
    mu_vec = apply(MAT,1,function(x) mean(x,na.rm=T))
    sd_vec = apply(MAT,1,function(x) sd(x,na.rm=T))
    MAT = sweep(MAT,1,mu_vec)
    for(k in 1:ncol(MAT)) MAT[,k] = MAT[,k]/sd_vec[k]
  }

  allG = unique(c(file$GeneA, file$GeneB))
  commonF = intersect(allG,row.names(MAT))
  pp = length(commonF)/length(allG)
  if(pp<prop) stop("Test data fails to satisfy the specified proprotion of features required to accurately assign samples to phenotypic groups!\n")
  if(!all(allG%in% row.names(MAT))){
    cat("Not all the required features are found in the Test data to create co-expression scores.\n")
    cat("Warning: The discriminant model will be built based on ", length(commonF)," out of ",length(allG), " features.\n")
    MAT_sub = MAT[match(commonF,row.names(MAT)),]
    cid1 = which(file$GeneA %in% commonF)
    cid2 = which(file$GeneB %in% commonF)
    file_sub = file[intersect(cid1,cid2),]
  }
  if(all(allG%in%row.names(MAT))){
    MAT_sub = MAT
    file_sub = file
  }
  SCORE= createCoexpressionScore(file_sub, MAT_sub)
  gg = gsub("_Xikbar","",colnames(file_sub)[-c(1:6)])
  XBARik = file_sub[,-c(1:6)]
  if(usePrior){
    priorM = prior[match(colnames(SCORE),prior$ID),]
    priorM = priorM[,match(gg,colnames(priorM))]
    row.names(priorM) = colnames(SCORE)
  }
  if(!usePrior) pp = 1/length(gg)

  denom = (file_sub$Si + file_sub$So)^2
  DISCRIMINANT = NULL
  for(i in 1:ncol(SCORE)){
    ss = SCORE[,i]
    dd = NULL
    for(j in 1:ncol(XBARik)){
        num = (ss - XBARik[,j])^2
        if(!usePrior){
          priorVal = pp
          if(priorVal==0)priorVal= 0.0000000001
          score = sum(num/denom, na.rm = T) - (2*log(pp))
        }
        if(usePrior){
          priorVal = priorM[i,j]
          if(priorVal==0)priorVal= 0.0000000001
          score = sum(num/denom, na.rm = T) - (2*log(priorVal))
        }
        dd = c(dd, score)
    }

    DISCRIMINANT =rbind(DISCRIMINANT, dd)
    colnames(DISCRIMINANT) = gg
  }

  Min = apply(DISCRIMINANT,1, min)
  DISCRIMINANT_new = sweep(DISCRIMINANT,1, FUN="-",Min)
  PROB = exp(-0.5*DISCRIMINANT_new)

  probSUM = apply(PROB,1,function(x) sum(x,na.rm=T))
  PROB = PROB/probSUM

  colnames(DISCRIMINANT) = paste0(gg,"_DiscriminantScore")
  colnames(PROB) = paste0(gg,"_ClassProbability")
  OUT = data.frame(colnames(SCORE), DISCRIMINANT,PROB, check.names=F)
  colnames(OUT)[1] = "ID"
  OUT$AssignedClass=""
  for(i in 1:nrow(OUT)){
    tmpVec =as.numeric(OUT[i,c((2+length(gg)):(1+2*length(gg)))])
    outcome =checkMax(tmpVec, gg)
    if(length(outcome)==1) OUT$AssignedClass[i] = outcome
    if(length(outcome)==0) OUT$AssignedClass[i] = NA
  }
  OUT
}



