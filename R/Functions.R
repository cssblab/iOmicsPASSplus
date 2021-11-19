#################
### Functions ###
#################

checkSPD<-function(MAT){
  ret = TRUE
  eigVal = eig(MAT)
  if(any(eigVal<0)) ret = FALSE
  ret
}

checkRange<- function(MAT){
  MAT = as.matrix(MAT)
  vec = as.numeric(MAT)
  cid_zeros = which(vec==0)
  if(length(cid_zeros)>0) vec = vec[-cid_zeros]
  quantile(vec)
}

cvFolds<-function (n, K = 5, R = 1, type = c("random", "consecutive", "interleaved")){
  if (!isTRUE((K > 1) && K <= n))
    stop("'K' outside allowable range")
  if (type == "random") {
    R <- round(rep(R, length.out = 1))
    if (!isTRUE(R > 0))
      R <- 1
    subsets <- replicate(R, sample(n))
  }
  else {
    R <- 1
    subsets <- as.matrix(seq_len(n))
  }
  which <- rep(seq_len(K), length.out = n)
  if (type == "consecutive")
    which <- rep.int(seq_len(K), tabulate(which))
  folds <- list(n = n, K = K, R = R, subsets = subsets, which = which)
  class(folds) <- "cvFolds"
  folds
}

log_likelihood <- function(precision, TEST) {

  p <- nrow(precision)
  n <- nrow(TEST)
  emp_cov <- cov(TEST)
  logdet <- determinant(precision, logarithm=T)$modulus
  loglik <- n * (logdet - sum(diag(emp_cov %*% precision)))
  #loglik <- logdet - sum(diag(emp_cov %*% precision))
  return(as.numeric(loglik))
}
CV_score <- function(precision, TEST) {
  p <- nrow(precision)
  n <- nrow(TEST)
  emp_cov <- cov(TEST)
  logdet <- determinant(precision, logarithm=T)$modulus

  loglik <- 0
  for(k in 1:nrow(TEST)){
    loglik <- loglik + (t(TEST[k,]) %*% precision %*% (TEST[k,]))
  }
  loglik <- n*logdet -loglik
  return(as.numeric(loglik))
}

glasso_cv <- function(ts, cor = FALSE, nfolds=5, rholist=NULL,nlambda = 30 ,lambda.min.ratio =0.1,verbose=T) {

  if(!cor) S <- cov(ts)
  if(cor) S <- cor(ts)
  p <-ncol(ts)
  n <-nrow(ts)
  if (is.null(rholist)) {
    lambda.max = max(max(S - diag(p)), -min(S - diag(p)))
    lambda.min = lambda.min.ratio * lambda.max
    rholist = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  }
  if (!is.null(rholist)) nlambda = length(rholist)
  folds <- cvFolds(n, nfolds, type="random")

  CVscore <- sapply(1:nfolds, function(ki) {
    if (verbose) cat("Fold ", ki, "\n")
    mat_train =ts[folds$which!=ki,]
    mat_test = ts[folds$which==ki,]

    if(!cor) {
      S_train <- cov(mat_train)
      S_test  <- cov(mat_test)
    }
    if(cor) {
      S_train <- cor(mat_train)
      S_test  <- cor(mat_test)
    }

    curr_n <-nrow(mat_train)
    p <- ncol(mat_train)

    if(checkSPD(S_train)) S <- S_train
    if(!checkSPD(S_train)) S<- nearest_spd(S_train)

    GLP   <- huge(S, rholist , method = "glasso", scr=F)

    #loglike <- lapply(GLP$icov,function(P_train) log_likelihood(P_train, mat_test))
    CV_k <- lapply(GLP$icov,function(P_train) CV_score(P_train, mat_test))
    CV_k
  })

  CV.score = matrix(unlist(CVscore),nrow=nlambda,ncol=nfolds)
  rowSum = apply(CV.score,1, sum)

  ind     <- which.max(rowSum)
  rhomax  <- rholist[ind]
  a       <- list()
  a$rholist <- rholist
  a$optLambda <- rhomax
  a$maxCVscore <- rowSum[ind]
  a$CV_mat <- CV.score
  a$CVvalue <- rowSum
  return(a)
}


glasso_Select<- function(S,sample.size, rholist=NULL,nlambda = 30 ,gamma=0.5,lambda.min.ratio =0.1,penalty ="lasso",verbose=T) {

    p <- ncol(S)
    if (is.null(rholist)) {
      lambda.max = max(max(S - diag(p)), -min(S - diag(p)))
      lambda.min = lambda.min.ratio * lambda.max
      rholist = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
    }
    if (!is.null(rholist)) nlambda = length(rholist)
    GLP   <- huge(S, rholist, method = "glasso", scr=F)

    eBIC = rep(NA,length(rholist))
    BIC = rep(NA,length(rholist))
    AIC = rep(NA,length(rholist))
    LL_vec = rep(NA,length(rholist))
    for(i in 1:length(rholist)){
      precision = GLP$icov[[i]]
      k = GLP$df[i]  # or df in glasso
      LL = GLP$loglik[i]
      ebic_score = (-1*sample.size*LL) + (k*log(sample.size)) + (4*k*gamma*log(p))  ## gamma=0 for BIC
      bic_score = (-1*sample.size*LL) + (k*log(sample.size))  ## gamma=0 for BIC
      aic_score = (-1*sample.size*LL) + (2*k)

      LL_vec[i] = LL
      eBIC[i] = ebic_score
      BIC[i] = bic_score
      AIC[i] = aic_score
    }
    cid = which(is.infinite(LL_vec))
    cid1 = which(is.infinite(AIC))
    cid2 = which(is.infinite(BIC))
    cid3 = which(is.infinite(eBIC))
    if(length(cid)>0) LL_vec[cid] = NA
    if(length(cid1)>0) AIC[cid1] = NA
    if(length(cid2)>0) BIC[cid2] = NA
    if(length(cid3)>0) eBIC[cid3] = NA
    ret = list()
    ret$rholist = rholist
    ret$LL = LL_vec
    ret$eBIC = eBIC
    ret$AIC = AIC
    ret$BIC = BIC
    return (ret)
}



### function to convert a precision matrix to partial correlation ##
convertPC<-function(MAT){

  PartialCor = MAT
  for(i in 1:ncol(MAT)){
    for(j in 1:ncol(MAT)){
      pcor = -1* MAT[i,j]/(sqrt(MAT[i,i])*sqrt(MAT[j,j]))
      PartialCor[i,j] = pcor
    }
  }
  diag(PartialCor)=1
  PartialCor
}

createMatrixB<-function(MAT, Omega, Prior){
  pp = ncol(MAT)
  #samplePREC = solve(MAT)
  PartialCor = convertPC(Omega)
  matA = matrix(rep(0,pp*pp),nrow = pp,ncol=pp)

  for(i in 1:pp){
    cid_pos = which(PartialCor [,i]>0)
    cid_neg = which(PartialCor [,i]<0)
    if(length(cid_pos)>0) matA[cid_pos,i] = 1
    if(length(cid_neg)>0) matA[cid_neg,i] = -1
  }
  diag(matA) = 0

  unionAdj = unionBinMat(matA, Prior)

  matrixB = MAT
  for(i in 1:ncol(MAT)){
    matrixB[,i] = unionAdj[,i]*MAT[,i]
  }
  RET = list()
  RET$matrixA = matA
  RET$matrixB = matrixB
  RET
}

unionBinMat<-function(matA, Prior){
  ret = matA+Prior
  ret = sign(ret)
  ret
}


computeScore<-function(U,phi, W, matrixB){

  tmpS = matrixB %*% W %*% t(matrixB)
  ff = max(abs(tmpS))
  tmpS = abs(tmpS)/ff
  for(i in 1:ncol(tmpS)){
    cid_neg= which(phi[,i] == 1)
    if(length(cid_neg)>0) tmpS[cid_neg,i] = (tmpS[cid_neg,i]) * (-1)
    tmpS[,i] = tmpS[,i] * U[,i]
  }
  tmpS
}


Network2Adj<-function(MAT, dir, ff=NULL,  PartialCor,reportCor=FALSE){

  ## direction = 1 for 3 column network to adjacency matrix and feature list (ff) is needed
  ## direction = 2 for adjacency matrix to 3 column network

  if(dir==1){
    RET = matrix(0,nrow = length(ff), ncol = length(ff))
    row.names(RET) = colnames(RET) = ff
    allF = unique(as.character(MAT[,1]))
    for(i in 1:length(allF)){
      currF = allF[i]
      cid_i = which(ff==currF)
      tmp_N = MAT[which(MAT[,1] == currF),]
      cid_j = match(tmp_N[,2],ff)
      if(length(cid_i)!=0 & length(cid_j)!=0){
        RET[cid_i, cid_j] = tmp_N[,3]
        RET[cid_j, cid_i] = tmp_N[,3]
      }
    }
  }
  if(dir==2){
    RET = NULL
    if(!is.numeric(MAT[,1])) {
      MATnew = as.matrix(MAT[,-1])
      row.names(MATnew) = MAT[,1]
    }
    else {
      MATnew = MAT
      row.names(MATnew) = ff
    }
    for(i in 1:ncol(MATnew)){
      cid = which(MATnew[,i]!=0)
      cid = cid[which(cid>i)]
      if(length(cid)>0){
        nodeA = row.names(MATnew)[cid]
        nodeB = colnames(MATnew)[i]
        pp = PartialCor[cid,i]
        dir = sign(pp)
        if(reportCor){
          tmpD = data.frame(nodeA, nodeB, dir, pp)
          colnames(tmpD) = c("nodeA","nodeB","dir","partialcor")
        }
        if(!reportCor){
          tmpD = data.frame(nodeA, nodeB, dir)
          colnames(tmpD) = c("nodeA","nodeB","dir")
        }
      }

      RET = rbind(RET, tmpD)
    }
  }

  RET
}


getsigns<- function(matrixA,matrixP){
  ret = matrix(0, nrow= nrow(matrixA), ncol= ncol(matrixA))
  for(i in 1:ncol(matrixA)){
    cidA = which(matrixA[,i]==0 & matrixP[,i]!=0)
    if(length(cidA)>0) ret[cidA,i] = 1
  }
  ret
}
convertBinary<- function(mat){
  nr = nrow(mat)
  nc = nrow(mat)
  BIN = matrix(rep(0, nr*nc), nrow = nr, ncol = nc)
  for(i in 1:nr) for(j in 1:nc) if(mat[i,j]!=0) BIN[i,j] = 1
  row.names(BIN) = row.names(mat)
  colnames(BIN) = colnames(mat)
  BIN
}

cutFusedmat<-function(Fusedmat, cutoff){
  newMat = Fusedmat
  for(i in 1:ncol(newMat)){
    cid = which(abs(Fusedmat[,i])<cutoff)
    if(length(cid)>0) newMat[cid,i] = 0
  }
  newMat
}

AssessAccuracy<-function(TRUEadj, Predadj){
  TP = 0
  TN = 0
  FP = 0
  FN = 0
  for(i in 1:(nrow(Predadj)-1)){
    currV = Predadj[i,c((i+1):ncol(Predadj))]
    cid1 = which(currV==1)
    cid0 = which(currV==0)
    if(length(cid1)>0){
      checksum = (TRUEadj[i,(cid1+i)]==1)
      right = sum(checksum)
      wrong = length(checksum) - right
      TP = TP +right
      FP = FP + wrong
    }
    if(length(cid0)>0){
      checksum = (TRUEadj[i,(cid0+i)]==0)
      right = sum(checksum)
      wrong = length(checksum) - right
      TN = TN +right
      FN = FN + wrong
    }
  }

  SEN = TP/(TP+FN)
  SPEC = TN/(TN+FP)
  Fallout = FP/(FP+TN)
  Precision = TP/(TP+FP)
  MCC = ((TP*TN) - (FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  RET = list()
  RET$Sensitivity = SEN
  RET$Specificity = SPEC
  RET$Fallout = Fallout
  RET$Precision = Precision
  RET$MCC = MCC
  cat("Sensitivity = ",round(SEN,3),"\nSpecificity = ",round(SPEC,3),"\nFall-out = ",round(Fallout,3),"\nPrecision = ",round(Precision,3),"\nMCC = ",round(MCC,3) ,"\n")

  RET
}

createCoexpressionScore<-function(ff, M){

  RET = matrix(NA,nrow=nrow(ff), ncol=ncol(M))
  for(i in 1:nrow(ff)){
    sign = ff$Sign[i]
    gA = ff$GeneA[i]
    gB = ff$GeneB[i]
    cid_A = which(row.names(M)==gA)
    cid_B = which(row.names(M)==gB)
    if(length(cid_A)==0|length(cid_B)==0) ss = rep(NA,ncol(RET))
    else ss = M[cid_A,] + sign*M[cid_B,]
    RET[i,] = ss
  }
  row.names(RET) = paste0(ff$GeneA,".....",ff$GeneB)
  colnames(RET) = colnames(M)
  RET
}

checkMax<-function(x, grp){
  nlen = length(x)
  pos.max = which.max(x)
  max.val = x[pos.max]
  ret= NA
  if(sum(x==max.val,na.rm=T)>1) ret = "tied"
  if(sum(x==max.val,na.rm=T)==1) ret = grp[pos.max]
  ret
}


checkOverlaps<- function(NETA,NETB){
  strA = strB = NULL
  for(i in 1:nrow(NETA)){
    tmpstr = NETA[i,c(1:2)]
    tmpstr =tmpstr[order(tmpstr)]
    strA = c(strA, paste0(tmpstr,collapse = "_"))
  }
  for(j in 1:nrow(NETA)){
    tmpstr = NETB[j,c(1:2)]
    tmpstr =tmpstr[order(tmpstr)]
    strB = c(strB,  paste0(tmpstr,collapse = "_"))
  }
  len = length(intersect(strA,strB))
  cat(paste0("There were ", len," overlapping network edges between the estimated glasso network and the fused network."),"\n")
}
