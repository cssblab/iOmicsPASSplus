## ---- echo=FALSE-------------------------------------------------------------------------------------------------------

htmltools::img(src = knitr::image_uri("C:/Users/mdckwlh/Desktop/iOmicsPASSv2plus_repo/vignettes/figures/iOmicsPASSv2plus_logo_small.png"), alt = 'logo', style="position:absolute; top:0; right:0; padding:5px;height:125px;width:200px")

## ----setup, include=FALSE----------------------------------------------------------------------------------------------

knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>",
  rmarkdown.html_vignette.check_title = FALSE
)


## ----fig.align="center", echo=FALSE, out.width = '90%'-----------------------------------------------------------------
knitr::include_graphics("figures/Overview_iOmicsPASSplus.png")

## ----eval=FALSE--------------------------------------------------------------------------------------------------------
#  if (!require(devtools)) install.packages("devtools")

## ----fig.align="center", echo=FALSE, out.width = '100%'----------------------------------------------------------------
knitr::include_graphics("figures/Cygwin_screenshot1.png")

## ----fig.align="center", echo=FALSE, out.width = '100%'----------------------------------------------------------------
knitr::include_graphics("figures/Cygwin_screenshot2.png")
knitr::include_graphics("figures/Cygwin_screenshot3.png")

## ----eval=FALSE--------------------------------------------------------------------------------------------------------
#  setwd("C:/PATH_TO_PROGRAM/iOmicsPASSplus/")
#  devtools::install_local("iOmicsPASSplus.tar.gz", dependencies =T)
#  ## Alternatively ##
#  devtools::install_github("CSSBlab/iOmicsPASSplus", build_vignettes = TRUE)
#  
#  # load the library
#  library(iOmicsPASSplus)
#  ## only need to be run once to compile iOmicsPASS, creating a program in /bin folder.
#  INSTALL.iOmicsPASS()
#  

## ----echo=FALSE--------------------------------------------------------------------------------------------------------
library(iOmicsPASSplus)

## ----------------------------------------------------------------------------------------------------------------------
## load the example data ##
data(Tulip_Protein)
data(Tulip_microRNA)
data(PhenotypeFile)

head(Tulip_Protein[,c(1:6)])

head(Tulip_microRNA[,c(1:6)])

head(PhenotypeFile)


## ----------------------------------------------------------------------------------------------------------------------

data(PPI_network)
data(TargetScan_network)

head(PPI_network)
head(TargetScan_network)

## ----------------------------------------------------------------------------------------------------------------------
data(bioPathways)
head(bioPathways)

## ---- eval=F-----------------------------------------------------------------------------------------------------------
#  ## creating an list object containing the two datasets and labeling them accordingly
#  row.names(Tulip_Protein) = Tulip_Protein$Protein
#  row.names(Tulip_microRNA) = Tulip_microRNA$miRNA
#  Tulip_Protein = Tulip_Protein[,-1]
#  Tulip_microRNA = Tulip_microRNA[,-1]
#  inputDat=list(Tulip_Protein, Tulip_microRNA)
#  names(inputDat) = c("Protein","microRNA")
#  
#  NetDeconvolute(inputDat, option=1,log.transform=TRUE, tag="supervised",criterion="eBIC",Calibration=TRUE, verbose=T)

## ----fig.align="center", echo=FALSE, out.width = '60%'-----------------------------------------------------------------
knitr::include_graphics("figures/DataQCplots.png")

## ----fig.align="center", echo=FALSE, out.width = '50%'-----------------------------------------------------------------
knitr::include_graphics("figures/Heatmap_CrossCovarianceMatrix_supervised2.png")

## ----fig.align="center", echo=FALSE, out.width = '70%'-----------------------------------------------------------------
knitr::include_graphics("figures/CalibrationPlots_glasso_supervised.png")


## ----eval=F------------------------------------------------------------------------------------------------------------
#  ## Refining a narrower lambda vector ##
#  lambda_new=exp(seq(log(0.5),log(0.01), length=30))
#  NetDeconvolute(inputDat, option=1,log.transform=TRUE, tag="supervised2",criterion="eBIC", Calibration=TRUE, lambda.vec=lambda_new, verbose=T)

## ----fig.align="center", echo=FALSE, out.width = '70%'-----------------------------------------------------------------
knitr::include_graphics("figures/CalibrationPlots_glasso_supervised2.png")

## ----eval=F------------------------------------------------------------------------------------------------------------
#  NetDeconvolute(inputDat, option=1,log.transform=TRUE,tag="supervised",criterion="eBIC",
#  Calibration=FALSE, optLambda=0.382,verbose=TRUE)

## ----fig.align="center", echo=FALSE, out.width = '70%'-----------------------------------------------------------------
knitr::include_graphics("figures/Plots_glasso_supervised.png")

## ---- eval=T-----------------------------------------------------------------------------------------------------------
data(TargetScan_network)
data(PPI_network)

head(TargetScan_network)
head(PPI_network)

## ---- eval=T-----------------------------------------------------------------------------------------------------------
TargetScan_network$NodeA_DT = "X"
TargetScan_network$NodeB_DT = "Y"
PPI_network$NodeA_DT = "X"
PPI_network$NodeB_DT = "X"
colnames(TargetScan_network) = c("NodeA","NodeB","Dir","NodeA_DT","NodeB_DT")
colnames(PPI_network) = c("NodeA","NodeB","Dir","NodeA_DT","NodeB_DT")
PriorNet = rbind(TargetScan_network,PPI_network)

head(PriorNet)

## ---- eval=F-----------------------------------------------------------------------------------------------------------
#  NetDeconvolute(inputDat, option=2, NetworkFile=PriorNet, tag="hybrid",criterion="eBIC",
#  log.transform=TRUE,Calibration=FALSE, optLambda=0.382, verbose=TRUE)

## ----fig.align="center", echo=FALSE, out.width = '85%'-----------------------------------------------------------------
knitr::include_graphics("figures/Plots_Hybridmethod_hybrid.png")

## ---- eval=F-----------------------------------------------------------------------------------------------------------
#  ## Using known biological networks ##
#  createInputParam(data.X=Tulip_Protein,data.Y=Tulip_microRNA,phenotype=PhenotypeFile,log.transform=TRUE,btw.net=TargetScan_network,within.net=PPI_network, Cross.Validate = TRUE, Enrichment=FALSE, tag="KnownNetwork")
#  
#  ## run iOmicsPASS to create the misclassification error plot for picking threshold
#  ## by Default, all results will be output to /iOmicsPASS/Output/.
#  iOmicsPASS.output<-iOmicsPASS.R(ff="input_param_KnownNetwork")

## ----fig.align="center", echo=FALSE, out.width = '70%'-----------------------------------------------------------------
knitr::include_graphics("figures/CVplot_Penalty_KnownNetwork.png")

## ---- eval=F-----------------------------------------------------------------------------------------------------------
#  ## rerun with threshold = 3.53 ##
#  createInputParam(data.X=Tulip_Protein,data.Y=Tulip_microRNA,phenotype =PhenotypeFile,
#  log.transform=TRUE,btw.net=TargetScan_network,within.net=PPI_network, Enrichment=FALSE, tag="KnownNetwork", min.thres=3.53,Cross.Validate=FALSE)
#  
#  ## run iOmicsPASS and set plotCV=FALSE ##
#  iOmicsPASS.output<-iOmicsPASS.R(ff="input_param_KnownNetwork", plotCV=FALSE, Cross.Validate = FALSE)
#  
#  ## Or, you can also carry out pathway enrichment on the selected proteins by setting Enrichment=TRUE ##
#  createInputParam(data.X=Tulip_Protein,data.Y=Tulip_microRNA,phenotype =PhenotypeFile, log.transform=TRUE,btw.net=TargetScan_network,within.net=PPI_network, pathway = bioPathways, Enrichment=TRUE, tag="withEnrichment", min.thres=3.53,Cross.Validate=FALSE)
#  
#  iOmicsPASS.output<-iOmicsPASS.R(ff="input_param_withEnrichment", plotCV=FALSE, Cross.Validate = FALSE)

## ----fig.align="center", echo=FALSE, out.width = '90%'-----------------------------------------------------------------
knitr::include_graphics("figures/Network_OIR.png")

## ---- eval=F-----------------------------------------------------------------------------------------------------------
#  ## Running with estimated glasso network from supervised approach ##
#  createInputParam(data.X="Combined_data.txt", within.net="Fusednetwork_hybrid.txt",
#  phenotype =PhenotypeFile,log.transform=FALSE, standardize.data = FALSE,Cross.Validate = TRUE,Enrichment=FALSE, tag="hybridNetwork")
#  
#  iOmicsPASS.supervisednet <-iOmicsPASS.R(ff="input_param_hybridNetwork")

## ----fig.align="center", echo=FALSE, out.width = '70%'-----------------------------------------------------------------
knitr::include_graphics("figures/CVplot_Penalty_Hybrid.png")

## ---- eval=F-----------------------------------------------------------------------------------------------------------
#  createInputParam(data.X="Combined_data.txt", within.net="Fusednetwork_hybrid.txt",
#  phenotype =PhenotypeFile,log.transform=FALSE, standardize.data = FALSE,min.thres=3.3, Cross.Validate = FALSE,pathway=bioPathways,Enrichment=TRUE, tag="hybridNetwork")
#  
#  iOmicsPASS.supervisednet <-iOmicsPASS.R(ff="input_param_hybridNetwork",Cross.Validate = F, plotCV = F)
#  

## ---- eval=T-----------------------------------------------------------------------------------------------------------
data(PhenotypeFile)
row.names(PhenotypeFile) = PhenotypeFile$TulipID
PhenotypeFile=PhenotypeFile[,-1]
head(PhenotypeFile)

## ---- eval=F-----------------------------------------------------------------------------------------------------------
#  priorProb = createPrior(PhenotypeFile, y = "Group",outputDir = "iOmicsPASS/inputFiles/")
#  
#  createInputParam(data.X=Tulip_Protein,data.Y=Tulip_microRNA,phenotype =PhenotypeFile, log.transform=TRUE,btw.net=TargetScan_network,within.net=PPI_network,tag="withPrior", Cross.Validate=T, usePrior = T, priorfile = priorProb)
#  
#  iOmicsPASS.output<-iOmicsPASS.R(ff="input_param_withPrior", plotCV=FALSE, Cross.Validate = FALSE)
#  

