# iOmicsPASS+
## Integrative -Omics Predictive Analysis of Subnetwork Signatures (Version II - *An R-package*)

To start, either download `iOmicsPASSplus.zip` file and unzip to local directory or use command line/Terminal to clone the entire github directory:
```{bash eval=FALSE}
> git clone https://github.com/cssblab/iOmicsPASSplus.git
```

### Introduction

**iOmicsPASS+** is a R-package incorporating **iOmicsPASS** (Koh et al., 2019), extended to other types of -omics data allowing for flexibility and increasing usability. It includes several module including a network inference module `NetDeconvolute()` using graphical LASSO (glasso) to estimate a sparse inverse covariance matrix, creating a confounding-free partial correlation network among features from up to three -omics datasets.

**iOmicsPASS** has been improved to **iOmicsPASS+** allowing for higher flexibility and enabling applications to different types of omics data. Improvements include:

* **Specification of direction of association**\
  Users may now specify the direction for every pair of interacting or co-varying molecule by adding an additional column in the network file. However, only molecules that show consistent sign of correlation in the empirical data as the user-specified direction will be considered.
  
* **Allows for a single network and input data**\
  Previously, at least two data and two networks were required as input. Now, users can input only one single data and create co-expressions among the variables in the data with a single network file.
  
* **Addition of a Network estimation module** `NetDeconvolute()`\
  Estimates a correlation network, linking the different features from up to three different data, using graphical LASSO (glasso) to estimate a sparse inverse covariance matrix, creating a confounding-free partial correlation network
  
* **New functions to help users compile and run iOmicsPASS using R**\
  Functions included in the R package facilitate users to build `INSTALL.iOmicsPASS()`, create input parameter file `createInputParam()`, create prior probabilities `createPrior()` and run the software `run.iOmicsPASS()` in the R-console.
  
* **Addition of a Prediction module** `Predict.iOmicsPASS()`\
  Uses the network signatures identified in the subentwork discovery module `run.iOmicsPASS()`to assign new samples to the phenotypic groups.
  
* **Adjustment for clinical information**\
  Users can incorporate clinical information such as age, gender and BMI, to modify the prior class probabilities used for assigning samples to the different groups.
  
  The figure below illustrates the overview of **iOmicsPASS+**
![Overview_iOmicsPASSplus](https://user-images.githubusercontent.com/37172948/142598004-93b4ed3a-42b9-428d-8620-92267601e840.png)
