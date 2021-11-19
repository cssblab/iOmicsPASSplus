\name{iOmicsPASSv2plus-package}
\alias{iOmicsPASSv2plus-package}
\alias{iOmicsPASSv2plus}
\docType{package}
\title{
Integrative -Omics for Prediction Analysis of Subnetwork Signatures -- version 2 extended
}
\description{
\code{iOmicsPASSv2plus} is a R-package incorporating iOmicsPASS (Koh et al., 2019), extended to other types of -omics data allowing for flexibility and increasing usability. It includes several module including a network inference module, \code{NetDeconvolute()} using graphical LASSO to estimate a sparse inverse covariance matrix, creating a confounding-free partial correlation network among features from up to three -omics datasets. The estimated network can be used to create co-expression scores in \code{iOmicsPASS.R()} to identify predictive signatures that best separates phenotypic outcomes. Those signatures can also be used to assign new samples with the same -omics data available into the phenotypic outcomes using \code{Predict.iOmicsPASS()}.
}

\author{
 Hiromi WL Koh, Hyungwon Choi
}
\note{
User needs to have a GNU compiler to compile the C++ program. Detailed steps can be found in the vignette.
}

\keyword{package}
\keyword{iOmicsPASS}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
