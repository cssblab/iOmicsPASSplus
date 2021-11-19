#' Phenotype file for the Tulip study
#'
#' A data frame with two columns containing the phenotype assigned to the 18 subjects where
#' 8 were obese insulin resistant (OIR, HOMA-IR >2.5) and 9 were lean insulin-sensitive
#' (LIS, HOMA-IR <1.0) normoglycemic males.
#'
#' @docType data
#'
#' @usage data(PhenotypeFile)
#'
#' @format A two-column data frame with 17 subjects (row) and their phenotype:
#' \describe{
#'   \item{TulipID}{Subject identifier in the Tuplip study}
#'   \item{Group}{Phenotypic group of each subject. OIR for obese insulin resistant and LIS for lean insulin sensitive.}
#'   \item{Age}{Age of the subject}
#'   \item{BMI}{Body mass index of each subject}
#' }
#' For further details, see \url{https://pubmed.ncbi.nlm.nih.gov/31024340/}
#'
"PhenotypeFile"

#' Biological pathways collated from multiple databases.
#'
#' A data frame with information of the biological pathways that each gene symbol
#' is associated with, collated from multiple databases including Gene Ontology (GO)
#' Consortium, KEGG, Reactome Pathway database and ConsensusPathDB.
#'
#' @docType data
#'
#' @usage data(bioPathways)
#'
#' @format A data frame with 17,250 unique gene symbols and 14,598 unique biological pathways:
#' \describe{
#'   \item{Genesym}{Gene symbol}
#'   \item{Pathwayid}{Pathway Identifier}
#'   \item{Function}{Description of the pathway}
#' }
#' For further details, see \url{https://pubmed.ncbi.nlm.nih.gov/31024340/}
#'
"bioPathways"

#' Protein-Protein interaction network file
#'
#' A data frame with three columns where the first two columns contain the pairs of protein identifiers
#' that chemically/physically interact with each other. The third column indicates the sign of the direction
#' of that interaction, which are assumed to be all positive interactions here.
#'
#' @docType data
#'
#' @usage data(PPI_network)
#'
#' @format A data frame with 1,201 unique proteins containing 1,499 interactions:
#' \describe{
#'   \item{geneA_genesym}{Gene symbol of Protein A}
#'   \item{geneB_genesym}{Gene symbol of Protein B that is interacting with Protein A}
#'   \item{sign}{Direction of protein-protein interaction, assumed to be all positive (indicated by 1)}
#' }
#' For further details, see \url{https://pubmed.ncbi.nlm.nih.gov/31024340/}
#'
"PPI_network"

#' microRNA to target gene network file
#'
#' A dataframe with three columns where the first two columns represent the miRNA probes and their genes targets.
#' The third column indicates the direction of interaction, where microRNAs that are translation
#' inhibitor of its target genes as negative and those that regulates gene translation as positive.
#'
#' @docType data
#'
#' @usage data(TargetScan_network)
#'
#' @format A data frame with 882 unique genes and 125 unique microRNAs with 8,533 microRNA-gene targets:
#' \describe{
#'   \item{GeneSym}{Gene symbol}
#'   \item{miRNA}{microRNA probe ID}
#'   \item{sign}{Direction of miRNA-gene regulation where positive regulation of gene translation are indicated by 1 and inhibition of translation by -1.}
#' }
#' For further details, see \url{https://pubmed.ncbi.nlm.nih.gov/31024340/}
#'
"TargetScan_network"


#' Protein expression data in Tulip study. Original data contains 1,499 proteins quantified and only
#' those that were different (p<0.1) between insulin resistant and insulin sensitive individuals
#' using 2-sample t-test were included in this example data.
#'
#' A dataset containing the 266 plasma protein abundance quantified by LC-MS across
#' 17 males.
#'
#' @docType data
#'
#' @usage data(Tulip_Protein)
#'
#' @format A data frame with 1,499 proteins (rows) across 17 subjects (columns):
#' \describe{
#'   \item{Protein}{Protein symbol}
#'   \item{TulipXX}{Sample ID in Tulip study}
#' }
#' For further details, see \url{https://pubmed.ncbi.nlm.nih.gov/31024340/}
#'
"Tulip_Protein"

#' microRNA data in Tulip study.
#'
#' A dataset containing 263 normalized microRNA copy number using multiplex RT-qPCR platform, MiRXES.
#' Original data contains 368 microRNA probes and only those that were different (p<0.1) between insulin
#' resistant and insulin sensitive individuals using 2-sample t-test were included in this example data.
#' @docType data
#'
#' @usage data(Tulip_microRNA)
#'
#' @format A data frame with 263 microRNA probes (rows) across 17 subjects (columns):
#' \describe{
#'   \item{miRNA}{microRNA probes}
#'   \item{TulipXX}{Sample ID in Tulip study}
#' }
#' For further details, see \url{https://pubmed.ncbi.nlm.nih.gov/31024340/}
#'
"Tulip_microRNA"



