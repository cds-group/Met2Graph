% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Met2Graph.R
\name{Met2RxnGraph}
\alias{Met2RxnGraph}
\title{Met2RxnGraph.R}
\usage{
Met2RxnGraph(
  infile,
  rmMets,
  recMets_list,
  outDir,
  outFormat,
  directed,
  add_rev_rxn
)
}
\arguments{
\item{infile}{Path to the input SBML file}

\item{rmMets}{logical value indicating whether the recurring metabolites must be removed. Default TRUE.}

\item{recMets_list}{a list of recurring metabolite names. The format must be the same of the model.
If not provided and rmMets=TRUE the list used is the one provided by the package "system.file("extdata", "metNames.txt", package = "Met2Graph")".}

\item{outDir}{path to output directory. If not provided outputs will be written in "system.file("extdata/output/", package = "Met2Graph")".}

\item{outFormat}{character string indicating the format of the graph to export. Default "ncol". For other formats and details please see write.graph function of igraph package.}

\item{directed}{Logical value, indicating if the graph must be created directed. Default TRUE.}

\item{add_rev_rxn}{Logical value, whether or not reversible reactions must be considered and edges written in both directions. If FALSE it takes into account the direction as present in the stoichiometric matrix, for general models it is usually forced from left to right.
In case of condition-specific models obtained through FBA this is not necessary. Default FALSE.}
}
\value{
A metabolic graph in tabular format ("from","to","Metabolite") and graph format as chosen by the user.
}
\description{
This function loads a metabolic model in SBML format and returns a metabolic graph having reactions as nodes and edges connecting two if one produces a metabolite consumed by the other one.
}
\details{
Builds a graph from the metabolic model.
The reactions are connected by metabolites.
Two reactions are connected if one produces a metabolite, which is a substrate of the other one.
These graphs are unweighted.
Recurring metabolites (e.g. ATP, CO2, NADP) can be removed to avoid unrealistic paths setting the parameter "rmMets"==TRUE and using the list of these compounds contained in the data of the package by default or one provided by the user through parameter "recMets".
If the parameter add_rev_rxn is set to TRUE, the edges will be assigned considering both directions of reversible reactions.
The resulting graph can be built directed or not. Given the rules behind the definition of connections it is naturally directed but the user can choose to not consider this aspect, setting the "directed" parameter to FALSE.
}
\examples{
\dontrun{
infile <- system.file("extdata", "kidney.xml", package = "Met2Graph")
outDir <- system.file("extdata/output/", package = "Met2Graph")
Met2RxnGraph(infile, rmMets = TRUE, outDir=outDir, outFormat = "graphml")}
}
