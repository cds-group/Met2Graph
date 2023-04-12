###############################################################################
#
# Met2Graph.R:     This file contains all functions related to the creation of graphs from metabolic models in SBML format.
# author: Ilaria Granata <ilaria.granata@icar.cnr.it>
#
#
# Documentation was created using roxygen
#
###############################################################################

#' MetExtract.R
#'
#'
#' This function loads a metabolic model in SBML format and returns a list object containing the model and data frames of reactants and products for each reaction.
#'
#' @param infile Path to the input SBML file
#' @import sybilSBML
#' @importFrom stats na.omit
#' @importFrom tools file_path_sans_ext
#' @return model and data frames of reactants and products for each reaction. The output values are stored in a list.
#' @export
#'
#' @examples
#' \dontrun{
#'out <- MetExtract(system.file("extdata", "kidney.xml", package = "Met2Graph"))}
MetExtract <- function(infile){
  #read SBML file and return stoichiometric matrix
  mod <- suppressWarnings(sybilSBML::readSBMLmod(infile, bndCond = FALSE))
  mod_name <- file_path_sans_ext(basename(infile))
  stoich_mat <- as.data.frame(as.matrix(mod@S))
  colnames(stoich_mat) <- mod@react_id
  rownames(stoich_mat) <- mod@met_id
  #from stoich matrix extract positive relationships between rxn and met
  pos_ind <- which(stoich_mat > 0, arr.ind=TRUE)
  pos_ind[,"row"] <- rownames(as.data.frame(stoich_mat))[as.numeric(pos_ind[,"row"])]
  pos_ind[,"col"] <- colnames(as.data.frame(stoich_mat))[as.numeric(pos_ind[,"col"])]
  #from stoich matrix extract negative relationships between rxn and met
  neg_ind <- which(stoich_mat < 0, arr.ind=TRUE)
  neg_ind[,"col"] <- colnames(as.data.frame(stoich_mat))[as.numeric(neg_ind[,"col"])]
  neg_ind[,"row"] <- rownames(as.data.frame(stoich_mat))[as.numeric(neg_ind[,"row"])]
  colnames(pos_ind)=c("product","Rxn")
  rownames(pos_ind)=NULL
  colnames(neg_ind)=c("reactant","Rxn")
  rownames(neg_ind)=NULL
  product_ind <- na.omit(pos_ind)
  reactant_ind <- na.omit(neg_ind)
  #output
  out <- list(mod, mod_name,product_ind,reactant_ind)              # Store output in list
  return(out)
}


#' Met2MetGraph.R
#'
#'
#' This function loads a metabolic model and returns a metabolic graph having metabolites as nodes and edges connecting reagents to products.
#'


#' Builds a metabolites-based graph from the metabolic model.
#' The metabolites are connected by reactions and, if present, by the relative catalyzing enzymes.
#' Two metabolites are connected if, in one or multiple reactions, one is a substrate and the other one is a product.
#' The gene-protein-reaction (GPR) associations are automatically extracted from the model. If not present, those extracted from the human generic GEM are used.
#'
#' If the parameter catalyzed is TRUE, only reactions with associated enzymes are kept.
#'
#' If gene expression data are provided each edge represented by the enzymes catalyzing the reactions is weighted by the expression values.
#'
#' The resulting graph can be obtained with multiple edges corresponding to the multiple enzymes connecting the two metabolites or can be simplified to single edges, the choice is indicated by the parameter "simpl", where a logical
#' value is required.
#'
#' Several methods of simplification are proposed, and can be indicated through the parameter "GPRparse".
#' For this parameter three options are possible:
#' 1) "meanSum": for each edge expression values of enzymes acting in the same reaction are averaged, and these averages are then summed up for enzymes acting in different reactions;
#' 2) "minMax": according to the GPR associations, enzymes catalyzing the same reaction are associated by "AND" or "OR" logical relationships based on complex or mutual exclusive activity respectively.
#' For each edge the weight is given by the minimum expression value of enzymes related by "AND" and the maximum value for "OR", again the values of enzymes acting in different reactions are summed up;
#' 3) "minSum": according to the GPR associations, enzymes catalyzing the same reaction are associated by "AND" or "OR" logical relationships based on complex or mutual exclusive activity respectively.
#' For each edge the weight is given by the minimum expression value of enzymes related by "AND" and the summed value for "OR", again the values of enzymes acting in different reactions are summed up.
#'
#' If the parameter add_rev_rxn is set to TRUE, each edge will be duplicated in both directions in case of reversible reactions.
#'
#' Recurring metabolites (e.g. ATP, CO2, NADP) can be removed to avoid unrealistic paths setting the parameter "rmMets"==TRUE and using the list of these compounds contained in the data of the package by default or one provided by the user through parameter "recMets".
#'
#' The resulting graph can be built directed or not. Given the rules behind the definition of connections it is naturally directed but the user can choose to not consider this aspect, setting the "directed" parameter to FALSE.
#'
#'
#'
#' @param infile Path to the input SBML file
#' @param catalyzed logical value indicating whether only reactions catalyzed by enzymes must be considered. Default TRUE.
#' @param rmMets logical value indicating whether the recurring metabolites must be removed. Default TRUE.
#' @param recMets_list a list of recurring metabolite names. The format must be the same of the model. If not provided and rmMets=TRUE the list used is the one provided by the package "system.file("extdata", "metNames.txt", package = "Met2Graph")".
#' @param exprDir path to directory containing expression data for each sample. Test data in "system.file("extdata/expression/", package = "Met2Graph")". Expression files must contain two columns, one for gene id (corresponding to the ids contained into the metabolic model, usually ENSEMBL ids) and one for expression values.
#' If the ENSEMBL id contains the version, this latter will be ignored.
#' @param delim column separator of expression files. Default tab.
#' @param simpl logical value indicating whether the simplification of multiple edges must be applied. If TRUE, GPRparse parameter needs to be set. Default TRUE.
#' @param GPRparse Character string indicating the method for GPR parsing. meanSum, minMax, minSum are implemented. Default meanSum.
#' @param outDir path to output directory. If not provided outputs will be written in "system.file("extdata/output/", package = "Met2Graph")".
#' @param outFormat character string indicating the format of the graph to export. Default "ncol". For other formats and details please see write.graph function of igraph package.
#' @param directed Logical value, indicating if the graph must be created directed. Default TRUE.
#' @param add_rev_rxn Logical value, whether or not reversible reactions must be considered and edges written in both directions. If FALSE it takes into account the direction as present in the stoichiometric matrix, for general models it is usually forced from left to right In case of condition-specific models obtained through FBA this is not necessary. Default FALSE.
#'
#' @importFrom Hmisc %nin%
#' @import stringr
#' @import stringi
#' @importFrom reshape2 melt
#' @importFrom tidyr separate replace_na
#' @importFrom igraph graph_from_data_frame write.graph simplify as_data_frame plot.igraph
#' @importFrom dplyr rename select inner_join
#' @importFrom tools file_path_sans_ext
#' @importFrom stats na.omit
#' @importFrom utils data read.table write.table
#' @return A metabolic graph with metabolites as nodes and reactions/enzymes as edges in the format chosen by the user. If expression values are present the edges are weighted by them, and one graph per each sample is returned. A tsv file with edges and expression values of each sample is also generated.
#' @export
#'
#' @examples
#' \dontrun{
#' infile <- system.file("extdata", "kidney.xml", package = "Met2Graph")
#' exprDir <- system.file("extdata/expression/", package = "Met2Graph")
#' outDir <- system.file("extdata/output/", package = "Met2Graph")
#' Met2MetGraph(infile, catalyzed = TRUE, rmMets = TRUE, exprDir=exprDir, simpl = TRUE , GPRparse = "meanSum", outDir=outDir, outFormat = "ncol", directed=TRUE, add_rev_rxn=FALSE, plot=FALSE)}
#'

Met2MetGraph <- function(infile, catalyzed, rmMets, recMets_list, exprDir, delim, simpl, GPRparse, outDir, outFormat, directed, add_rev_rxn, plot) {

  # check for necessary specs

  if(!exists("infile")) {
    print("ERROR: no file in input")
    stop()
  }

  if (missing(catalyzed)) {
    print ("WARNING: catalyzed is not specified, default TRUE")
    catalyzed=TRUE
  } else {  catalyzed=catalyzed
  }

  if(missing(rmMets)) {
    print("WARNING: rmMets is missing, default TRUE")
    rmMets= TRUE
  } else {
    rmMets=rmMets }

  if (missing(recMets_list) && rmMets==FALSE) {
    recMets=NULL
  }
  else if (missing(recMets_list) && rmMets==TRUE) {
    print ("using the list of recurring metabolite names provided by the package")
    recMets_list= read.delim(system.file("extdata", "metNames.txt", package = "Met2Graph"))
    recMets_list=recMets_list }
  else if (!missing(recMets_list) && rmMets==TRUE) {
    recMets_list=recMets_list
  }


  if(missing(exprDir)) {
    print("WARNING: exprDir is missing, the edges will not be weighted and one graph will be generated")
    exprDir=NULL
    simpl=NULL
    GPRparse=NULL
  } else {
    exprDir==exprDir
    if (missing(simpl)) {
      print ("WARNING: simpl parameter is missing, default FALSE")
      simpl=TRUE}
    else if (!missing(simpl)) {
      simpl=simpl
    }
    if(missing(GPRparse) && simpl==TRUE) {
      print ("WARNING: GPRparse parameter is missing, default meanSum")
      GPRparse="meanSum" }
    else if(missing(GPRparse) && simpl==FALSE) {
      GPRparse=NULL}
    else if (!missing(GPRparse)){
      GPRparse=GPRparse }
  }

  if(!is.null(exprDir) && missing(delim)) {
    print("WARNING: delim parameter for expression files is not indicated, default TAB")
    delim= " "}
  else if (!missing(exprDir) && !missing(delim)){
    delim=delim}
  else if (missing(exprDir)){
    delim=NULL}
  else if (missing(delim)){
    delim=NULL}


  if(missing(outDir)) {
    print ("WARNING: No outDir path provided. Output directory will be written in extdata/output")
    outDir=system.file("extdata/output/", package = "Met2Graph")
  } else {
    outDir=outDir
  }

  if(missing(outFormat)) {
    print("WARNING: outFormat is missing, default ncol")
    outFormat= "ncol"
  } else {
    outFormat=outFormat
  }

  if(missing(directed)) {
    print("WARNING: directed parameter is not indicated, default TRUE")
    directed= TRUE
  } else {
    directed=directed
  }
  if(missing(add_rev_rxn)) {
    print("WARNING: add_rev_rxn parameter is not indicated, default FALSE")
    add_rev_rxn=FALSE
  } else {
    add_rev_rxn=add_rev_rxn
  }

  if(missing(plot)) {
    print("WARNING: plot parameter is not indicated, default FALSE")
    plot=FALSE
  } else {
    plot=plot
  }

  #extract reactants and products
  print("Metabolites-based graph building...")
  my_mod<-MetExtract(infile)
  subdir<-paste(outDir,"MetGraphs",sep="/")
  dir.create(subdir)
  mod=my_mod[[1]]
  mod_name=my_mod[[2]]
  product_ind=as.data.frame(my_mod[[3]])
  reactant_ind=as.data.frame(my_mod[[4]])
  directed=directed
  delim=delim
  #extract GPR
  print ("automatic extraction of GPR from the model")
  gpr=data.frame(Rxn=mod@react_id,GPR=mod@gpr)
  if (all(sapply(gpr[,2], function(x)all(x==""))) == TRUE)  {
    print("gpr values are missing in the model, those extracted from human HMR2 will be used")
    gpr=read.delim(system.file("extdata", "gpr.txt", package = "Met2Graph"),sep="\t")
  }
  gpr[,2]=as.character(gpr[,2])
  gpr[,2]=gsub("[()]", "", gpr[,2])
  #merge reactants and products by reaction id
  merge_neg_pos=merge(reactant_ind,product_ind,by="Rxn")
  gpr_split=data.frame(gpr[,1],str_split_fixed(gpr[,2],"\\ and |\\ or ",n=Inf))
  gpr_split=as.data.frame(apply(gpr_split,2,function(x)gsub('\\s+', '',x)))
  inpGraph=merge(merge_neg_pos,gpr_split,by.x="Rxn",by.y=colnames(gpr_split)[1])
  inpGraph=suppressWarnings(melt(inpGraph, id.vars=c("Rxn", "reactant", "product")))
  inpGraph=select(inpGraph, -c("variable"))
  inpGraph=inpGraph[,c(2,3,1,4)]
  inpGraph=inpGraph %>%
    rename(from = reactant,
           to=product,
           Reaction = Rxn,
           gene=value)
  #reversible reactions
  if (add_rev_rxn==TRUE){
    inpGraph$Reversible=mod@react_rev[match(inpGraph$Reaction,mod@react_id)]
    inpGraph=rbind(inpGraph,
                   inpGraph %>%
                     filter(Reversible==TRUE) %>%
                     mutate(from=to,
                            to=from))
    inpGraph_unique=distinct(inpGraph, from, to, Reaction, gene, Reversible)
  }
  #catalyzed reactions
  if (catalyzed==TRUE){
    inpGraph=inpGraph[!(inpGraph$gene==""), ]
    inpGraph=na.omit(inpGraph)
  }
  #remove recurring metabolites
  if (rmMets==TRUE){
    print("removing recurring metabolites...")
    mets=data.frame(metsName=mod@met_name,metsID=mod@met_id)
    recMets=subset.data.frame(mets,tolower(mets$metsName) %in% tolower(recMets_list[,1]))
    print("the following recurring metabolites have been removed:")
    print(unique(recMets$metsName))
    inpGraph = inpGraph[which(inpGraph$from  %nin% recMets[,2]), ]
    inpGraph = inpGraph[which(inpGraph$to  %nin% recMets[,2]), ]
  }
  #write output
  if (is.null(exprDir)){
    fileNameout <- file.path(subdir, paste0(mod_name,"_metabolites_based_graph.tsv"))
    write.table(inpGraph, fileNameout, sep="\t")
    graphMet <- graph_from_data_frame(inpGraph, directed=directed)
    fileNameout <- file.path(subdir, paste0(mod_name,"_metabolites_based_graph", '.',outFormat))
    write.graph(graphMet, fileNameout, outFormat)
  }

  if (!is.null(exprDir)) {
    print("adding expression data...")
    # read expression data:
    expr_files <- list.files(path=exprDir,
                             pattern = "*", full.names = T, recursive = FALSE)
    expr_files_trim <- list.files(path=exprDir,
                                  pattern = "*", full.names = F, recursive = FALSE)
    basenames<-file_path_sans_ext(basename(expr_files_trim))

    fileNumbers <- seq(expr_files_trim)


    graph_list=list()
    for (x in fileNumbers) {
      sample <- read.delim(expr_files[x],
                           header = FALSE, row.names = NULL,sep=delim)
      colnames(sample)[1]="gene"
      colnames(sample)[2]="expr"
      #ignore version of ENSEMBL ids if present
      sample=suppressWarnings(sample %>%
                                separate(gene, c("gene", "version"), "\\."))

      inpGraph$expr=sample$expr[match(inpGraph$gene,sample$gene)]
      inpGraph$expr[inpGraph$expr==0] <- 0.00000001
      inpGraph=inpGraph %>% replace_na(list(expr = 0))

      #multiple-edges graphs
      if (simpl==FALSE) {
        graph_list[[x]] <-inpGraph
        graph_df <- suppressWarnings(Reduce(function(...) merge(..., by = c('from', 'to','Reaction','gene')), graph_list))
        suppressWarnings(names(graph_df)[-(1:4)]<-basenames)
        write.table(graph_df,
                    paste(subdir,"Weights_df.tsv",sep="/"),
                    append = FALSE,
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)
        newFileName <-  paste0("expr_",basenames[x])
        graphMet <- graph_from_data_frame(inpGraph, directed = directed)
        fileNameout <- file.path(subdir, paste0(newFileName,".", outFormat))
        write.graph(graphMet, fileNameout, outFormat)
      }


      #meanSum simplification
      if (simpl==TRUE) {
        if (GPRparse=="meanSum") {
          print("multiple edges simplification by meanSum")
          # Average enzymes acting in the same reactions and sum up enzymes acting in different reactions
          inpGraph$met_nodes <- paste(inpGraph$from, inpGraph$to,sep = "_")
          inpGraph$HMR_met_nodes <- paste(inpGraph$Reaction, inpGraph$met_nodes ,sep = "_")
          graph <- data.frame(mean=tapply(inpGraph$expr, inpGraph$HMR_met_nodes, mean))
          graph$nodes <- substring(rownames(graph), 10)
          graph=suppressWarnings(graph %>% separate(nodes, c("from", "to"), "_"))
          rownames(graph)=NULL
          graph=graph[,c(2,3,1)]
          newFileName_simpl <-  paste("meanSum_",basenames[x], sep = "")
        }

        #minMax simplification
        if (GPRparse=="minMax") {
          print("multiple edges simplification by minMax")
          #inpGraph <- inpGraph[order(inpGraph$gene),]
          gpr_exp=stri_replace_all_fixed(gpr[,2], inpGraph$gene, inpGraph$expr, vectorize_all=FALSE)
          gpr_expression=as.data.frame(gpr_exp)
          rownames(gpr_expression)=gpr$Rxn.name
          number.of.or <- max(str_count(gpr_expression$gpr_exp, "or"))
          gpr_expr_split_or=suppressWarnings(separate(data=gpr_expression, col=gpr_exp, into=paste0("expr", 1:(number.of.or+1)), sep="\\ or ", extra = "drop"))
          colsOfInterest <- colnames(gpr_expr_split_or[,1:ncol(gpr_expr_split_or)])
          gpr_expr_split_or[colsOfInterest] <- suppressWarnings(lapply(gpr_expr_split_or[colsOfInterest], function(x)
            sapply(strsplit(as.character(x), "and", fixed = TRUE),
                   function(x) min(as.numeric(x)))))
          gpr_expr_split_or$max <- do.call(pmax, c(gpr_expr_split_or, na.rm = TRUE))
          gpr_minmax=as.data.frame(gpr_expr_split_or$max)
          rownames(gpr_minmax)=rownames(gpr_expr_split_or)
          colnames(gpr_minmax)[1]="gpr_minMax"
          inpGraph$gpr_minmax=gpr_minmax$gpr_minMax[match(inpGraph$Reaction,rownames(gpr_minmax))]
          inpGraph$met_nodes <- paste(inpGraph$from, inpGraph$to,sep = "_")
          inpGraph$HMR_met_nodes <- paste(inpGraph$Reaction, inpGraph$met_nodes ,sep = "_")
          graph <- subset.data.frame(inpGraph,select=c("from","to","gpr_minmax"))
          newFileName_simpl <-  paste("minMax_",basenames[x], sep = "")
        }

        #minSum simplification
        if (GPRparse=="minSum") {
          print("multiple edges simplification by minSum")
          #inpGraph <- inpGraph[order(inpGraph$gene),]
          gpr_exp=stri_replace_all_fixed(gpr[,2], inpGraph$gene, inpGraph$expr, vectorize_all=FALSE)
          gpr_expression=as.data.frame(gpr_exp)
          rownames(gpr_expression)=gpr$Rxn.name
          number.of.or <- max(str_count(gpr_expression$gpr_exp, "or"))
          gpr_expr_split_or=suppressWarnings(separate(data=gpr_expression, col=gpr_exp, into=paste0("expr", 1:(number.of.or+1)), sep="\\ or ", extra = "drop"))
          colsOfInterest <- colnames(gpr_expr_split_or[,1:ncol(gpr_expr_split_or)])
          gpr_expr_split_or[colsOfInterest] <- suppressWarnings(lapply(gpr_expr_split_or[colsOfInterest], function(x)
            sapply(strsplit(as.character(x), "and", fixed = TRUE),
                   function(x) min(as.numeric(x)))))
          gpr_expr_split_or$sum <- rowSums(gpr_expr_split_or,na.rm = TRUE)
          gpr_minsum=as.data.frame(gpr_expr_split_or$sum)
          rownames(gpr_minsum)=rownames(gpr_expr_split_or)
          colnames(gpr_minsum)[1]="gpr_minsum"
          inpGraph$gpr_minsum=gpr_minsum$gpr_minsum[match(inpGraph$Reaction,rownames(gpr_minsum))]
          inpGraph$met_nodes <- paste(inpGraph$from, inpGraph$to,sep = "_")
          inpGraph$HMR_met_nodes <- paste(inpGraph$Reaction, inpGraph$met_nodes ,sep = "_")
          graph <- subset.data.frame(inpGraph,select=c("from","to","gpr_minsum"))
          newFileName_simpl <-  paste("minSum_",basenames[x], sep = "")
        }
        #conversion to graph object
        colnames(graph) <- c("from", "to", "weight")
        graphMet <- graph_from_data_frame(graph, directed = directed)
        #additional simplification by sum for multiple edges
        graphMet <- simplify(graphMet, remove.loops = TRUE, remove.multiple = TRUE, edge.attr.comb = "sum")
        graph_simpl=igraph::as_data_frame(graphMet)
        graph_list[[x]] <-graph_simpl
        fileNameout_simpl <- file.path(subdir, paste0(newFileName_simpl, '.', outFormat))
        write.graph(graphMet, fileNameout_simpl, outFormat)
        #output dataframe of simplified weights
        graph_simpl_df <- suppressWarnings(Reduce(function(...) merge(..., by = c('from', 'to')), graph_list))
        suppressWarnings(names(graph_simpl_df)[-(1:2)]<-basenames)
        write.table(graph_simpl_df,
                    paste(subdir,"simplified_weights_df.tsv",sep="/"),
                    append = FALSE,
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)

      }
    }

  }
}



#' Met2RxnGraph.R
#'
#'
#' This function loads a metabolic model in SBML format and returns a metabolic graph having reactions as nodes and edges connecting two if one produces a metabolite consumed by the other one.
#'
#' Builds a graph from the metabolic model.
#' The reactions are connected by metabolites.
#' Two reactions are connected if one produces a metabolite, which is a substrate of the other one.
#' These graphs are unweighted.
#' Recurring metabolites (e.g. ATP, CO2, NADP) can be removed to avoid unrealistic paths setting the parameter "rmMets"==TRUE and using the list of these compounds contained in the data of the package by default or one provided by the user through parameter "recMets".
#' If the parameter add_rev_rxn is set to TRUE, the edges will be assigned considering both directions of reversible reactions.
#' The resulting graph can be built directed or not. Given the rules behind the definition of connections it is naturally directed but the user can choose to not consider this aspect, setting the "directed" parameter to FALSE.

#' @param infile Path to the input SBML file
#' @param rmMets logical value indicating whether the recurring metabolites must be removed. Default TRUE.
#' @param recMets_list a list of recurring metabolite names. The format must be the same of the model.
#' If not provided and rmMets=TRUE the list used is the one provided by the package "system.file("extdata", "metNames.txt", package = "Met2Graph")".
#' @param outDir path to output directory. If not provided outputs will be written in "system.file("extdata/output/", package = "Met2Graph")".
#' @param outFormat character string indicating the format of the graph to export. Default "ncol". For other formats and details please see write.graph function of igraph package.
#' @param directed Logical value, indicating if the graph must be created directed. Default TRUE.
#' @param add_rev_rxn Logical value, whether or not reversible reactions must be considered and edges written in both directions. If FALSE it takes into account the direction as present in the stoichiometric matrix, for general models it is usually forced from left to right.
#'  In case of condition-specific models obtained through FBA this is not necessary. Default FALSE.
#' @importFrom Hmisc %nin%
#' @importFrom igraph graph_from_data_frame write.graph
#' @importFrom stats na.omit
#' @importFrom utils data read.table write.table
#' @return A metabolic graph in tabular format ("from","to","Metabolite") and graph format as chosen by the user.
#' @export
#'
#' @examples
#' \dontrun{
#' infile <- system.file("extdata", "kidney.xml", package = "Met2Graph")
#' outDir <- system.file("extdata/output/", package = "Met2Graph")
#' Met2RxnGraph(infile, rmMets = TRUE, outDir=outDir, outFormat = "graphml")}

Met2RxnGraph <- function(infile, rmMets, recMets_list, outDir, outFormat, directed, add_rev_rxn) {


  # check for necessary specs
  if(!exists("infile")) {
    print("ERROR: no file in input")
    stop()
  }

  if(missing(rmMets)) {
    print("WARNING: rmMets is missing, default TRUE")
    rmMets= TRUE
  } else {
    rmMets=rmMets }

  if (missing(recMets_list) && rmMets==FALSE) {
    recMets=NULL
  }
  else if (missing(recMets_list) && rmMets==TRUE) {
    print ("using the list of recurring metabolite names provided by the package")
    recMets_list= read.delim(system.file("extdata", "metNames.txt", package = "Met2Graph"))
    recMets_list=recMets_list }
  else if (!missing(recMets) && rmMets==TRUE) {
    recMets_list=recMets_list
  }

  if(missing(outDir)) {
    print ("WARNING: No outDir path provided. Output directory will be written in extdata/output")
    outDir=system.file("extdata/output/", package = "Met2Graph")
  } else {
    outDir=outDir
  }

  if(missing(outFormat)) {
    print("WARNING: outFormat is missing, default ncol")
    outFormat= "ncol"
  } else {
    outFormat=outFormat
  }
  if(missing(directed)) {
    print("WARNING: directed parameter is not indicated, default TRUE")
    directed= TRUE
  } else {
    directed=directed
  }
  if(missing(add_rev_rxn)) {
    print("WARNING: add_rev_rxn parameter is not indicated, default FALSE")
    add_rev_rxn=FALSE
  } else {
    add_rev_rxn=add_rev_rxn
  }

  print("Reactions-based graph building...")
  my_mod<-MetExtract(infile)
  subdir<-paste(outDir,"RxnGraphs",sep="/")
  dir.create(subdir)
  mod=my_mod[[1]]
  mod_name=my_mod[[2]]
  product_ind=as.data.frame(my_mod[[3]])
  reactant_ind=as.data.frame(my_mod[[4]])
  #reversible reactions
  if (add_rev_rxn==TRUE){
    colnames(product_ind)[1]="met"
    colnames(reactant_ind)[1]="met"
    product_ind$Reversible=mod@react_rev[match(as.data.frame(product_ind)$Rxn,mod@react_id)]
    reactant_ind$Reversible=mod@react_rev[match(as.data.frame(reactant_ind)$Rxn,mod@react_id)]
    reactant_ind=rbind(as.data.frame(reactant_ind),
                       as.data.frame(product_ind) %>%
                         filter(Reversible==TRUE))
    product_ind=rbind(as.data.frame(product_ind),
                      as.data.frame(reactant_ind) %>%
                        filter(Reversible==TRUE)) }
  #merge reactions by metabolite
  merge_pos_neg=merge(product_ind,reactant_ind,by.x="product",by.y="reactant")
  inpGraph=merge_pos_neg %>%
    rename(
      Metabolite = product,
      from = Rxn.x,
      to = Rxn.y
    )
  inpGraph=inpGraph[,c(2,3,1)]
  #remove recurring metabolites
  if (rmMets==TRUE){
    print("removing recurring metabolites...")
    mets=data.frame(metsName=mod@met_name,metsID=mod@met_id)
    recMets=subset.data.frame(mets,tolower(mets$metsName) %in% tolower(recMets_list[,1]))
    print("the following recurring metabolites have been removed:")
    print(unique(recMets$metsName))
    inpGraph = inpGraph[which(inpGraph$Metabolite  %nin% recMets[,2]), ]
  }
    fileNameout <- file.path(subdir, paste0(mod_name,"_reactions_based_graph.tsv"))
    write.table(inpGraph, fileNameout, sep="\t")
    graphMet <- graph_from_data_frame(inpGraph, directed = TRUE)
    fileNameout <- file.path(subdir, paste0(mod_name,"_reactions_based_graph", '.',outFormat))
    write.graph(graphMet, fileNameout, outFormat)
}


#' Met2EnzGraph
#'
#'
#' This function loads a metabolic model in SBML format and returns a metabolic graph having enzymes as nodes and edges connecting them if they catalyze two reactions which produce and consume a metabolite, respectively.
#'
#' The enzymes catalyzing reactions are connected by metabolites.
#' Two enzymes are connected if they catalyze reactions which produce and consume a metabolite, respectively.
#' The graphs are unweighted.
#' Enzymes for each reaction are annotated according to GPR rules. Enzymes related by AND are complexes and the considered as a single node, while OR relationships are split as different nodes.
#' The gene-protein-reaction (GPR) associations are needed. It can be provided by the user through the parameter "gpr". If not provided by the user, it is automatically extracted form the model, if not present in the model, the one from the human generic GEM is used.
#' Recurring metabolites (e.g. ATP, CO2, NADP) can be removed to avoid unrealistic paths setting the parameter "rmMets"==TRUE and using the list of these compounds contained in the data of the package by default or one provided by the user through parameter "recMets".
#' If the parameter add_rev_rxn is set to TRUE, the edges will be assigned considering both directions of reversible reactions.
#' The resulting graph can be built directed or not. Given the rules behind the definition of connections it is naturally directed but the user can choose to not consider this aspect, setting the "directed" parameter to FALSE.
#'

#' @param infile Path to the input SBML file
#' @param rmMets rmMets logical value indicating whether the recurring metabolites must be removed. Default TRUE.
#' @param recMets_list a list of recurring metabolite names. The format must be the same of the model. If not provided and rmMets=TRUE the list used is the one provided by the package "system.file("extdata", "metNames.txt", package = "Met2Graph")".
#' @param outDir path to output directory. If not provided outputs will be written in "system.file("extdata/output/", package = "Met2Graph")".
#' @param outFormat character string indicating the format of the graph to export. Default "ncol". For other formats and details please see write.graph function of igraph package.
#' @param directed Logical value, indicating if the graph must be created directed. Default TRUE.
#' @param add_rev_rxn Logical value, whether or not reversible reactions must be considered and edges written in both directions. If FALSE it takes into account the direction as present in the stoichiometric matrix, for general models it is usually forced from left to right.
#' In case of condition-specific models obtained through FBA this is not necessary. Default FALSE.
#' @importFrom igraph graph_from_data_frame write.graph
#' @importFrom Hmisc %nin%
#' @import stringi
#' @import stringr
#' @importFrom tidyr separate gather
#' @importFrom dplyr rename select inner_join
#' @importFrom stats na.omit
#' @importFrom utils data read.table write.table
#' @return A metabolic graph in tabular format ("from","to","Metabolite") and igraph format chosen.
#' @export
#'
#' @examples
#' \dontrun{
#' infile <- system.file("extdata", "kidney.xml", package = "Met2Graph")
#' outDir <- system.file("extdata/output/", package = "Met2Graph")
#' Met2EnzGraph(rmMets=TRUE, outDir=outDir, outFormat="graphml")}

Met2EnzGraph <- function(infile, rmMets, recMets_list, outDir, outFormat, directed, add_rev_rxn) {


  # check for necessary specs
  if(!exists("infile")) {
    print("ERROR: no file in input")
    stop()
  }

  if(missing(rmMets)) {
    print("WARNING: rmMets is missing, default TRUE")
    rmMets= TRUE
  } else {
    rmMets=rmMets }

  if (missing(recMets_list) && rmMets==FALSE) {
    recMets=NULL
  }
  else if (missing(recMets_list) && rmMets==TRUE) {
    print ("using the list of recurring metabolite names provided by the package")
    recMets_list= read.delim(system.file("extdata", "metNames.txt", package = "Met2Graph"))
    recMets_list=recMets_list }
  else if (!missing(recMets_list) && rmMets==TRUE) {
    recMets_list=recMets_list
  }

  if(missing(outDir)) {
    print ("WARNING: No outDir path provided, the output directory will be written in extdata/output/")
    outDir=system.file("extdata/output/", package = "Met2Graph")
  } else {
    outDir=outDir
  }

  if(missing(outFormat)) {
    print("WARNING: outFormat is missing, default ncol")
    outFormat= "ncol"
  } else {
    outFormat=outFormat
  }

  if(missing(directed)) {
    print("WARNING: directed parameter is not indicated, default TRUE")
    directed= TRUE
  } else {
    directed=directed
  }
  if(missing(add_rev_rxn)) {
    print("WARNING: add_rev_rxn parameter is not indicated, default FALSE")
    add_rev_rxn=FALSE
  } else {
    add_rev_rxn=add_rev_rxn
  }

  print("Enzymes-based graph building...")
  my_mod<-MetExtract(infile)
  subdir<-paste(outDir,"EnzGraphs",sep="/")
  dir.create(subdir)
  mod=my_mod[[1]]
  mod_name=my_mod[[2]]
  product_ind=as.data.frame(my_mod[[3]])
  reactant_ind=as.data.frame(my_mod[[4]])
  colnames(product_ind) [1] = "met"
  colnames(reactant_ind) [1] = "met"
  #extract GPR
  print("extracting gpr rules")
  print ("automatic extraction of GPR from the model")

  if (length(mod@gpr) == 0){
    print("gpr values are missing in the model, those extracted from human HMR2 will be used")
    gpr = read.delim(system.file("extdata", "gpr.txt", package = "Met2Graph"),
                     sep = "\t")
  }else{
    gpr=data.frame(Rxn=mod@react_id,GPR=mod@gpr)
    if (all(sapply(gpr[,2], function(x)all(x==""))) == TRUE)  {
      print("gpr values are missing in the model, those extracted from human HMR2 will be used")
      gpr=read.delim(system.file("extdata", "gpr.txt", package = "Met2Graph"),sep="\t")
    }
  }
  
  gpr[,2]=as.character(gpr[,2])
  gpr[,2]=gsub("[()]", "", gpr[,2])
  if (add_rev_rxn==TRUE){
    product_ind$Reversible=mod@react_rev[match(as.data.frame(product_ind)$Rxn,mod@react_id)]
    reactant_ind$Reversible=mod@react_rev[match(as.data.frame(reactant_ind)$Rxn,mod@react_id)]
    product_ind_Rev=product_ind[product_ind$Reversible=="TRUE",]
    reactant_ind_Rev=reactant_ind[reactant_ind$Reversible=="TRUE",]
    reactant_ind=rbind(reactant_ind,
                       product_ind_Rev)
    product_ind=rbind(product_ind,
                      reactant_ind_Rev)}
  #merge reactions by metabolite
  merge_pos_neg=merge(product_ind,reactant_ind,by="met")
  inpGraph=merge_pos_neg %>%
    rename(
      Metabolite = met,
      Rxn1 = Rxn.x,
      Rxn2 = Rxn.y
    )
  #replace reactions with corresponding enzymes
  gpr_enzymes_rxn1 <- stri_replace_all_fixed(inpGraph$Rxn1, gpr[,1], gpr[,2], vectorize_all=FALSE)
  gpr_enzymes_rxn2<- stri_replace_all_fixed(inpGraph$Rxn2, gpr[,1], gpr[,2], vectorize_all=FALSE)
  gpr_enzymes=cbind(gpr_enzymes_rxn1,gpr_enzymes_rxn2)
  rm(gpr_enzymes_rxn1,gpr_enzymes_rxn2)
  rownames(gpr_enzymes)=make.unique(as.character(inpGraph$Metabolite), sep = ".")
  gpr_enzymes[gpr_enzymes==""] <- NA
  gpr_enzymes<-na.omit(gpr_enzymes)
  colnames(gpr_enzymes)=c("from","to")
  #split edges if enzymes are in OR relationship
  number.of.or <- max(str_count(gpr_enzymes[,"from"], "or"))
  print("splitting enzymes by OR relationships")
  gpr_enzymes_split_or_rxn1=as.data.frame(gpr_enzymes[,"from"], row.names=row.names(gpr_enzymes))
  colnames(gpr_enzymes_split_or_rxn1)="from"
  gpr_enzymes_split_or_rxn1<-suppressWarnings(separate(data=gpr_enzymes_split_or_rxn1, col=from, into=paste0("from", 1:(number.of.or+1)), sep=" or ", extra = "drop"))
  gpr_enzymes_split_or_rxn1$Metabolite=row.names(gpr_enzymes_split_or_rxn1)
  gpr_enzymes_split_or_rxn1_gathered <- gather(gpr_enzymes_split_or_rxn1, "key", "from", -Metabolite)
  rm(gpr_enzymes_split_or_rxn1)
  gpr_enzymes_split_or_rxn1_gathered=gpr_enzymes_split_or_rxn1_gathered[,]
  gpr_enzymes_split_or_rxn1_gathered=dplyr::select(gpr_enzymes_split_or_rxn1_gathered,-key)
  gpr_enzymes_split_or_rxn1_gathered[gpr_enzymes_split_or_rxn1_gathered==""] <- NA
  gpr_enzymes_split_or_rxn1_gathered<-gpr_enzymes_split_or_rxn1_gathered[!is.na(gpr_enzymes_split_or_rxn1_gathered$from),]
  gpr_enzymes_split_or_rxn2=as.data.frame(gpr_enzymes[,"to"],row.names=row.names(gpr_enzymes))
  colnames(gpr_enzymes_split_or_rxn2)="to"
  gpr_enzymes_split_or_rxn2=suppressWarnings(separate(data=gpr_enzymes_split_or_rxn2, col=to, into=paste0("to", 1:(number.of.or+1)), sep=" or ", extra = "drop"))
  gpr_enzymes_split_or_rxn2$Metabolite=row.names(gpr_enzymes_split_or_rxn2)
  gpr_enzymes_split_or_rxn2_gathered <- gather(gpr_enzymes_split_or_rxn2, "key", "to", -Metabolite)
  rm(gpr_enzymes_split_or_rxn2)
  gpr_enzymes_split_or_rxn2_gathered=select(gpr_enzymes_split_or_rxn2_gathered,-key)
  gpr_enzymes_split_or_rxn2_gathered[gpr_enzymes_split_or_rxn2_gathered==""] <- NA
  gpr_enzymes_split_or_rxn2_gathered<-gpr_enzymes_split_or_rxn2_gathered[!is.na(gpr_enzymes_split_or_rxn2_gathered$to),]
  print("building edges..")
  gpr_enzymes_edges<-inner_join(gpr_enzymes_split_or_rxn1_gathered,gpr_enzymes_split_or_rxn2_gathered, by="Metabolite")
  gpr_enzymes_edges=gpr_enzymes_edges[!duplicated(gpr_enzymes_edges[,c('from', 'to')]),]
  met <- sapply(strsplit(gpr_enzymes_edges$Metabolite, '[.]'),"[",1)
  gpr_enzymes_edges$Metabolite=met
  gpr_enzymes_edges[] <- lapply(gpr_enzymes_edges, gsub, pattern="\\(|)", replacement='')
  gpr_enzymes_edges[] <- lapply(gpr_enzymes_edges, gsub, pattern="\\(|)", replacement='')
  gpr_enzymes_edges=gpr_enzymes_edges[!grepl("HMR",gpr_enzymes_edges$from),]
  gpr_enzymes_edges=gpr_enzymes_edges[!grepl("HMR",gpr_enzymes_edges$to),]
  gpr_enzymes_edges = gpr_enzymes_edges[,c(2,3,1)]
  inpGraph = gpr_enzymes_edges
  #remove recurring metabolites
  if (rmMets==TRUE){
    print("removing recurring metabolites...")
    mets=data.frame(metsName=mod@met_name,metsID=mod@met_id)
    recMets=subset.data.frame(mets,tolower(mets$metsName) %in% tolower(recMets_list[,1]))
    print("the following recurring metabolites have been removed:")
    print(unique(recMets$metsName))
    inpGraph = inpGraph[which(gpr_enzymes_edges$Metabolite  %nin% recMets[,2]), ]
  }
  #conversion to graph object
    fileNameout <- file.path(subdir, paste0(mod_name,"_enzymes_based_graph.tsv"))
    write.table(inpGraph, fileNameout, sep="\t")
    graphMet <- graph_from_data_frame(inpGraph, directed = TRUE)
    fileNameout <- file.path(subdir, paste0(mod_name,"_enzymes_based_graph", '.',outFormat))
    write.graph(graphMet, fileNameout, outFormat)
}

#' Plot_graphs
#'
#'
#' This function plots graphs to pdf.

#' @param inDir path to the input directory.
#' @param outDir path to the output directory.
#' @param format Character string indicating the format of the graph to read. Default "ncol". For other formats and details please see read.graph function of igraph package.
#' @param ... Additional plotting parameters. See igraph.plotting for the complete list.
#' @import igraph
#' @importFrom tools file_path_sans_ext
#' @return A pdf file.
#' @export
#'
#' @examples
#' \dontrun{
#' inDir <- system.file("extdata/output/MetGraphs/", package = "Met2Graph")
#' outDir <- system.file("extdata/output/", package = "Met2Graph")
#' Plot_graphs(inDir=inDir, outDir=outDir, format="ncol", vertex.size= 1,edge.arrow.size=0.1,vertex.label=NA, layout=layout_with_graphopt)}

Plot_graphs <- function(inDir, outDir, format,  ...) {


  # check for necessary specs
  if (missing(inDir)) {
    print("ERROR: please indicate input directory")
    stop()}
  if (missing(outDir)) {
    print("ERROR: please indicate output directory")
    stop()}
  if(missing(format)) {
    print("WARNING: format is missing, default ncol")
    format= "ncol"
  } else {
    format=format
  }

  subdir<-paste(outDir,"Plots",sep="/")
  dir.create(subdir)

  # read graphs:
  graphs <- list.files(path=inDir,
                       pattern = paste0("*.",format), full.names = T, recursive = FALSE)
  count <- 0
  for (x in graphs) {
    count <- count + 1
    print(count)
    g <- read.graph(x, format=format)
    pdf(file.path(subdir, paste0(file_path_sans_ext(basename(x)),".pdf")),
        width = 4,
        height = 4)
    plot.igraph(g, ...)
    dev.off()
  }

}

