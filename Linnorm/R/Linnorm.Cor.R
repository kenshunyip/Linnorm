#' Linnorm-gene correlation network analysis.
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will perform correlation network analysis on the dataset.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param RowSamples	Logical. In the datamatrix, if each row is a sample and each row is a feature, set this to TRUE so that you don't need to transpose it. Linnorm works slightly faster with this argument set to TRUE, but it should be negligable for smaller datasets. Defaults to FALSE.
#' @param input	Character. "Raw" or "Linnorm". In case you have already transformed your dataset with Linnorm, set input into "Linnorm" so that you can put the Linnorm transformed dataset into the "datamatrix" argument. Defaults to "Raw".
#' @param method	Character. "pearson", "kendall" or "spearman". Method for the calculation of correlation coefficients. Defaults to "pearson".
#' @param MZP	Double >=0, <= 1. Minimum non-Zero Portion Threshold for this function. Genes not satisfying this threshold will be removed for correlation calculation. For exmaple, if set to 0.3, genes without at least 30 percent of the samples being non-zero will be considered for this study. Defaults to 0.5.
#' @param sig.q	Double >=0, <= 1. Only gene pairs with q values less than this threshold will be included in the "Results" data frame. Defaults to 0.05.
#' @param plotNetwork	Logical. Should the program output the network plot to a file? An "igraph" object will be included in the output regardless. Defaults to TRUE. 
#' @param plotNumPairs	Integer >= 50. Number of gene pairs to be used in the network plot. Defaults to 5000.
#' @param plotdegree	Integer >= 0. In the network plot, genes (vertices) without at least this number of degree will be removed. Defaults to 0.
#' @param plotname	Character. Name of the network plot. File extension will be appended to it. Defaults to "networkplot".
#' @param plotformat	Character. "pdf" or "png". Network plot output format. Defaults to "png".
#' @param plotVertexSize	Double >0. Controls vertex Size in the network plot. Defaults to 1.
#' @param plotFontSize	Double >0. Controls font Size in the network plot. Defaults to 1.
#' @param plot.Pos.cor.col	Character. Color of the edges of positively correlated gene pairs. Defaults to "red".
#' @param plot.Neg.cor.col	Character. Color of the edges of negatively correlated gene pairs. Defaults to "green".
#' @param vertex.col	Character. "cluster" or a color. This controls the color of the vertices. Defaults to "cluster".
#' @param plotlayout	Character. "kk" or "fr". "kk" uses Kamada-Kawai algorithm in igraph to assign vertex and edges. It scales edge length with correlation strength. However, it can cause overlaps between vertices. "fr" uses Fruchterman-Reingold algorithm in igraph to assign vertex and edges. It prevents overlatps between vertices better than "kk", but edge lengths are not scaled to correlation strength. Defaults to "kk".
#' @param clusterMethod	Character. "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_infomap", "cluster_label_prop", "cluster_leading_eigen", "cluster_louvain", "cluster_optimal", "cluster_spinglass" or "cluster_walktrap". These are clustering functions from the igraph package. Defaults to "cluster_edge_betweenness".
#' @param ... arguments that will be passed into Linnorm's transformation function.
#' @details  This function performed gene correlated study in the dataset by using Linnorm transformation
#' @return This function will output a list with the following objects:
##' \itemize{
##'  \item{Results:}{ A data frame containing the results of the analysis, showing only the significant results determined by "sig.q" (see below).}
##'  \item{Cor.Matrix:}{ The resulting correlation matrix between each gene. }
##'  \item{q.Matrix:}{ A matrix of q values of each of the correlation coefficient from Cor.Matrix. }
##'  \item{Cluster:}{ A data frame that shows which gene belongs to which cluster.}
##'  \item{igraph:}{ The igraph object for users who want to draw the network plot manually. }
##'  \item{Linnorm:}{ Linnorm transformed data matrix.}
##' }
#' @return The "Results" data frame has the following columns:
##' \itemize{
##'  \item{Gene1:}{ Name of gene 1.}
##'  \item{Gene2:}{ Name of gene 2.}
##'  \item{Cor:}{ Correlation coefficient between the two genes.}
##'  \item{p.value:}{ p value of the correlation coefficient.}
##'  \item{q.value:}{ q value of the correlation coefficient.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric correlation coefficient kendall pearson spearman
#' @export
#' @examples
#' data(Islam2011)
#' #Analysis on Islam2011 embryonic stem cells
#' results <- Linnorm.Cor(Islam2011[,1:48], plotNetwork=FALSE)
Linnorm.Cor <- function(datamatrix,RowSamples = FALSE, input="Raw", method = "pearson", MZP=0.5, sig.q=0.05, plotNetwork=TRUE, plotNumPairs=5000, plotdegree=0, plotname="networkplot", plotformat = "png", plotVertexSize=1, plotFontSize=1, plot.Pos.cor.col="red", plot.Neg.cor.col="green", vertex.col="cluster", plotlayout="kk", clusterMethod = "cluster_edge_betweenness", ...) {
	#Correlation network analysis by using Linnorm transformed data
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	if (input != "Raw" && input != "Linnorm") {
		stop("input argument is not recognized.")
	}
	if (method != "pearson" && method != "spearman"&& method != "kendall") {
		stop("method is not recognized.")
	}
	if (sig.q < 0 || sig.q > 1) {
		stop("Invalid sig.q value.")
	}
	if (plotNumPairs <50) {
		stop("plotNumPairs is too small.")
	}
	if (plotdegree < 0) {
		stop("Invalid plotdegree value.")
	}
	if (plotformat != "pdf" && plotformat != "png") {
		stop("plotformat is not recognized.")
	}
	if (plotVertexSize <= 0) {
		stop("Invalid plotVertexSize.")
	}
	if (plotFontSize <= 0) {
		stop("Invalid plotFontSize.")
	}
	if (plotlayout != "kk" && plotlayout != "fr") {
		stop("Invalid plotlayout.")
	}
	if (!areColors(vertex.col) && vertex.col != "cluster") {
		stop("Invalid vertex.col.")
	}
	igraphFun <- c("cluster_edge_betweenness", "cluster_fast_greedy", "cluster_infomap", "cluster_label_prop", "cluster_leading_eigen", "cluster_louvain", "cluster_optimal", "cluster_spinglass", "cluster_walktrap")
	if (!(clusterMethod %in% igraphFun)) {
		stop("Invalid clusterMethod.")
	}
	if (!is.logical(RowSamples)){
		stop("Invalid RowSamples.")
	}
	if (!is.logical(plotNetwork)){
		stop("Invalid plotNetwork.")
	}
	#Check data format
	if (!RowSamples) {
		datamatrix <- t(datamatrix)
	}
	#Linnorm transformation
	if (input == "Raw") {
		datamatrix <- Linnorm(datamatrix, RowSamples = TRUE, ...)
	}
	#Backup data that will be filtered, so that we can include them in the output
	Backup <- colSums(datamatrix != 0) < nrow(datamatrix) * MZP
	Backup2 <- 0
	if (sum(Backup) != 0) {
		Backup2 <-  datamatrix[,Backup]
	}
	#Filter zeroes based on MZP threshold
	datamatrix <- datamatrix[,colSums(datamatrix != 0) >= nrow(datamatrix) * MZP]
	
	#Calculate correlation coefficients
	correlation <- cor(datamatrix, method=method)
	datamatrix <- datamatrix[,colnames(correlation)]
	correlations <- correlation[upper.tri(correlation,diag=FALSE)]
	
	#Index for locating genes with index in "correlations"
	#Note that if you reverse column 1 and column 2 in index, it becomes lower index.
	index <- createUpperIndex(ncol(correlation), length(correlations))
	
	#Calculate p and q values for each correlation coefficient
	pvalues <- r.sig(correlations, nrow(datamatrix))
	qvalues <- p.adjust(pvalues,"BH")
	
	#Construct q value matrix
	qvaluematrix <- UpperToMatrix(qvalues,index)
	
	#Result matrix
	wanted <- which(qvalues <= sig.q)
	resultmatrix <- data.frame(Gene1=rownames(correlation)[index[wanted,1]],Gene2=rownames(correlation)[index[wanted,2]],Cor=correlations[wanted],p.value=pvalues[wanted],q.value=qvalues[wanted])
	
	
	#Network for the top "plotNumPairs" significantly correlated gene pairs.
	AllGene1 <- as.character(resultmatrix[,1])
	AllGene2 <- as.character(resultmatrix[,2])
	
	if (plotNumPairs > length(resultmatrix[,5])) {
		plotNumPairs <- length(resultmatrix[,5])
	}
	
	orderbyCor <- order(resultmatrix[,4],decreasing=FALSE)
	nodes <- vector(mode="character",length(plotNumPairs) * 2)
	index <- 1
	for (i in 1:plotNumPairs) {
		nodes[index] <- AllGene1[orderbyCor[i]]
		index <- index + 1
		nodes[index] <- AllGene2[orderbyCor[i]]
		index <- index + 1
	}
	
	g1 <- graph(nodes,directed=FALSE)
	Thislayout <- 0
	
	positive <- which(resultmatrix[orderbyCor[1:plotNumPairs],3] > 0)
	negative <- which(resultmatrix[orderbyCor[1:plotNumPairs],3] < 0)
	E(g1)[positive]$color <- "red"
	E(g1)[negative]$color <- "green"
	if (plotlayout == "kk") {
		E(g1)$weight <- (resultmatrix[orderbyCor[1:plotNumPairs],3] - 3)^2
		Thislayout <- layout_with_kk(g1)
	}
	if (plotlayout == "fr") {
		E(g1)$weight <- resultmatrix[orderbyCor[1:plotNumPairs],3] + 1
		Thislayout <- layout_with_fr(g1)
	}
	
	#Clustering
	clustering <- 0
	if (clusterMethod == "cluster_edge_betweenness") { clustering <- cluster_edge_betweenness(g1)}
	if (clusterMethod == "cluster_fast_greedy") { clustering <- cluster_fast_greedy(g1)}
	if (clusterMethod == "cluster_infomap") { clustering <- cluster_infomap(g1)}
	if (clusterMethod == "cluster_label_prop") { clustering <- cluster_label_prop(g1)}
	if (clusterMethod == "cluster_leading_eigen") { clustering <- cluster_leading_eigen(g1)}
	if (clusterMethod == "cluster_louvain") { clustering <- cluster_louvain(g1)}
	if (clusterMethod == "cluster_optimal") { clustering <- cluster_optimal(g1)}
	if (clusterMethod == "cluster_spinglass") { clustering <- cluster_spinglass(g1)}
	if (clusterMethod == "cluster_walktrap") { clustering <- cluster_walktrap(g1)}
	
	if (vertex.col == "cluster") {
		vertex.col <- clustering$membership
	}
	
	#Cluster results
	Clust.res <- data.frame(Gene=clustering$names,Cluster=clustering$membership)
	
	#Plot the correlation network
	if (plotNetwork) {
		if(plotformat == "pdf") {
			pdf(paste(plotname,".pdf",sep=""),width = 10, height = 10)
		}
		if (plotformat == "png"){
			png(paste(plotname,".png",sep=""),res=2000, width = 10, height = 10, units = 'in')
		}
		plot1 <- plot(g1, vertex.color=vertex.col, vertex.size=0.8 * plotVertexSize,vertex.frame.color="transparent", vertex.label.color="black",vertex.label.cex=0.03 * plotFontSize, edge.curved=0.1,edge.width=0.05 * plotVertexSize, layout=Thislayout, margin=0)
		print(plot1)
		dev.off()
	}
	#Reconstruct Linnorm transformed matrix
	if (sum(Backup) != 0) {
		datamatrix <- cbind(datamatrix, Backup2)
	}
	if (!RowSamples) {
		datamatrix <- t(datamatrix)
	}
	#Prepare for result output
	listing <- list(resultmatrix,correlation,qvaluematrix,Clust.res,g1,datamatrix)
	result <- setNames(listing, c("Results", "Cor.Matrix", "q.Matrix", "Cluster", "igraph", "Linnorm"))
	return (result)
}
