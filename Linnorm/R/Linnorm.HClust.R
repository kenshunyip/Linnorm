#' Linnorm-hierarchical clustering analysis.
#'
#' This function first performs Linnorm transformation on the dataset. Then, it will perform hierarchical clustering analysis.
#' @param datamatrix	The matrix or data frame that contains your dataset. Each row is a feature (or Gene) and each column is a sample (or replicate). Raw Counts, CPM, RPKM, FPKM or TPM are supported. Undefined values such as NA are not supported. It is not compatible with log transformed datasets.
#' @param RowSamples	Logical. In the datamatrix, if each row is a sample and each row is a feature, set this to TRUE so that you don't need to transpose it. Linnorm works slightly faster with this argument set to TRUE, but it should be negligable for smaller datasets. Defaults to FALSE.
#' @param MZP Double >=0, <= 1. Minimum non-Zero Portion Threshold for this function. Genes not satisfying this threshold will be removed from HVG anlaysis. For exmaple, if set to 0.3, genes without at least 30 percent of the samples being non-zero will be removed. Defaults to 0.
#' @param DataImputation	Logical. Perform data imputation on the dataset after transformation. Defaults to TRUE.
#' @param input	Character. "Raw" or "Linnorm". In case you have already transformed your dataset with Linnorm, set input into "Linnorm" so that you can input the Linnorm transformed dataset into the "datamatrix" argument. Defaults to "Raw".
#' @param method_hclust	Charcter. Method to be used in hierarchical clustering. (From hclust {fastcluster}: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".) Defaults to "ward.D".
#' @param method_dist	Charcter. Method to be used in hierarchical clustering. (From Dist {amap}: the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "correlation", "spearman" or "kendall". Any unambiguous substring can be given.) Defaults to "pearson".
#' @param  Group	Character vector with length equals to sample size. Each character in this vector corresponds to each of the columns (samples) in the datamatrix. If this is provided, sample names will be colored according to their group. Defaults to NULL.
#' @param num_Clust	Integer >= 0. Number of clusters in hierarchical clustering. No cluster will be highlighted if this is set to 0. Defaults to 0.
#' @param Color	Character vector. Color of the groups/clusters in the plot. This vector must be as long as num_Clust, or Group if it is provided. Defaults to "Auto".
#' @param ClustRect	Logical. If num_Clust > 0, should a rectangle be used to highlight the clusters? Defaults to TRUE.
#' @param RectColor	Character. If ClustRect is TRUE, this controls the color of the rectangle. Defaults to "red".
#' @param fontsize	Numeric. Font size of the texts in the figure. Defualts to 0.5.
#' @param linethickness	Numeric. Controls the thickness of the lines in the figure. Defaults to 0.5.
#' @param plot.title	Character. Set the title of the plot. Defaults to "Hierarchical clustering".
#' @param ... arguments that will be passed into Linnorm's transformation function.
#' @details  This function performs PCA clustering using Linnorm transformation.
#' @return It returns a list with the following objects:
##' \itemize{
##'  \item{Results:}{ If num_Clust > 0, this outputs a named vector that contains the cluster assignment information of each sample. Else, this outputs a number 0.}
##'  \item{plot:}{ Plot of hierarchical clustering.}
##'  \item{Linnorm:}{ Linnorm transformed data matrix.}
##' }
#' @keywords Linnorm RNA-seq Raw Count Expression RPKM FPKM TPM CPM normalization transformation Parametric hierarchical Clustering
#' @export
#' @examples
#' #Obtain example matrix:
#' data(Islam2011)
#' #Example:
#' HClust.results <- Linnorm.HClust(Islam2011, Group=c(rep("ESC",48), rep("EF",44), rep("NegCtrl",4)))

Linnorm.HClust <- function(datamatrix, RowSamples = FALSE, MZP = 0, DataImputation = TRUE, input="Raw", method_hclust="ward.D", method_dist="pearson", Group=NULL, num_Clust=0, Color="Auto", ClustRect=TRUE, RectColor="red", fontsize=0.5, linethickness=0.5, plot.title="Hierarchical clustering", ...) {
	#Hierarchical clustering with Linnorm transformed dataset
	#Author: (Ken) Shun Hang Yip <shunyip@bu.edu>
	#Note from http://stackoverflow.com/questions/24140339/tree-cut-and-rectangles-around-clusters-for-a-horizontal-dendrogram-in-r last accessed Feb 9th, 2017
	if (input != "Raw" && input != "Linnorm") {
		stop("input argument is not recognized.")
	}
	if (MZP > 1 || MZP < 0) {
		stop("Invalid MZP.")
	}
	if (length(Group) != 0) {
		if (length(Group) != length(datamatrix[1,])) {
			stop("Group must be a vector with the same length as sample size.")
		}
	}
	if (num_Clust < 0) {
		stop("Invalid number of clusters.")
	}
	if (!is.logical(RowSamples)){
		stop("Invalid RowSamples.")
	}
	if (!is.logical(DataImputation)){
		stop("Invalid DataImputation.")
	}
	if (!is.logical(ClustRect)){
		stop("Invalid ClustRect.")
	}
	if (!RowSamples) {
		datamatrix <- t(datamatrix)
	}
	#Linnorm transformation
	if (input == "Raw") {
		datamatrix <- Linnorm(datamatrix, DataImputation=DataImputation, RowSamples = TRUE,...)
	}
	#Backup data that will be filtered, so that we can include them in the output
	Backup <- colSums(datamatrix != 0) < nrow(datamatrix) * MZP
	Backup2 <- 0
	if (sum(Backup) != 0) {
		Backup2 <-  datamatrix[,Backup]
	}
	#Filter zeroes based on MZP threshold
	datamatrix <- datamatrix[,colSums(datamatrix != 0) >= nrow(datamatrix) * MZP]
	
	#Clustering
	hc <- hclust(Dist(datamatrix, method = method_dist), ,method = method_hclust)
	dendr <- dendro_data(hc, type = "rectangle")
	
	#plot object
	render_plot <- 0
	
	#Render Color of plot
	colorCode <- 0
	if (length(Color) == 0) {
		Color = "Auto"
	}
	if (length(Color) == 1) {
		if (Color != "Auto") {
			if (!areColors(Color)) {
				stop("Invalid Color.")
			}
			if (length(unique(Group)) > 1) {
				stop("Number of Color provided does not equal to the number of groups.")
			}
			if (num_Clust != 1 && length(unique(Group)) == 0) {
				stop("Number of Color provided does not equal to num_Clust.")
			}
			colorCode <- c("grey60", Color)
		}
		if (Color == "Auto") {
			if (length(unique(Group)) > num_Clust) {
				colorCode <- c("grey60", rainbow(length(unique(Group))))
			} else {
				colorCode <- c("grey60", rainbow(num_Clust))
			}
		}
	} else {
		if (!areColors(Color)) {
			stop("Invalid Color.")
		}
		if (length(unique(Group)) > num_Clust) {
			if (length(Color) != length(unique(Group))) {
				stop("Number of Color provided does not equal to the number of groups.")
			}
		} else {
			if (length(Color) != num_Clust) {
				stop("Number of Color provided does not equal to num_Clust.")
			}
		}
		colorCode <- c("grey60", Color)
	}
	
	#Cluster object
	clust <- 0
	
	if (num_Clust > 0) {
		clust <- cutree(hc, k = num_Clust)
		# Split dendrogram into upper grey section and lower coloured section
		height <- unique(dendr$segments$y)[order(unique(dendr$segments$y), decreasing = TRUE)]
		cut.height <- mean(c(height[num_Clust], height[num_Clust-1]))
		dendr$segments$line <- ifelse(dendr$segments$y == dendr$segments$yend &
		   dendr$segments$y > cut.height, 1, 2)
		dendr$segments$line <- ifelse(dendr$segments$yend  > cut.height, 1, dendr$segments$line)

		# Number the clusters
		dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
		change <- which(dendr$segments$cluster == 1)
		for (i in 1:num_Clust) dendr$segments$cluster[change[i]] = i + 1
		dendr$segments$cluster <-  ifelse(dendr$segments$line == 1, 1, 
					 ifelse(dendr$segments$cluster == 0, NA, dendr$segments$cluster))
		dendr$segments$cluster <- na.locf(dendr$segments$cluster)
		
		#plotting
		render_plot <- ggplot() + 
		geom_segment(data = segment(dendr),
			aes(x=x, y=y, xend=xend, yend=yend, size=factor(line), colour=factor(cluster)), size=linethickness, lineend = "square", show.legend = FALSE
		) +
		scale_colour_manual(values = colorCode) +
		scale_y_reverse(expand = c(0.2, 0)) + 
		labs(x = NULL, y = NULL) +
		coord_flip() +
		ggtitle(plot.title) +
		theme(axis.line.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.title.y = element_blank(),
			panel.background = element_rect(fill = "white"),
			panel.grid = element_blank()
		)
	} else {
		#Plotting without clusters
		render_plot <- ggplot() + 
		geom_segment(data = segment(dendr),
			aes(x=x, y=y, xend=xend, yend=yend, size=factor(line)), size=linethickness, lineend = "square", show.legend = FALSE
		) +
		scale_y_reverse(expand = c(0.2, 0)) + 
		labs(x = NULL, y = NULL) +
		coord_flip() +
		theme(axis.line.y = element_blank(),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.title.y = element_blank(),
			panel.background = element_rect(fill = "white"),
			panel.grid = element_blank()
		)
	}
	
	if (length(Group) != 0) {
		#Set label color
		labcol <- as.numeric(label(dendr)[,3])
		unilab <- as.character(unique(as.character(Group)))
		index <- 1
		for (i in length(unilab):1) {
			labcol[which(as.character(dendr$label$label) %in% rownames(datamatrix)[which(Group == unilab[i])])] <- index+1
			index <- index + 1
		}
		render_plot <- render_plot + geom_text(data = label(dendr), 
			aes(x, y, label = label, colour=factor(labcol)),
			hjust = -0.2, size = fontsize, show.legend = FALSE
		)
	} else {
		render_plot <- render_plot + geom_text(data = label(dendr), 
			aes(x, y, label = label),
			hjust = -0.2, size = fontsize, show.legend = FALSE
		)
	}
	
	if (ClustRect) {
		if (length(rownames(datamatrix)) > length(unique(rownames(datamatrix)))) {
			warning("Duplicate sample names found. Rectangle not drawn.")
		} else {
			#rectangle
			clust.df <- data.frame(label=rownames(datamatrix), cluster=factor(clust))
			dendr2 <- merge(dendr[["labels"]],clust.df, by="label")
			rect <- aggregate(x~cluster,dendr2,range)
			rect <- data.frame(rect$cluster,rect$x)
			ymax <- mean(hc$height[length(hc$height)-((num_Clust-2):(num_Clust-1))])
			render_plot <- render_plot + geom_rect(data=rect, 
			aes(xmin=X1-.5, xmax=X2+.5, ymin=0, ymax=ymax), 
			color=RectColor, fill=NA, size=linethickness)
		}
	}
	render_plot <- ggplot_build(render_plot)
	#Reconstruct Linnorm transformed matrix for output
	if (sum(Backup) != 0) {
		datamatrix <- cbind(datamatrix, Backup2)
	}
	if (!RowSamples) {
		datamatrix <- t(datamatrix)
	}
	#Results for output
	listing <- list(clust, render_plot, datamatrix)
	results <- setNames(listing, c("Results", "plot", "Linnorm"))
	return (results)
}
