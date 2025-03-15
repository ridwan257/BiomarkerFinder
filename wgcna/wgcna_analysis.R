library(dplyr)
library(tibble)
library(readr)
library(limma)
library(WGCNA)


allowWGCNAThreads()
setwd("~/Reuben/Biomarker Finding/trial4/wgcna")
source('../../libexpr/analysis2.R')
## =============================================================================
## ---------------------------- DATA LOADING GSE -------------------------------
## =============================================================================
raw_edata <- read_edata('../mutual-data/merged_edata.tsv')
mdata <- read_mdata('../mutual-data/merged_mdata.tsv')



## =============================================================================
## ---------------------- filtering low expr genes -----------------------------
## =============================================================================
edata_filt <- raw_edata
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Variance filters
# m_mad <- apply(raw_edata, 1, mad)
# summary(m_mad)
# keep1 <- m_mad > quantile(m_mad, 0.20)
# table(keep1)


# var_expr <- rowVars(raw_edata)
# summary(var_expr)
# keep1 <- var_expr > quantile(var_expr, 0.50)
# table(keep1)
# m_edata_filt <- m_edata_filt[keep1, ]


# mean_expr <- rowMeans(raw_edata)
# summary(mean_expr)
# keep3 <- mean_expr > 4#quantile(mean_expr, 0.20)
# table(keep3)
# 
# sum(keep1 & keep3)

# edata_filt <- raw_edata[keep1 & keep3, ]
# edata_filt <- edata_filt[keep1, ]

## =============================================================================
## --------------------------- OUTLIER SAMPLES -----------------------------
## =============================================================================
par(mfrow=c(1,1))
plotDensities(edata_filt, group = mdata$condition)
boxplot(edata_filt)

diss_sample <- 1 - adjacency(edata_filt, type = 'distance')

par(mfrow=c(1, 1))
clst <- hclust(as.dist(diss_sample), method = 'average')
dend <- as.dendrogram(clst, hang = 0.1)
dend <- dendextend::color_labels(dend, col = ifelse(mdata$condition[clst$order]=='SC', 1, 2))
plot(dend, cex=0.8)
legend('topright', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)

par(mfrow=c(1, 1), mar=c(6, 5, 4, 3))
plotMDS(edata_filt, top=Inf, col=as.numeric(as.factor(mdata$condition)))


adj_sample <- adjacency(edata_filt, type = 'distance')
simm_sample <- colSums(adj_sample) - 1
simm_sample <- scale(simm_sample)
simm_sample
find_outliar_samples(edata_filt, mdata$condition)
## ------------------------ excluding samples
# rownames(mdata)
# mdata <- mdata[-c(32, 33),]
# edata_filt <- edata_filt[, rownames(mdata)]

## =============================================================================
## --------------------------- TANSPOSING THE DATA -----------------------------
## =============================================================================
edata <- t(edata_filt)
rownames(edata)
dim(edata)
table(mdata$condition)
# edata <- t(raw_edata[sig_genes_meta,])
# all(rownames(edata) == rownames(mdata))
## =============================================================================
## ----------------------------- SOFT thresholds -------------------------------
## =============================================================================
power <- 5:20
sft <- pickSoftThreshold(edata, corFnc=bicor, networkType='signed', powerVector=power)

sft$powerEstimate
par(mfrow = c(1, 2), mar=c(6, 4, 3, 1))
plot(
	sft$fitIndices[, 1],
	-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
	xlab = "Soft Threshold (power)",
	ylab = "SFT, signed R^2",
	type = "p",
	main = "Scale independence",
	pch='•'
)

text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
	 labels = power, col = "red", cex=1)
abline(h = 0.8, col = "red")
abline(h = 0.9, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)",
	 ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = power, col = "red")


## =============================================================================
## ------------------------------- BUILD NETWORK -------------------------------
## =============================================================================
ADJ1 <- adjacency(edata, type='signed', power=14, corFnc = "bicor")
TOMdiss <- TOMdist(ADJ1)

gene_tree <- hclust( as.dist(TOMdiss), method='average' )

tree_dynamic <- cutreeDynamic(
	gene_tree, distM = TOMdiss, deepSplit = 2,  minClusterSize = 50,
	pamStage=TRUE, pamRespectsDendro=FALSE
)
mcolors <- labels2colors(tree_dynamic)
plotDendroAndColors(
	gene_tree, colors = mcolors, dendroLabels = FALSE, 
	addGuide = TRUE, guideHang = 0.05, hang=0.03
)
table(mcolors)
length(unique(mcolors))


## --------------------- MERGE NETWORK
MEs <- moduleEigengenes(edata, mcolors)$eigengenes

MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method='average')

par(mar = c(6, 4, 3, 3), mfrow=c(1, 1))
plot(METree)
abline(h=0.20, col='red')


merge <- mergeCloseModules(edata, mcolors, cutHeight=0.20, verbose = 3)
table(merge$colors)
length(unique(merge$colors))

plotDendroAndColors(
	gene_tree, colors = data.frame(original=mcolors, merged=merge$colors),
	dendroLabels = F, hang=0.03, addGuide=T, guideHang = 0.05
)


## =============================================================================
## ---------------------------------- ANALYSIS ---------------------------------
## =============================================================================
module_col <- merge$colors
# module_col <- mcolors
length(unique(module_col))
table(module_col)
MEs <- moduleEigengenes(edata, module_col)$eigengenes
MEs <- orderMEs(MEs)

## ------------------------- TRAIT CORRELATION
# trait_data <- model.matrix(~0+condition, data=mdata)
# t_race <- model.matrix(~0+race, data=mdata)
# 
trait_data <- data.frame(
	age = mdata$age,
	sex = ifelse(mdata$sex=='Male', 1, 0),
	conditionSNC = ifelse(mdata$condition=='SNC', 1, 0),
	batch = ifelse(mdata$batch==1, 1, 0)
)

moduleTraitCor <- cor(MEs, trait_data)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(mdata))

textMatrix <- paste(
	signif(moduleTraitCor, 2), " (",
	signif(moduleTraitPvalue, 1), ")",
	sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 10, 3, 3), mfrow=c(1, 1))
labeledHeatmap(
	Matrix = moduleTraitCor,
	xLabels = colnames(trait_data),
	yLabels = names(MEs),
	ySymbols = names(MEs),
	colorLabels = FALSE,
	colors = blueWhiteRed(50),
	textMatrix = textMatrix,
	setStdMargins = FALSE,
	cex.text = 0.8,
	zlim = c(-1,1),
	main = paste("Module-trait relationships")
)

table(module_col)

# cat("Module", "Size", "Smoker_healthy (pval)", sep='\t')
# for(i in 1:nrow(textMatrix)){
# 	cat(names(MEs)[i], sum(module_col==substring(names(MEs)[i], 3)), "\t", sep = '\t')
# 	cat(textMatrix[i,3], sep = '\t')
# 	cat('\n')
# }


## ------------------------------ JUST CHECKING GENES
i_genes <- c(
	'ADH7', 'NQO1', 'CYP1A1', 'CYP1B1', 'AKR1B10', 'ALDH1A1', 'ALDH3A1', 'NFE2L2', 
	'NRF1', 'MAFG', 'TP53'
)
for (g in i_genes){
	print(paste(
		g, '->',
		module_col[ which(colnames(edata) == g) ]
	))
}

module_col[ which(colnames(edata) == 'LOXL2') ]

## ------------------------------ WRITTING MODULE GENES
length( colnames(edata)[module_col=="magenta"] )
writeLines(
	colnames(edata)[module_col=="magenta"],
	# con='module.txt'
)

# head(colnames(edata)[module_col=="red"], n=100)
write_clip(colnames(edata)[module_col=="magenta"])



# genes_of_interest <- c("NQO1", "CYP1A1", "CYP1B1", "AKR1B10", "ALDH1A1", "ALDH3A1")
# 
# edata_subset <- t(edata[, genes_of_interest])  # Subset expression data
# edata_melt <- reshape2::melt(edata_filt) %>% 
# 	left_join(
# 		data.frame(Var2=rownames(mdata), condition=mdata$condition),
# 		by='Var2'
# 	)
# 
# ggplot(edata_melt, aes(x=Var2, y=value, fill=condition)) +  
# 	geom_boxplot() + 
# 	theme(axis.text.x = element_text(angle=90, hjust=1)) +
# 	ggtitle("Expression Levels of Key Genes")




## =============================================================================
## --------------------------- CALCULATING METRICS -----------------------------
## =============================================================================
GS1 <- cor(trait_data[,3], edata)
Alldegrees1 <- intramodularConnectivity(ADJ1, module_col)
KMEs <- signedKME(edata, MEs, outputColumnName='MM.')
geneTraitPvalue <- corPvalueStudent(GS1, nrow(mdata))

table(module_col)

# colorlevels <- c(
# 	"yellow", "red", "cyan", 
# 	
# )

colorlevels <- rownames(moduleTraitPvalue)[moduleTraitPvalue[,'conditionSNC'] < 0.05] %>% 
	substring(3)

colorlevels


par( mfrow=c(3, 4), mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))){
	which_module <- colorlevels[[i]]
	restrict1 <- (module_col==which_module)
	verboseScatterplot(
		Alldegrees1$kWithin[restrict1],
		abs(GS1[restrict1]),
		col=module_col[restrict1],
		main=which_module,
		xlab = "Connectivity", ylab = "Gene Significance",
		abline = TRUE
	)
}

## =============================================================================
## ------------------------------ GS vs MM --------------------------------
## =============================================================================

par( mfrow=c(2, 3), mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))){
	which_module=colorlevels[[i]]
	restrict1 = (module_col==which_module)
	# restrict1 <- TRUE
	gene_pval <- geneTraitPvalue[restrict1]
	# pch_list <- ifelse(gene_pval < 0.05, 2, 1)
	sig_genes <- abs(KMEs[restrict1,  paste("MM", which_module, sep='.')]) > 0.8 &
		abs(GS1[restrict1]) > 0.2
	pch_list <- ifelse(sig_genes, 2, 1)
	# pch_list <- ifelse(module_col == which_module, 2, 1)
	col_list <- module_col[restrict1]
	if(which_module=='lightyellow') { col_list <- rep('gold', length(module_col[restrict1])) }
	col_list[sig_genes] <- 'red'
	
	verboseScatterplot(
		x=KMEs[restrict1,  paste("MM", which_module, sep='.')],
		y=GS1[restrict1],
		col=col_list,
		main=which_module,
		xlab = "Module Membership", ylab = "Gene Significance",
		abline = T,
		pch = pch_list
	)
	if(which_module=='darkturquoise'){
		abline(v=0.80, h=0.20, col='red')
	} else{
		abline(v=0.80, h=-0.20, col='red')
	}
	
	# abline(v=c(0.80, -0.80), h=c(0.20, -0.20), col='red')
}

## =============================================================================
## ---------------------------- WRITE THE RESULT -------------------------------
## =============================================================================

selected_modules <- colorlevels
# selected_modules <- c("darkturquoise", "darkgrey", "lightyellow", "black", "yellow")

gene_names_ml <- c()

for(i in 1:length(selected_modules)){
	which_module <- selected_modules[i]
	restrict1 <- module_col == which_module
	
	filtered_gene <- ( abs(GS1) > 0.20 ) & 
		( abs(KMEs[,  paste("MM", which_module, sep='.')]) > 0.80 ) #& restrict1
	
	# table(module_col[filtered_gene])
	# table(filtered_gene)
	print(paste("Module", which_module, '-', sum(filtered_gene)))
	
	gene_names_ml <- c( gene_names_ml, colnames(filtered_gene)[filtered_gene] )
}
length(gene_names_ml)
gene_names_ml <- unique(gene_names_ml)
length(gene_names_ml)

sort( gene_names_ml )



sort( intersect(sig_genes_meta, gene_names_ml) )
# sort( intersect(sig_genes, gene_names_ml) )
# 
# writeLines( intersect(sig_genes_meta, gene_names_ml), con='./upreg_genes.txt' )
# writeLines( intersect(sig_genes, gene_names_ml), con='./upreg_genes_rra.txt' )
writeLines( intersect(sig_genes_meta, gene_names_ml), con='./common_genes.txt' )
# writeLines( intersect(sig_genes_meta, gene_names_ml), con='./common_genes_excluded.txt' )

# writeLines(sort( intersect(sig_genes_meta, gene_names_ml) ))


ggvenn::ggvenn(
	list(
		DGE_analysis = sig_genes_meta,
		WGCNA = gene_names_ml
	),
	fill_color = c("#A7C7E7", "#F7A072"),
	stroke_size = 0.5,
	set_name_size = 4,
	text_size = 4
)







## =============================================================================
## ---------------------------- DRAW NETWORK HUBS -------------------------------
## =============================================================================

# ----------------------------- selecting hub based on gene significance module of interest
which_module <- "lightyellow"
# restrict1 <- module_col == which_module
# restrict1 <- TRUE
# restrict1 <- module_col == which_module


filtered_gene <- ( abs(GS1) > 0.30 ) & 
	( abs(KMEs[,  paste("MM", which_module, sep='.')]) > 0.843 )
table( filtered_gene )

hub_gene_names <- colnames(edata)[filtered_gene]

hub_gene_degree <- Alldegrees1[hub_gene_names, ]


gene_degree <- hub_gene_degree %>% 
	arrange(desc(kWithin)) %>% 
	slice_head(n=10)

# top_hub_gene <- rownames(gene_degree)
top_hub_gene <- hub_gene_names
top_hub_gene

hub_gene_adjacency <- ADJ1[top_hub_gene, top_hub_gene] * 2
# hub_gene_adjacency[hub_gene_adjacency < 1e-5] <- 0
hub_gene_adjacency[hub_gene_adjacency < 0.2] <- 0
hub_graph <- igraph::graph_from_adjacency_matrix(
	hub_gene_adjacency,
	mode = "undirected",
	weighted = TRUE,
	diag = FALSE
)
# igraph::V(hub_graph)$size <- igraph::degree(hub_graph) * 3
# igraph::V(hub_graph)$size <- gene_degree[top_hub_gene, 'kWithin'] * 3
igraph::V(hub_graph)$size <- abs(GS1[,top_hub_gene]) * 60
par(mfrow=c(1, 1), mar=c(4, 3, 2, 3))
plot(
	hub_graph,
	layout = igraph::layout_with_fr,
	vertex.label = igraph::V(hub_graph)$name,
	vertex.label.cex = 1.2,
	vertex.label.font = 1.8,
	vertex.color = adjustcolor("white", alpha.f = 1),
	vertex.frame.color = adjustcolor("black", alpha.f = 0.3),
	edge.color = adjustcolor("black", alpha.f = 0.5),
	edge.width = igraph::E(hub_graph)$weight * 5
)







