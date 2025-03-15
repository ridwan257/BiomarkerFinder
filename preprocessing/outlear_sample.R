library(dendextend)
library(readr)
library(tibble)

setwd("~/Reuben/Biomarker Finding/trial4/dge")
source('../../libexpr/analysis2.R')


find_outliar_samples <- function(t_edata, groups, f=2){
	outliear_list <- character(0)
	unique_grp <- unique(groups)
	
	for(grp in unique_grp){
		tmp_df <- t_edata[, groups==grp]
		while(TRUE){
			iac_matrix <- cor(tmp_df, method = "pearson")
			avg_iac <- rowMeans(iac_matrix, na.rm = TRUE)
			mean_iac <- mean(avg_iac, na.rm = TRUE)
			sd_iac <- sd(avg_iac, na.rm = TRUE)
			outliers <- names(avg_iac)[avg_iac < (mean_iac - f * sd_iac)]
			if(length(outliers) == 0){ break }
			outliear_list <- c(outliear_list, outliers)
			tmp_df <- tmp_df[, ! colnames(tmp_df) %in% outliers]
		}
	}
	
	return(outliear_list)
}
## =============================================================================
## ---------------------------- DATA LOADING -------------------------------
## =============================================================================
g_edata <- read_edata('../mutual-data/gse_raw_gene_collapsed.tsv')
g_mdata <- read_mdata('../mutual-data/gse_raw_mdata.tsv')

m_edata <- read_edata('../mutual-data/emtab_raw_gene_collapsed.tsv')
m_mdata <- read_mdata('../mutual-data/emtab_raw_mdata.tsv')



## ---------- OUTLIEAR FINDING EMTAB
par(mfrow=c(2, 2))
clst <- hclust(dist(t(scale(m_edata))), method='average')
dnd <- as.dendrogram(clst, hang = 0.1)
dnd <- color_labels(dnd, col = ifelse(m_mdata$condition[clst$order]=='SC', 1, 2))
plot(dnd, cex = 0.1, ylab='Height', main='Cluster Dendogram E-MTAB-1690 (Before)')
legend('topleft', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)

# plotMDS(m_edata, col=as.numeric(as.factor(m_mdata$condition)), top = Inf, cex = 0.9)
# legend('bottomright', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)
# plot_pca(
# 	m_edata, m_mdata$condition, cex=0.8, lpos = 'bottomleft',
# 	main='PCA plot (Before)'
# )

m_outlier <- find_outliar_samples(m_edata, m_mdata$condition, f=2)
m_outlier
m_mdata_filt <- m_mdata[! rownames(m_mdata) %in% m_outlier, ]
table(m_mdata_filt$condition)

m_edata_filt <- m_edata[, rownames(m_mdata_filt)]


clst <- hclust(dist(t(scale(m_edata_filt))), method='average')
dnd <- as.dendrogram(clst, hang = 0.1)
dnd <- color_labels(dnd, col = ifelse(m_mdata_filt$condition[clst$order]=='SC', 1, 2))
plot(dnd, cex = 0.1, ylab='Height', main='Cluster Dendogram E-MTAB-1690 (After)')
legend('topright', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)

# plotMDS(m_edata_filt, col=as.numeric(as.factor(m_mdata_filt$condition)), top = Inf, cex = 0.9)
# legend('bottomright', c('SC', 'SNC'), col=c(1, 2), pch=19)

# plot_pca(
# 	m_edata_filt, m_mdata_filt$condition, cex=0.8, lpos = 'bottomright', 
# 	main='PCA plot (After)'
# )




## ---------- OUTLIEAR FINDING GSE
# par(mfrow=c(2, 2))
clst <- hclust(dist(t(scale(g_edata))), method='average')
dnd <- as.dendrogram(clst, hang = 0.1)
dnd <- color_labels(dnd, col = ifelse(g_mdata$condition[clst$order]=='SC', 1, 2))
plot(dnd, cex = 0.1, ylab='Height', main='Cluster Dendogram GSE19027 (Before)')
legend('topleft', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)

# plotMDS(m_edata, col=as.numeric(as.factor(m_mdata$condition)), top = Inf, cex = 0.9)
# legend('bottomright', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)
# plot_pca(
# 	g_edata, g_mdata$condition, cex=0.8, lpos = 'topleft',
# 	main='PCA plot (Before)'
# )

# plotMDS(g_edata, col=as.numeric(as.factor(g_mdata$condition)), top = Inf)


g_outlier <- find_outliar_samples(g_edata, g_mdata$condition, f=2)
g_outlier
g_outlier <- c(g_outlier, "GSM470865", "GSM470866")

g_mdata_filt <- g_mdata[! rownames(g_mdata) %in% g_outlier, ]
table(g_mdata_filt$condition)

g_edata_filt <- g_edata[, rownames(g_mdata_filt)]


clst <- hclust(dist(t(scale(g_edata_filt))), method='average')
dnd <- as.dendrogram(clst, hang = 0.1)
dnd <- color_labels(dnd, col = ifelse(g_mdata_filt$condition[clst$order]=='SC', 1, 2))
plot(dnd, cex = 0.1, ylab='Height', main='Cluster Dendogram GSE19027 (After)')
legend('topright', c('SC', 'SNC'), col=c(1, 2), pch=19, cex = 0.8)



plot_pca(
	g_edata_filt, g_mdata_filt$condition, cex=0.8, lpos = 'bottomleft', 
	main='PCA plot (After)'
)

# plotMDS(g_edata_filt, col=as.numeric(as.factor(g_mdata_filt$condition)), top = Inf)


## ---------- OUTLIEAR FINDING GSE
# all(rownames(m_edata) == rownames(g_edata))
# merge_edata <- cbind(m_edata, g_edata)
# g_mdata$batch <- 2
# m_mdata$batch <- 1
# merge_mdata <- rbind(m_mdata[,c('condition', 'batch')], g_mdata[,c('condition', 'batch')])
# all(colnames(merge_edata) == rownames(merge_mdata))

# outliear_samples <- find_outliar_samples(merge_edata, merge_mdata$condition)
# outliear_samples
# merge_mdata_filt <- merge_mdata[!rownames(merge_mdata) %in% outliear_samples, ]
# merge_edata_filt <- merge_edata[,rownames(merge_mdata_filt)]

# table(merge_mdata_filt$condition)
# table(merge_mdata$condition)

# design_mat <- model.matrix(~condition, data = merge_mdata)
# merge_edata_combat <- sva::ComBat(
# 	merge_edata,
# 	batch = merge_mdata$batch,
# 	mod = design_mat
# )
# 
# plotDensities(merge_edata_combat, legend = F)
# 
# outliear_samples <- find_outliar_samples(merge_edata_combat, merge_mdata$condition)
# outliear_samples
# merge_mdata_filt <- merge_mdata[!rownames(merge_mdata) %in% outliear_samples, ]
# merge_edata_filt <- merge_edata[,rownames(merge_mdata_filt)]
# table(merge_mdata_filt$condition)

## ---------- WRITE DATA
g_edata_filt %>%
	as.data.frame() %>%
	rownames_to_column('Gene') %>%
	write_tsv('../mutual-data/gse_prepr_edata.tsv')

m_edata_filt %>%
	as.data.frame() %>%
	rownames_to_column('Gene') %>%
	write_tsv('../mutual-data/emtab_prepr_edata.tsv')

g_mdata_filt %>%
	rownames_to_column('Sample') %>%
	write_tsv('../mutual-data/gse_prepr_mdata.tsv')

m_mdata_filt %>%
	rownames_to_column('Sample') %>%
	write_tsv('../mutual-data/emtab_prepr_mdata.tsv')





















