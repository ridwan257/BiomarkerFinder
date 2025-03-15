library(sva)
library(dplyr)
library(readr)
library(tibble)
library(limma)

setwd("~/Reuben/Biomarker Finding/trial4/preprocessing")
source('../../libexpr/analysis2.R')
## =============================================================================
## ---------------------------- DATA LOADING GSE -------------------------------
## =============================================================================
g_edata <- read_edata('../mutual-data/gse_prepr_edata.tsv')
g_mdata <- read_mdata('../mutual-data/gse_prepr_mdata.tsv')

m_edata <- read_edata('../mutual-data/emtab_prepr_edata.tsv')
m_mdata <- read_mdata('../mutual-data/emtab_prepr_mdata.tsv')


all(rownames(m_edata) == rownames(g_edata))

table(m_mdata$condition)
table(g_mdata$condition)
## ---------- remove outlier
# rownames(g_mdata)
# g_mdata <- g_mdata[-2,]
# g_edata <- g_edata[,rownames(g_mdata)]
# 
# rownames(m_mdata)
# m_mdata <- m_mdata[-c(30, 36),]
# m_edata <- m_edata[,rownames(m_mdata)]

## ------------- normalize the data
par(mfrow=c(1,1))
boxplot(g_edata)
plotDensities(g_edata, group = g_mdata$condition)
plotDensities(m_edata, group = m_mdata$condition)

# g_edata <- normalizeBetweenArrays(g_edata, method = 'quantile')
# m_edata <- normalizeBetweenArrays(m_edata, method = 'quantile')

## ------------- Adding batch info
m_mdata$batch <- 1
g_mdata$batch <- 2


## ------------- merge the data
edata <- cbind(m_edata, g_edata)
mdata <- rbind(
	m_mdata[,c('age', 'sex', 'condition', 'batch')],
	g_mdata[,c('age', 'sex', 'condition', 'batch')]
)
table(mdata$condition)

all(colnames(edata) == rownames(mdata))
par(mar=c(6, 4, 3, 2))
boxplot(edata, col = mdata$batch+1, las=2, main = "Boxplot of Merged Data (Before)")
plotDensities(edata, group = mdata$batch)


# plotMDS(edata, pch = mdata$batch, col=as.numeric(as.factor(mdata$condition)), top=Inf)
par(mfrow=c(2, 2), mar=c(6, 5, 2, 1))
plot_pca(edata, mdata$condition, xlim=c(-65, 110), lpos = 'topleft', cex=0.8,
		 main = "PCA of Merged Data (Before)")
# legend("topright", c("GSE", "EMTAB", "SNC", "SC"), pch=c(2, 1, 19, 19), col=c(1, 1, 2, 1))
## --------------------- design
design_mat <- model.matrix(~condition, data = mdata)

# edata_limma <- removeBatchEffect(
# 	edata_norm,
# 	batch = mdata$batch,
# 	design = design_mat
# )
# plotDensities(edata_limma, group = mdata$batch)
# plotMDS(edata_limma, pch = mdata$batch, col=as.numeric(as.factor(mdata$condition)))


edata_combat <- ComBat(
	dat = edata,
	batch = mdata$batch, 
	mod = design_mat, 
	par.prior = TRUE,
	prior.plots = FALSE
)

plotDensities(edata_combat, group = mdata$condition)
boxplot(edata_combat, col = mdata$batch+1, las=2, main = "Boxplot of Merged Data (After)")


plot_pca(edata_combat, mdata$condition, xlim=c(-100, 100), lpos = 'topleft',
		 main = "PCA of Merged Data (After)", cex=0.8)
# plotMDS(edata_combat, pch = mdata$batch, col=as.numeric(as.factor(mdata$condition)), top = Inf)
# legend("topright", c("GSE", "EMTAB", "SNC", "SC"), pch=c(2, 1, 19, 19), col=c(1, 1, 2, 1))



clst <- hclust( dist(t(edata_combat)), method = 'average' )
plot(clst)

h_genes <- c('GAPDH', 'ACTB', 'B2M', 'TBP')
par(mfrow=c(2, 4), mar=c(3, 5, 4, 3))
for(g in h_genes){
	boxplot(edata[g, ] ~ mdata$batch, main = paste(g, " - before"))
	boxplot(edata_combat[g, ] ~ mdata$batch, main = paste(g, " - after"))
}




## ------------- differential gene expression
w <- arrayWeights(edata_combat, design = design_mat)
fit <- lmFit(edata_combat, design = design_mat)
fit2 <- eBayes(fit)



result <- topTable(fit2, coef = 2, number = Inf) %>% 
	rownames_to_column('Gene')


plot_volcano(result, lab=T, minFC = -0.8, maxFC = 0.8, by='P.Value')

par(mfrow=c(1, 1))
plotMD(fit2)
abline(h=0, lwd=2, col='red')

addmargins( table(mdata$batch ,mdata$condition) )

## --------------------------------- Write The Results
edata_combat %>% 
	as.data.frame() %>% 
	rownames_to_column('Gene') %>% 
	write_tsv('../mutual-data/merged_edata.tsv')

mdata %>% 
	as.data.frame() %>% 
	rownames_to_column('Sample') %>% 
	write_tsv('../mutual-data/merged_mdata.tsv')









# death_data <- readxl::read_excel("~/Reuben/Thesis_final/data/death_rate.xlsx")
# 
# death_data %>% 
# 	tidyr::pivot_longer(
# 		cols = c('Male', 'Female'),
# 		names_to = 'Gender',
# 		values_to = 'Death Rate'
# 	) %>% 
# 	ggplot(aes(x=Year, y=`Death Rate`, colour = Gender)) +
# 	geom_line(linewidth=1) +
# 	# geom_point() +
# 	theme_bw() +
# 	scale_x_continuous(
# 		limits = c(1930, 2025),
# 		breaks = seq(1930, 2030, by = 10),
# 		# minor_breaks = seq(1930, 2030, by = 5),
# 		expand = expansion(add = c(0, 0), mult = c(0, 0.015))
# 	) +
# 	theme(
# 		# panel.grid.major.x = element_line(color = "grey80"),
# 		# panel.grid.minor.x = element_line(color = "grey90"),
# 	)
















