library(limma)
library(tibble)

setwd("~/Reuben/Biomarker Finding/trial4/dge")
source('../../libexpr/analysis2.R')
## =============================================================================
## ---------------------------- DATA LOADING GSE -------------------------------
## =============================================================================
g_edata <- read_edata('../mutual-data/gse_prepr_edata.tsv')
g_mdata <- read_mdata('../mutual-data/gse_prepr_mdata.tsv')

m_edata <- read_edata('../mutual-data/emtab_prepr_edata.tsv')
m_mdata <- read_mdata('../mutual-data/emtab_prepr_mdata.tsv')


all(rownames(m_edata) == rownames(g_edata))

## ---------- DAMPLES
par(mfrow=c(1, 1))
clst <- hclust(dist(t(m_edata)), method='average')
plot(clst)

clst <- hclust(dist(t(g_edata)), method='average')
plot(clst)

plotMDS(g_edata, col=as.numeric(as.factor(g_mdata$condition)), top = Inf)

## ------------- normalize the data
# g_edata <- normalizeBetweenArrays(g_edata, method = 'quantile')
# m_edata <- normalizeBetweenArrays(m_edata, method = 'quantile')

## --------------- gene filter
# --------- EMTAB
m_edata_filt <- m_edata
# m_mean_expr <- rowMeans(m_edata)
# summary(m_mean_expr)
# keep1 <- m_mean_expr > quantile(m_mean_expr, 0.25)
# table(keep1)
# m_edata_filt <- m_edata[keep1, ]

# m_median <- apply(m_edata, 1, median)
# summary(m_median)
# keep1 <- rowSums(m_edata > quantile(m_median, 0.25)) >= 1
# table(keep1)
# m_edata_filt <- m_edata[keep1, ]

# m_median <- median(m_edata)
# m_median
# keep1 <- rowSums(m_edata > m_median) >= 1
# table(keep1)
# m_edata_filt <- m_edata[keep1, ]


# --------- GSE
g_edata_filt <- g_edata
# g_mean_expr <- rowMeans(g_edata)
# summary(g_mean_expr)
# keep1 <- g_mean_expr > quantile(g_mean_expr, 0.25)
# table(keep1)
# g_edata_filt <- g_edata[keep1, ]

# g_median <- apply(g_edata, 1, median)
# summary(g_median)
# keep1 <- rowSums(g_edata > quantile(g_median, 0.25)) >= 1
# table(keep1)
# g_edata_filt <- g_edata[keep1, ]

# g_median <- median(g_edata)
# g_median
# keep1 <- rowSums(g_edata > g_median) >= 1
# table(keep1)
# g_edata_filt <- g_edata[keep1, ]



# ----------- common
length( intersect(rownames(g_edata_filt), rownames(m_edata_filt)) )


## ---------------- differential gene expression
g_design_mat <- model.matrix(~condition, data = g_mdata)
m_design_mat <- model.matrix(~condition, data = m_mdata)

## --------------- EMTAB
m_w <- arrayWeights(m_edata_filt, design = m_design_mat)
# m_fit <- lmFit(m_edata, design=m_design_mat, weights = m_w)
m_fit <- lmFit(m_edata_filt, design=m_design_mat)
m_fit2 <- eBayes(m_fit)


m_result <- topTable(m_fit2, coef = 2, number = Inf) %>% 
	rownames_to_column('Gene')

## --------------- GSE
g_w <- arrayWeights(g_edata_filt, design = g_design_mat)
# g_fit <- lmFit(g_edata, design=g_design_mat, weights = g_w)
g_fit <- lmFit(g_edata_filt, design=g_design_mat)
g_fit2 <- eBayes(g_fit)


g_result <- topTable(g_fit2, coef = 2, number = Inf) %>% 
	rownames_to_column('Gene')



## --------------- plot data
plot_volcano(m_result, highlight=15, show_line = T, legend = T) +
	ggtitle("Dataset E-MTAB-1690")
limma::plotMA(m_fit2, status=case_when(
	m_fit2$coefficients[,2] > 1 & m_fit2$p.value[,2] < 0.05 ~ 'up',
	m_fit2$coefficients[,2] < -1 & m_fit2$p.value[,2] < 0.05 ~ 'down',
	TRUE ~ 'same'
), main = 'MA plot E-MTAB-1690', legend='topleft')
abline(h=0, col='red', lwd=2)


plot_volcano(g_result, highlight=15, show_line = T, legend = T) +
	ggtitle("Dataset GSE19027")

limma::plotMA(g_fit2, status=case_when(
	g_fit2$coefficients[,2] > 1 & g_fit2$p.value[,2] < 0.05 ~ 'up',
	g_fit2$coefficients[,2] < -1 & g_fit2$p.value[,2] < 0.05 ~ 'down',
	TRUE ~ 'same'
), main = 'MA plot GSE19027', legend='topleft')
abline(h=0, col='red', lwd=2)


## ----------------- store the results
write_tsv(g_result, "../mutual-data/gse_top_table.tsv")
write_tsv(m_result, "../mutual-data/emtab_top_table.tsv")


## ----------------- checking genes

m_deg <- m_result %>% 
	filter( abs(logFC) > 0.5 & P.Value < 0.05 ) %>% 
	pull('Gene')

m_deg

g_deg <- g_result %>% 
	filter( abs(logFC) > 0.5 & P.Value < 0.05 ) %>% 
	pull('Gene')

g_deg

sort( intersect(m_deg, g_deg) )


































