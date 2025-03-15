library(affy)
library(hgu133acdf)
library(hgu133plus2cdf)
library(limma)
library(readr)
library(tibble)

setwd("~/Reuben/Biomarker Finding/trial4/preprocessing/")

compute_score <- function(calls, group) {
	apply(calls[, group], 1, function(row) {
		sum(ifelse(row == "P", 1, ifelse(row == "M", 0.51, 0)))
	})
}

threshold <- 0.50
## ====================================================================================
## ----------------------------------- GSE DATA ---------------------------------------
## ====================================================================================
g_cel_files <- list.files(path="../mutual-data/gse", pattern="*.CEL", full.names=TRUE)
g_raw_data <- ReadAffy(filenames=g_cel_files)
g_mdata <- read_tsv('../mutual-data/gse_raw_mdata.tsv')

g_mas5_data <- mas5calls(g_raw_data)
g_pma_matrix <- exprs(g_mas5_data)
colnames(g_pma_matrix)
colnames(g_pma_matrix) <- gsub('\\..*', '', colnames(g_pma_matrix))

all(g_mdata$Sample == colnames(g_pma_matrix))

g_sc_scores <- compute_score(g_pma_matrix, g_mdata$condition=='SC')
g_snc_scores <- compute_score(g_pma_matrix, g_mdata$condition=='SNC')

g_sc_threshold <- threshold * sum(g_mdata$condition == "SC")
g_snc_threshold <- threshold * sum(g_mdata$condition == "SNC") 

g_passing_probes <- g_sc_scores >= g_sc_threshold | g_snc_scores >= g_snc_threshold
table(g_passing_probes)


## ==================exprs()## ====================================================================================
## ----------------------------------- EMTAB DATA ---------------------------------------
## ====================================================================================
m_cel_files <- list.files(path="../mutual-data/emtab", pattern="*.CEL", full.names=TRUE)
m_raw_data <- ReadAffy(filenames=m_cel_files)
m_mdata <- read_tsv('../mutual-data/emtab_raw_mdata.tsv')

m_mas5_data <- mas5calls(m_raw_data)
m_pma_matrix <- exprs(m_mas5_data)
colnames(m_pma_matrix)
colnames(m_pma_matrix) <- gsub('\\..*', '', colnames(m_pma_matrix))

all(m_mdata$Sample == colnames(m_pma_matrix))

m_sc_scores <- compute_score(m_pma_matrix, m_mdata$condition=='SC')
m_snc_scores <- compute_score(m_pma_matrix, m_mdata$condition=='SNC')

m_sc_threshold <- threshold * sum(m_mdata$condition == "SC")
m_snc_threshold <- threshold * sum(m_mdata$condition == "SNC") 

m_passing_probes <- m_sc_scores >= m_sc_threshold | m_snc_scores >= m_sc_threshold
table(m_passing_probes)




## ====================================================================================
## ----------------------------------- RMA normalization ---------------------------------------
## ====================================================================================
## ----------------- GSE
g_rma_data <- rma(g_raw_data)

g_rma_data_filt <- g_rma_data[g_passing_probes, ]
g_edata <- exprs(g_rma_data_filt)

## ----------------- EMTAB
m_rma_data <- rma(m_raw_data)

m_rma_data_filt <- m_rma_data[m_passing_probes, ]
m_edata <- exprs(m_rma_data_filt)

## ------------------------------------------------
dim(m_edata)
dim(g_edata)
length(intersect(rownames(m_edata), rownames(g_edata)))

## ----------------------------------- Feature data
# GSE
g_probe_names <- featureNames(g_rma_data_filt)
g_gene_map <- AnnotationDbi::select(
	hgu133a.db::hgu133a.db, 
	keys = g_probe_names, 
	columns = "SYMBOL", 
	keytype = "PROBEID"
)
g_gene_map <- g_gene_map[match(g_probe_names, g_gene_map$PROBEID),]
g_gene_map <- g_gene_map[!is.na(g_gene_map$SYMBOL), ]
sum(g_gene_map$SYMBOL == '')
g_edata <- g_edata[g_gene_map$PROBEID, ]

# EMTAB
m_probe_names <- featureNames(m_rma_data_filt)
m_gene_map <- AnnotationDbi::select(
	hgu133plus2.db::hgu133plus2.db,
	keys = m_probe_names, 
	columns = "SYMBOL", 
	keytype = "PROBEID"
)
m_gene_map <- m_gene_map[match(m_probe_names, m_gene_map$PROBEID),]
sum(is.na(m_gene_map$SYMBOL))
m_gene_map <- m_gene_map[!is.na(m_gene_map$SYMBOL), ]
sum(m_gene_map$SYMBOL == '')
m_edata <- m_edata[m_gene_map$PROBEID, ]


length(unique(g_gene_map$SYMBOL))
length(unique(m_gene_map$SYMBOL))

length(intersect(m_gene_map$SYMBOL, g_gene_map$SYMBOL))

## ----------------------------------- Some plots
# plotDensities(g_edata, group = g_mdata$condition)
# plotDensities(m_edata, group = m_mdata$condition)
# boxplot(m_edata)
## ----------------------------------- Collapsing genes
final_edata_gse <- avereps(g_edata, ID=g_gene_map$SYMBOL)
# agg_res <- WGCNA::collapseRows(
# 	g_edata,
# 	rowGroup=g_gene_map$SYMBOL,
# 	rowID=g_gene_map$PROBEID,
# 	method="MaxMean"
# )
# final_edata_gse <- agg_res$datETcollapsed

final_edata_emtab <- avereps(m_edata, ID=m_gene_map$SYMBOL)
# agg_res <- WGCNA::collapseRows(
# 	m_edata,
# 	rowGroup=m_gene_map$SYMBOL,
# 	rowID=m_gene_map$PROBEID,
# 	method="MaxMean"
# )
# final_edata_emtab <- agg_res$datETcollapsed


print(length(intersect(rownames(final_edata_emtab), rownames(final_edata_gse))))
## ----------------------------------- Common Genes
common_genes <- intersect(rownames(final_edata_emtab), rownames(final_edata_gse))

final_edata_emtab <- final_edata_emtab[common_genes, ]
final_edata_gse <- final_edata_gse[common_genes, ]

colnames(final_edata_emtab) <- gsub('\\..*', '', colnames(final_edata_emtab))
colnames(final_edata_gse) <- gsub('\\..*', '', colnames(final_edata_gse))
## ---------------------------------- Writing the result
final_edata_gse |>
	as.data.frame() |>
	rownames_to_column('Gene') |>
	write_tsv('../mutual-data/gse_raw_gene_collapsed.tsv')

final_edata_emtab |>
	as.data.frame() |>
	rownames_to_column('Gene') |>
	write_tsv('../mutual-data/emtab_raw_gene_collapsed.tsv')











