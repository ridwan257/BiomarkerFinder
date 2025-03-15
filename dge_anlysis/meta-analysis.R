library(dplyr)
library(readr)


setwd("~/Reuben/Biomarker Finding/trial4/dge")


m_result <- read_tsv('../mutual-data/emtab_top_table.tsv')
g_result <- read_tsv('../mutual-data/gse_top_table.tsv')


## ------------------------- Fisher’s Combined Probability Test
merged_df <- inner_join(m_result, g_result, by = "Gene", suffix = c("_m", "_g"))
merged_df <- merged_df %>%
	# filter(
	# 	abs(logFC_m) > 0.1 | abs(logFC_g) > 0.1
	# ) %>% 
	rowwise() %>%
	mutate(
		# Combine p-values using Fisher's method
		fisher_stat = -2 * (log(P.Value_m) + log(P.Value_g)),
		fisher_pval = pchisq(fisher_stat, df = 4, lower.tail = FALSE)  # df = 2*2=4
	) %>%
	ungroup()


## ------------------------ Stouffer’s Z-Score Method
# merged_df <- merged_df %>%
# 	mutate(
# 		# Convert p-values to Z-scores with sign from logFC
# 		z1 = qnorm(1 - P.Value_m/2) * sign(logFC_m),
# 		z2 = qnorm(1 - P.Value_g/2) * sign(logFC_g),
# 		
# 		# Combine Z-scores
# 		stouffer_z = (z1 + z2) / sqrt(2),
# 		stouffer_pval = 2 * pnorm(-abs(stouffer_z))  # Two-tailed p-value
# 	)
merged_df <- merged_df %>%
	mutate(
		# Directional z-scores
		z_dir_m = sign(logFC_m) * qnorm(1 - P.Value_m/2),
		z_dir_g = sign(logFC_g) * qnorm(1 - P.Value_g/2),
		
		# Weights (magnitude of logFC)
		weight_m = 1, #3.26, #abs(logFC_m),
		weight_g = 1, #2.97, #abs(logFC_g),
		
		# Combine using Stouffer's with explicit weights
		stouffer_z = (weight_m * z_dir_m + weight_g * z_dir_g) / sqrt(weight_m^2 + weight_g^2),
		stouffer_pval = 2 * pnorm(-abs(stouffer_z))
	)
## -------------------------- Adjust for Multiple Testing
merged_df <- merged_df %>%
	mutate(
		fisher_fdr = p.adjust(fisher_pval, method = "fdr"),
		stouffer_fdr = p.adjust(stouffer_pval, method = "fdr")
	)

merged_df <- merged_df %>% 
	mutate(
		logFC_mean = (logFC_m + logFC_g) / 2
	)
## -------------------------- Filter
sig_genes_meta <- merged_df %>%
	# filter(stouffer_fdr < 0.1) %>%
	filter(stouffer_pval < 0.05) %>%
	pull(Gene)

sort(sig_genes_meta)

# hist(merged_df$stouffer_pval, breaks=100)


# par(mfrow=c(1,1))
# hist(merged_df$P.Value_m, breaks = 20)
# hist(merged_df$P.Value_g, breaks = 20)
# hist(merged_df$stouffer_fdr, breaks = 20)

# merged_df %>%
# 	select(
# 		Gene, logFC_m, logFC_g, P.Value_m, P.Value_g, logFC_mean,
# 		stouffer_z, stouffer_pval, stouffer_fdr
# 		#z_dir_m, z_dir_g, weight_m, weight_g,
# 		#stouffer_fdr#, logFC_mean
# 	) %>%
# 	arrange(stouffer_pval) %>%
# 	mutate(
# 		across(where(is.numeric), function(x){round(x, 4)})
# 	) %>%
# 	writexl::write_xlsx('stouffer_result.xlsx')














