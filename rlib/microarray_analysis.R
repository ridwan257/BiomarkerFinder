library(ggrepel)
library(ggplot2)

## outliar finding using median absolute deviation
find_outliar_mad <- function(values, threshold=3){
	median_val <- median(values)
	abs_diff <- abs(values - median_val)
	mad_ <- median(abs_diff)
	outlier_threshold <- threshold * mad_
	
	return(abs_diff > outlier_threshold)
}

## merge duplicate probe
outliarMerge <- function(edata, rowGroup){
	new_edata <- edata %>% 
		as.data.frame() %>% 
		mutate(
			row_mean = rowMeans(.),
			Gene = rowGroup
		) %>%
		tibble::rownames_to_column('SPOT_ID') %>%
		group_by(Gene) %>% 
		filter(
			if( n() == 2 ){
				if ( sum(grepl('^\\d+_at$', SPOT_ID)) == 1 ) { grepl('^\\d+_at$', SPOT_ID) }
				else{ c(T,T) }
			} 
			else if ( n() > 2 ){
				( ! find_outliar_mad(row_mean) )
			}
			else {
				TRUE
			}
		) %>% 
		summarize( across(-SPOT_ID, mean) ) %>% 
		select(everything(), -row_mean) %>% 
		tibble::column_to_rownames('Gene') %>% 
		as.matrix()
	
	return(new_edata)
}



## ---------------- plot valocano
plot_volcano <- function(
		results, minFC=-1, maxFC=1, alpha=0.05, by="P.Value", lab=FALSE,
		highlight = NA, show_line = FALSE, byBvalue = FALSE, a=0.8, 
		fsize = 3, gene_list = NA, show_grid=TRUE, msize=2.5,
		mcolor='black', legend=FALSE, n_overlap = 10
)
{
	# by = "P.Value"
	# alpha = 0.01
	
	min_fc_value <- min(results$logFC)
	max_fc_value <- max(results$logFC)
	max_log_pval <- max(-log10(results[[by]]))
	
	plot_data <- results %>%
		dplyr::select(Gene, logFC, !!sym(by), B) %>%
		dplyr::mutate(logFC := as.numeric(logFC), !!sym(by) := as.numeric(!!sym(by))) %>% 
		dplyr::mutate(
			state = case_when(
				logFC >= maxFC & !!sym(by) <= alpha ~ 'up',
				logFC <= minFC & !!sym(by) <= alpha ~ 'down',
				TRUE ~ 'same'
			)
		)
	
	figure <- ggplot(plot_data, aes(
		x = logFC, 
		y = -log10( !!sym(by) ), 
		color = state)) +
		geom_hline(yintercept = 0, color = "black", linewidth=0.5)
	
	if(show_line){
		figure <- figure +
			geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "red") +
			geom_vline(xintercept = c(minFC, maxFC), linetype = "dashed", color = "red")
	}
	
	figure <- figure +
		geom_point(show.legend = legend, alpha=a) +
		scale_color_manual(values = c(
			'up'= 'blue',
			'down' = 'red',
			'same' = 'gray'),
			breaks = c("down", "same", "up"),
			labels = c("Down-regulated", "Not-Significant", "Up-regulated")
		)
	
	
	if(!is.na(highlight)){
		if(byBvalue) { plot_data <- arrange(plot_data, desc(B)) }
		else { plot_data <- arrange( plot_data, !!sym(by), desc(abs(logFC)) ) }
		
		figure <- figure + 
			geom_point(
				data = plot_data[1:highlight, ],
				# color = mcolor,
				# size = msize
			) +
			geom_text_repel(
				data = plot_data[1:highlight, ],
				aes(label = Gene),
				size = fsize,
				fontface = 'bold',
				max.overlaps = n_overlap,
				colour = 'black',
				show.legend = NA
			)
	}
	
	else if (!all(is.na(gene_list))){
		mark_gene_df <- dplyr::filter(plot_data, Gene %in% gene_list)
		
		figure <- figure + 
			geom_point(
				data = mark_gene_df,
				color = mcolor,
				size = msize
			) +
			geom_text_repel(
				data = mark_gene_df,
				aes(label = Gene),
				size = fsize,
				fontface = 'bold',
				max.overlaps = n_overlap,
				colour = 'black',
				show.legend = FALSE
			)
	}
	
	else if (lab){
		figure <- figure + geom_text_repel(
			data = . %>% dplyr::filter(state !=  'same'),
			aes(label = Gene),
			size = fsize,
			fontface = 'bold',
			max.overlaps = n_overlap,
			colour = 'black',
			show.legend = FALSE
		)
	}
	
	
	figure <- figure +
		labs(x = "logFC", y = paste0("-log10(", by, ")", collapse = '')) +
		scale_y_continuous(limits = c(0, NA), breaks = seq(0, ceiling(max_log_pval), 1)) + 
		scale_x_continuous(breaks = seq(floor(min_fc_value), ceiling(max_fc_value), 1)) +
		theme_bw() +
		theme(
			legend.position = "top",
			legend.title = element_blank(),
			legend.text = element_text(size = 12),
			panel.border = element_rect(color = 'black', fill = NA, linewidth = 1.5)
		)
	
	if(!show_grid){
		figure <- figure +
			theme(
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank()
			)
	}
	
	return(figure)
}




