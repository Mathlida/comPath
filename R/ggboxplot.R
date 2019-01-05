ggboxplot <- function(pathdata,pathname){
	pathdata %>% 
		group_by(type) %>%
		mutate(outlier=ifelse(is_outlier(fc),genes,as.numeric(NA))) %>%
		ggplot(aes(x=type, fc)) + 
		xlab("Data type") + ylab("Fold chage") +
		ggtitle(substr(pathname,1,8)) +
		theme_bw() +
		theme(plot.title = element_text(hjust = 0.5)) +
		geom_boxplot(outlier.colour = NA) +
		geom_beeswarm(aes(color=fc),size=1.5) +
		geom_text_repel(data=. %>% filter(!is.na(outlier)), aes(label=genes)) +
		scale_colour_gradient(low="blue",high="red") +
		geom_hline(yintercept = 0,color="red") 	
}