create.box.gsva <- function(dat, meta, genelist.ct, hist, sig, type){
  
  # get sample ids
  s <- meta[which(meta$disease %in% hist),'sample_id']
  
  # subset to gene signature and samples per histology
  rownames(dat) <- trimws(rownames(dat))
  dat <- dat[rownames(dat) == sig, colnames(dat) %in% s]
  
  # add metadata
  dat <- melt(dat, value.name = "gsva.score")
  dat$sample_id <- rownames(dat)
  dat <- merge(dat, meta[,c('sample_id', 'study_id', 'disease', 'disease_name', 'label')], by = 'sample_id')
  
  # plot gsva scores
  name <- sig
  n <- genelist.ct[which(genelist.ct$x %in% name),'freq']
  dat <- dat %>% group_by(study_id, disease, disease_name, label) %>% 
    mutate(gsva.score = round(gsva.score, digits = 2)) %>%
    mutate(median = round(median(gsva.score), digits = 2)) %>%
    as.data.frame()
  tt <- unique(as.character(dat[order(dat$study_id, dat$median, decreasing = F),'label']))
  dat$label <- factor(dat$label, levels = tt)
  p <- ggplot(dat, aes(x = label, y = gsva.score, alpha = 0.5)) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = NA, width = 0.5, aes(fill = study_id)) +
    theme_bw() + theme_Publication() + 
    xlab('') + ylab('GSVA scores') +
    ggtitle(paste0("GSVA scores: ", name, "\n(Geneset size = ",n,")")) + guides(alpha = FALSE)
  p <- plotly_build(p)
  
  # format data for display
  dat$signature <- paste0(name,' (n = ',n,')')
  dat <- dat[,c('signature','sample_id','study_id','disease','disease_name','label','gsva.score')]
  
  newList <- list(p, dat)
  return(newList)
  
}