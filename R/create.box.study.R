create.box.study <- function(total.sub, gene, type, hist){
  
  # subset by gene and histology
  for.box <- total.sub[which(total.sub$Var1 == gene & total.sub$disease %in% hist),]
  
  # calculate log FPKM and z-scores
  for.box <- for.box %>% group_by(study_id, disease) %>%
    mutate(log2fpkm = round(log2(value + 1), digits = 2)) %>%
    mutate(z.score = round((log2fpkm-mean(log2fpkm))/sd(log2fpkm), digits = 2)) %>%
    mutate(median.fpkm = round(median(value), digits = 2),
           median.log2fpkm = round(median(log2fpkm), digits = 2),
           median.zscore = round(median(z.score), digits = 2)) %>%
    as.data.frame()
  
  if(type == "logFPKM"){
    val <- "log2fpkm"
    tt <- unique(as.character(for.box[order(for.box$study_id, for.box$median.log2fpkm, decreasing = F),'label']))
    for.box$label <- factor(for.box$label, levels = tt)
  } else if(type == "FPKM"){
    val <- "value"
    tt <- unique(as.character(for.box[order(for.box$study_id, for.box$median.fpkm, decreasing = F),'label']))
    for.box$label <- factor(for.box$label, levels = tt)
  } else if(type == "z-score"){
    val <- "z.score"
    tt <- unique(as.character(for.box[order(for.box$study_id, for.box$median.zscore, decreasing = F),'label']))
    for.box$label <- factor(for.box$label, levels = tt)
  }
  
  # plot
  p <- ggplot(for.box, aes_string(x = 'label', y = val, fill = 'study_id')) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    theme_bw() + theme_Publication() + 
    xlab('') + ylab(type) +
    ggtitle("Expression Boxplot")
  
  p <- plotly_build(p)
  return(p)
  
}
