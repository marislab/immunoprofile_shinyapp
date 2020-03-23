create.expr.hm <- function(dataset, total.sub, genelist.ct, sig, hist, type){
  
  # subset to gene signature
  genes <- unique(dataset[which(dataset$Signature %in% sig),'newName'])
  for.heatmap <- total.sub[which(total.sub$gene %in% genes),]
  
  # subset to histology
  for.heatmap <- for.heatmap[which(for.heatmap$disease %in% hist),]
  
  # summarise data for heatmap
  for.heatmap <- for.heatmap %>%
    mutate(log2fpkm = round(log2(fpkm + 1), digits = 2)) %>%
    mutate(z.score = round((log2fpkm-mean(log2fpkm))/sd(log2fpkm), digits = 2)) %>%
    unique() %>%
    as.data.frame()
  
  # create matrix
  mat <- dcast(for.heatmap, gene~sample_id, value.var = type)
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  
  # add annotation
  # annot <- unique(for.heatmap[,c('label','study_id')])
  # rownames(annot) <- annot$label
  # annot <- annot[order(annot$study_id, annot$label),]
  # mat <- mat[,match(annot$label, colnames(mat))]
  # annot$label <- NULL
  
  # get info on clusters
  # annotation for column
  dend <- hclust(dist(t(mat), method = "euclidean"), method = "complete")
  annot <- cutree(dend, k = 3)
  dend <- melt(annot, value.name = "cluster")
  for.heatmap <- merge(for.heatmap, dend, by.x = "sample_id", by.y = "row.names")
  
  # create heatmap
  p <- heatmaply(x = mat, main = sig, plot_method = 'plotly', showticklabels = c(FALSE, TRUE),
                 column_text_angle = -90, 
                 col_side_colors = annot,
                 k_col = 3, dist_method = "euclidean", hclust_method ="complete",
                 fontsize_row = 14, fontsize_col = 14)
  
  # create boxplot
  for.box <- for.heatmap %>%
    group_by(sample_id, cluster) %>%
    summarise(fpkm = mean(fpkm),
              z.score = mean(z.score),
              log2fpkm = mean(log2fpkm)) %>%
    mutate(cluster = paste0("Cluster ",as.character(cluster)))
  my_comparisons <- list(c("1", "2"), c("2", "3"), c("1", "3"))
  q <- ggplot(for.box, aes_string(x = 'cluster', y = type)) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = NA, width = 0.5, aes(fill = cluster)) +
    theme_bw() + theme_Publication2() + 
    xlab('') + ylab(type) + guides(alpha = FALSE, fill = FALSE) +
    stat_compare_means(comparisons = my_comparisons) +
    stat_compare_means(method = "anova")
  q <- plotly_build(q)
  
  # format data for display
  name <- sig
  n <- genelist.ct[which(genelist.ct$x %in% name),'freq']
  for.heatmap$signature <- paste0(name,' (n = ',n,')')
  for.heatmap <- for.heatmap[,c('sample_id','cluster','signature','gene','study_id','disease','disease_name','label','fpkm','log2fpkm','z.score')]
  
  newList <- list(p, q, for.heatmap)
  return(newList)
}