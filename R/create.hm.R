create.hm <- function(dataset, total.sub, genelist.ct, sig, hist, type){
  
  # subset to gene signature
  genes <- unique(dataset[which(dataset$Signature %in% sig),'newName'])
  for.heatmap <- total.sub[which(total.sub$gene %in% genes),]
  
  # subset to histology
  for.heatmap <- for.heatmap[which(for.heatmap$disease %in% hist),]
  
  # summarise data for heatmap
  for.heatmap <- for.heatmap %>%
    mutate(log2fpkm = round(log2(fpkm + 1), digits = 2)) %>%
    mutate(z.score = round((log2fpkm-mean(log2fpkm))/sd(log2fpkm), digits = 2)) %>%
    group_by(study_id, disease, disease_name, label, gene) %>%
    dplyr::summarise(mean.fpkm = round(mean(fpkm), digits = 2), 
                     median.fpkm = round(median(fpkm), digits = 2),
                     median.log2fpkm = round(median(log2fpkm), digits = 2),
                     mean.zscore = round(mean(z.score), digits = 2),
                     median.zscore = round(median(z.score), digits = 2),
                     zscore.25 = round(quantile(z.score, probs = 0.25), digits = 2),
                     zscore.75 = round(quantile(z.score, probs = 0.75), digits = 2),
                     zscore.95 = round(quantile(z.score, probs = 0.95), digits = 2)) %>%
    unique() %>%
    as.data.frame()
  
  # create matrix
  mat <- dcast(for.heatmap, gene~label, value.var = type)
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  
  # add annotation
  annot <- unique(for.heatmap[,c('label','study_id')])
  rownames(annot) <- annot$label
  annot <- annot[order(annot$study_id, annot$label),]
  mat <- mat[,match(annot$label, colnames(mat))]
  annot$label <- NULL
  
  # create heatmap
  if(type == "mean.fpkm"){
    # pheatmap(mat = mat, annotation_col = annot, main = sig, cellwidth = 8, cellheight = 8, scale = "row", angle_col = 90)
    p <- heatmaply(x = mat, main = sig, col_side_colors = annot, plot_method = 'plotly', 
                   scale = "row", column_text_angle = -90, k_col = 4, fontsize_row = 14, fontsize_col = 14)
  } else {
    # pheatmap(mat = mat, annotation_col = annot, main = sig, cellwidth = 8, cellheight = 8, angle_col = 90)
    p <- heatmaply(x = mat, main = sig, col_side_colors = annot, plot_method = 'plotly',
                   column_text_angle = -90, k_col = 4, fontsize_row = 14, fontsize_col = 14)
  }
  
  # format data for display
  name <- sig
  n <- genelist.ct[which(genelist.ct$x %in% name),'freq']
  for.heatmap$signature <- paste0(name,' (n = ',n,')')
  for.heatmap <- for.heatmap[,c('signature','gene','study_id','disease','disease_name','label','mean.fpkm','mean.zscore','median.zscore','zscore.25','zscore.75','zscore.95')]
  
  newList <- list(p, for.heatmap)
}
