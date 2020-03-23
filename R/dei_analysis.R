# function to do Teffector analysis
# study  <- "TARGET"
# hist <- "NBL"
genes <- c("CD8A", "GZMA", "GZMB", "PRF1", "CXCL9", "CXCL10", "TBX21")
dei <- function(total.sub, study, hist){
  if(!file.exists('data/TCGA_upper_lower_quartiles.txt')){
    tcga <- total.sub %>%
      filter(study_id  == "TCGA",
             gene %in% genes)
    
    # tcga quartile
    q <- tcga %>%
      group_by(label) %>%
      summarise(inflamed_75 = round(quantile(log2(fpkm + 1), probs = 0.75,  names = F), 2),
                desert_25 = round(quantile(log2(fpkm + 1), probs = 0.25, names = F), 2))
    
    # write out
    write.table(q, file = "data/TCGA_upper_lower_quartiles.txt", quote = F, row.names = F, sep = "\t")
  }  else {
    # read quartile information from TCGA
    q <- read.delim('data/TCGA_upper_lower_quartiles.txt', stringsAsFactors = F)
  }
  
  # for the study and histology in question, do the analysis
  query <- total.sub %>%
    filter(study_id %in% study,
           disease %in% hist,
           gene %in% genes) %>%
    group_by(label, sample_id) %>%
    summarise(avg.fpkm = mean(fpkm))
  
  # assign colors to bars
  query$Type <- ifelse(query$avg.fpkm > mean(q$inflamed_75), "Inflamed",
                            ifelse(query$avg.fpkm < mean(q$desert_25), "Desert", "Excluded"))
  
  # barplot
  query$sample_id <- reorder(query$sample_id, query$avg.fpkm)
  p <- ggplot(query, aes(x = sample_id, y = log2(avg.fpkm + 1), fill = Type)) + 
    geom_bar(stat = "identity") +
    theme_bw() +
    xlab("Samples") + ylab("Average Expression in FPKM (CD8+ Teffector Signature)") +
    theme(axis.text.x = element_blank()) +
    ggtitle(paste0("T effector analysis: ", hist)) +
    geom_hline(yintercept = log2(mean(q$inflamed_75) + 1), linetype = "dashed", color = "red") +
    geom_hline(yintercept = log2(mean(q$desert_25) + 1), linetype = "dashed", color = "blue") +
    scale_fill_manual(values = c("Desert" = "#FFCC00",
                                 "Inflamed" = "#FF0000",
                                 "Excluded" = "#FF9900"))
  p <- plotly_build(p)
  newList <- list(p, query)
  return(newList)
}
