# just for NBL currently
# this script will do a diff expr analysis on Inflamed vs Desert samples

# group1 = "Inflamed"
# group2 = "Desert"
# study <- "TARGET"
# hist <- "NBL"
# query <- dei.res[[2]]

dei.diffexpr <- function(query, study, hist, group1 = "Inflamed", group2 = "Desert"){
  
  fname <- file.path('data', 'teffector_diffexpr', paste0(study, '_', hist, '.txt'))
  if(!file.exists(fname)){
    # counts file for all pediatric datasets
    if(!exists('hist.csv')){
      hist.csv <- 'data/pediatric_datasets_counts.RDS'
      hist.csv <<- readRDS(hist.csv)
    }
    
    # get design matrix
    query.sub <- query %>%
      filter(Type %in% c(group1, group2)) %>%
      as.data.frame()
    rownames(query.sub) <- query.sub$sample_id
    
    # subset to query samples
    hist.csv <- hist.csv[,which(colnames(hist.csv) %in% rownames(query.sub))]
    query.sub  <- query.sub[colnames(hist.csv),]
    
    # remove 0 counts
    hist.csv <- hist.csv[apply(hist.csv, 1, function(x) !all(x==0)),]
    
    # voom + limma
    if(identical(rownames(query.sub), colnames(hist.csv))){
      print("Pass")
      query.sub$Type <- factor(query.sub$Type)
      design.matrix <- model.matrix(~0+query.sub$Type)
      rownames(design.matrix) <- colnames(hist.csv)
      colnames(design.matrix) <- levels(query.sub$Type)
      tmp.voom <- voom(counts = as.matrix(hist.csv), design = design.matrix)
      tmp.voom <- tmp.voom$E
      
      fit <- lmFit(object = tmp.voom, design = design.matrix)
      contrast.matrix <- makeContrasts(contrasts = c("Inflamed-Desert"), levels = design.matrix)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      res <- topTable(fit2, number = Inf, p.value = 0.01)
      res <- res %>%
        rownames_to_column("gene_symbol") %>%
        filter(abs(logFC) > 1)  %>%
        select(gene_symbol, logFC, adj.P.Val)
        
      write.table(res, file = fname, quote = F, sep = "\t", row.names = F)
    } else {
      print("Fail")
    }
  } else {
    # read in already generated diff expression result
    res <- read.delim(fname, stringsAsFactors = F)
    return(res)
  }
}