# pathway analysis on inflamed vs desert genes results
# diffexpr.res <- read.delim('data/teffector_diffexpr/TARGET_NBL.txt')
dei.pathway <- function(diffexpr.res, pathway = "H"){
  ranks <- diffexpr.res$logFC
  names(ranks) <- diffexpr.res$gene_symbol
  m_df = msigdbr(species = "Homo sapiens", category = pathway)
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  fgseaRes <- fgsea(pathways = m_list, 
                    stats = ranks,
                    minSize=15,
                    maxSize=500,
                    nperm=10000) 
  fgseaRes <- fgseaRes %>%
    filter(padj < 0.05) %>%
    arrange(padj)
}