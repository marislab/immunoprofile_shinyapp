# check code andboxplots etc for each 
setwd('~/Projects/Maris-lab/immunoprofile_shinyapp/')
load('data/filtered_expr_with_meta.RData')
load('data/histology_labs.RData')
gsva.zscore <- readRDS('data/Immunoprofiling_signatures.RDS')
gsva.log2counts <- readRDS('data/Immunoprofiling_signatures_log2counts.RDS')
dataset <- read.delim('data/ImmuneSig_GeneList_OBD_final_2_19_2019.txt', stringsAsFactors = F)
dataset$Signature <- trimws(dataset$Signature)
genelist.ct <- plyr::count(dataset$Signature)
genes <- unique(dataset$newName)
sigs <- unique(dataset$Signature)

# histology breakdown by study
load('data/histology_labs.RData')
studies <- unique(lb$study_id)

# gene
gene <- 'MYCN'
type <- 'FPKM'

# gsva
meta <- x
dat <- gsva.zscore
genelist.ct <- genelist.ct
hist <- lb$disease
sig <- sigs[1]
type <- 'z-scores'

# heatmap
type <- 'mean.fpkm'

# sunburst
load('data/histology_labs.RData')
sb.dat <- plyr::count(x[,c('study_id','disease_name','definition')])
colnames(sb.dat) <- c("level1","level2","level3","size")
sn <- sunburst(
  d3_nest(sb.dat, value_cols = "size"),
  count = T, percent = F, 
  legend = F, withD3 = T)
sn

# pie charts
sig.plot <- dataset %>% 
  group_by(Signature) %>% 
  summarise(count = n(), newName = toString(newName)) %>%
  as.data.frame()
sig.plot$Signature <- reorder(sig.plot$Signature, sig.plot$count)
p <- plot_ly(data = sig.plot, 
             x = ~Signature,
             y = ~count, 
             text = ~newName,
             type = "bar") %>%
  layout(xaxis = list(title = "", tickangle = -90),
         yaxis = list(title = "Number of Genes"),
         title = "Gene Signatures")
p
