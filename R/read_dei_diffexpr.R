# function to do limma analysis between desert and inflamed types
# this has to be done on the toil server or separately as it cannot be done on the fly
# once the results are ready, read those in
# study <- 'TARGET'
# hist <- 'NBL'
# comp <- c("Inflamed-Desert")
read.dei.diffexpr <- function(study, hist) {
  fname <- paste0('data/teffector_diffexpr/', study, '_', hist, '.txt')
  if(!file.exists(fname)){
    # do limma analysis
  } else {
    res <- read.delim(fname, stringsAsFactors = F)
    return(res)
  }
}