library(data.table)
library(dplyr)

###############
## READ DATA ##
###############

Max <- fread("~/Downloads/MaxData.tsv")
Sum <- fread("~/Downloads/SumData.tsv")
ScoreMetadata <- fread("~/Git_Repos/metabolomics_spectral_similarity_score/Metadata/Score_Metadata.txt")

##############################
## MAX: PEARSON CORRELATION ##
##############################

MaxPearson <- do.call(rbind, lapply(ScoreMetadata$Score, function(Score1) {
  
  message(Score1)
  
  do.call(rbind, lapply(ScoreMetadata$Score, function(Score2) {
    
    if (!is.numeric(Max[[Score1]])) {Max[[Score1]] <- as.numeric(Max[[Score1]])}
    if (!is.numeric(Max[[Score2]])) {Max[[Score2]] <- as.numeric(Max[[Score2]])}
    
    cor(Max[[Score1]], Max[[Score2]], method = "pearson")
    
  }))
}))

data.frame(
  Score1 = rep(ScoreMetadata$Score, each = 74),
  Score2 = rep(ScoreMetadata$Score, 74),
  Pearson = MaxPearson
) %>% fwrite("~/Git_Repos/metabolomics_spectral_similarity_score/Data/Max_Pearson.txt", quote = F, row.names = F, sep = "\t")

##############################
## MAX: SPEARMAN CORRELATION ##
##############################


