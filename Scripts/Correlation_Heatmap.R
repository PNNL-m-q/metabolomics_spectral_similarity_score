library(data.table)
library(dplyr)

###############
## READ DATA ##
###############

Max <- fread("~/Downloads/MaxData.tsv")
Sum <- fread("~/Downloads/SumData.tsv")
ScoreMetadata <- fread("~/Git_Repos/metabolomics_spectral_similarity_score/Metadata/Score_Metadata.txt")

# Remove the duplicates
# Jensen Differences Divergence  --> Jensen Difference Distance
# Vicis Wave Hadges Distance --> VW1 
# Vicis Symmetric Chi Squared 1 Distance --> VW2 
# Vicis Symmetric Chi Squared 2 Distance --> VW3
# Vicis Symmetric Chi Squared 3 Distance --> VW4
# Max Symmetric Chi Squared Distance --> VW5
# Min Symmetric Chi Squared Distance --> VW6
# Spectral Similarity Score --> Cosine Correlation

ScoreMetadata <- ScoreMetadata %>%
  filter(Score %in% c("Vicis Wave Hadges Distance", "Vicis-Symmetric X2 1 Distance",
                      "Vicis-Symmetric X2 2 Distance", "Vicis-Symmetric X2 3 Distance",
                      "Max Symmetric Chi Squared Distance", "Min Symmetric Chi Squared Distance",
                      "Spectral Similarity Score", "Kumar Johnson Distance") == FALSE)


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
  Score1 = rep(ScoreMetadata$Score, each = 66),
  Score2 = rep(ScoreMetadata$Score, 66),
  Pearson = MaxPearson
) %>% fwrite("~/Git_Repos/metabolomics_spectral_similarity_score/Data/Max_Pearson.txt", quote = F, row.names = F, sep = "\t")

##############################
## SUM: PEARSON CORRELATION ##
##############################

SumPearson <- do.call(rbind, lapply(ScoreMetadata$Score, function(Score1) {
  
  message(Score1)
  
  do.call(rbind, lapply(ScoreMetadata$Score, function(Score2) {
    
    if (!is.numeric(Sum[[Score1]])) {Sum[[Score1]] <- as.numeric(Sum[[Score1]])}
    if (!is.numeric(Sum[[Score2]])) {Sum[[Score2]] <- as.numeric(Sum[[Score2]])}
    
    cor(Sum[[Score1]], Sum[[Score2]], method = "pearson")
    
  }))
}))

data.frame(
  Score1 = rep(ScoreMetadata$Score, each = 66),
  Score2 = rep(ScoreMetadata$Score, 66),
  Pearson = SumPearson
) %>% fwrite("~/Git_Repos/metabolomics_spectral_similarity_score/Data/Sum_Pearson.txt", quote = F, row.names = F, sep = "\t")



