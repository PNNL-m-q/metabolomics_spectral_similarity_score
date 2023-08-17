library(data.table)
library(dplyr)
library(randomForest)

ScoreMetadata <- fread("~/Git_Repos/metabolomics_spectral_similarity_score/Metadata/Score_Metadata.txt")
Sum <- fread("~/Downloads/SumData.tsv")

# Make sure all columns are numeric
data.frame(
  Score = ScoreMetadata$Score,
  Class = lapply(ScoreMetadata$Score, function(x) {class(Sum[[x]])}) %>% unlist()
) %>%
  filter(Class != "numeric")
Sum$`Kumar Johnson Distance` <- as.numeric(Sum$`Kumar Johnson Distance`)
Sum$`Kumar Johnson Divergence` <- as.numeric(Sum$`Kumar Johnson Divergence`)

# Z Score normalize
for (x in ScoreMetadata$Score) {
  message(x)
  Sum[[x]] <- (Sum[[x]] - mean(Sum[[x]], na.rm = T)) / sd(Sum[[x]], na.rm = T)
}

# Make a sample size
n_samps <- 1000

# Subset down to true positives 
TPs <- Sum %>% filter(Truth.Annotation == "True.Positive")
TPs <- TPs[sample(1:nrow(TPs), n_samps)]

# Subset down to true negatives and randomly sample 
TNs <- Sum %>% filter(Truth.Annotation == "True.Negative")
TNs <- TNs[sample(1:nrow(TNs), n_samps)]

toModel <- rbind(TPs, TNs) %>%
  dplyr::select(Truth.Annotation, ScoreMetadata$Score) %>%
  dplyr::mutate(Truth.Annotation = factor(ifelse(Truth.Annotation == "True.Positive", "1", "0"), levels = c("0", "1"))) %>%
  data.frame()

# As determined by the other model
RF <- randomForest(Truth.Annotation ~ ., data = toModel, importance = T, trees = 1562, mtry = 15)



