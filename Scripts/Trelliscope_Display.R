library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(trelliscopejs)
library(cowplot)
library(doParallel)
library(foreach)

# 1. Pull score metadata
ScoreMetadata <- fread("~/Git_Repos/metabolomics_spectral_similarity_score/Metadata/Score_Metadata.txt")

# 2. SUM: Load all data
SUM <- fread("~/Downloads/SumData.tsv")

# 3. MAX: Load all data 
MAX <- fread("~/Downloads/MaxData.tsv")

# 4. Write a function to extract the score for sum and max datasets
extractScore <- function(Dataset, Score) {
  
  message(Score)
  
  # Pull all the data associated with the score
  Data <- Dataset %>% dplyr::select(!!Score, Truth.Annotation)
  Data$Truth.Annotation <- gsub(".", " ", Data$Truth.Annotation, fixed = T)
  Data$Truth.Annotation[Data$Truth.Annotation == ""] <- "Unknown"
  Data$Truth.Annotation <- factor(Data$Truth.Annotation, levels = c("True Positive", "True Negative", "Unknown"))
  Data[[Score]] <- as.numeric(Data[[Score]])
  
  attr(Data, "Transformed") <- "Not transformed"
  
  # If score is greater than 1000, than take 1/log10(Score)
  if (max(Data[[Score]], na.rm = T) > 250) {
    Data[[Score]] <- lapply(Data[[Score]], function(x) {
      if (is.na(x)) {return(NA)} else if (x == 0) {return(0)} else {log10(1/x)}
    }) %>% unlist()
    attr(Data, "Transformed") <- "Transformed"
  }

  return(Data)
}

# 5. Extend score metadata to hold cases where max should be held or not
ScoreData_Sum <- lapply(unlist(ScoreMetadata$Score), function(x) {extractScore(SUM, x)})
names(ScoreData_Sum) <- unlist(ScoreMetadata$Score)

# 6. Create an overlap score function
calc_overlap <- function(data) {

  # Rename first column to make it easier
  colnames(data)[[1]] <- "Score"
  
  # Pull score
  score <- data %>% filter(Truth.Annotation %in% c("True Negative", "Unknown"))
  max <- max(score$Score, na.rm = T)
  min <- min(score$Score, na.rm = T)
  third <- quantile(score$Score, probs = 0.75, na.rm = T)
  first <- quantile(score$Score, probs = 0.25, na.rm = T)
  
  # Get overlap
  overlap <- data %>%
    dplyr::filter(Truth.Annotation == "True Positive") %>%
    dplyr::mutate(
      Overlap = Score <= third & Score >= first
    )
  
  # Calculate overlap 
  return(
    round(sum(overlap$Overlap) / nrow(overlap), 8)
  )

}

# 7. Add overlap score to ScoreMetadata
OS <- lapply(ScoreData_Sum, function(x) {calc_overlap(x)}) %>% unlist()
ScoreMetadata$Overlap <- OS
fwrite(ScoreMetadata, "~/Downloads/ScoreMetadata_Store.csv", row.names = F, quote = F)

# 8. Make and hold plots 
lapply(1:length(ScoreData_Sum), function(n) {
  message(ScoreMetadata$Score[n])
  data <- ScoreData_Sum[[n]]
  colnames(data) <- c("Score", "Truth Annotation")
  plot <- ggplot(data, aes(x = `Truth Annotation`, fill = `Truth Annotation`, y = Score)) +
    geom_boxplot() + 
    theme_bw() + 
    ylab(ifelse(attr(data, "Transformed") == "Transformed", "Transformed Score", "Score")) +
    ggtitle(ScoreMetadata$Score[n])
  ggsave(file.path("~/Downloads/HoldPlots/SUM", paste0(ScoreMetadata$Score[n], ".png")), plot = plot)
})

# 9. Make trelliscope display      
ScoreMetadata %>%      
  mutate(
    panel = map_plot(Score, function(x) {
      ggdraw() + draw_image(paste0("~/Downloads/HoldPlots/SUM/", x, ".png"))
    }),
    cogs = map_cog(Score, function(x) {
      tibble(
        Family = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(Family) %>% unlist(), desc = "Family"),
        Type = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(Type) %>% unlist(), desc = "Type"),
        `Theoretical Bounds` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(TheoreticalBounds) %>% unlist(), desc = "Theoretical Bounds"),
        Cluster = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(SumCluster) %>% unlist() %>% as.factor(), desc = "Cluster"),
        TStatistic = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(TStatisticSum) %>% unlist() %>% as.numeric() %>% round(8), desc = "TStatistic"),
        Overlap = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(Overlap) %>% unlist() %>% as.numeric() %>% round(8), desc = "Overlap")
      )})
    ) %>%
  dplyr::select(Score, panel, cogs) %>%
  trelliscope(name = "Sum", path = "~/Git_Repos/metabolomics_spectral_similarity_score/Trelliscopes/SS_Scores/")
  



registerDoParallel(detectCores())
foreach(n = 1:length(ScoreData_Sum)) %dopar% {
  data <- ScoreData_Sum[[n]]
  plot <- ggplot(data, aes(x = `Truth.Annotation`, y = data[[1]], fill = `Truth.Annotation`)) +
    geom_boxplot() + theme_bw() + theme(legend.position = "none") + xlab("Truth Annotation") + 
    ggtitle(attr(data, "Transformed")) + ylab(colnames(data)[1])
  ggsave(paste0("~/Desktop/HoldPlots/Sum/", colnames(data)[[1]], ".png"), plot)
}

fix1 <- ScoreData[[50]]
fix2[fix2$VW1 > 50, "VW1"] <- log10(1 / (fix2[fix2$VW1 > 50, "VW1"] %>% unlist()))
fix2 <- fix2[!is.infinite(fix2$VW1),]
fix2 <- ScoreData[[53]]

plot <- ggplot(fix2, aes(x = `Truth.Annotation`, y = VW1, fill = `Truth.Annotation`)) +
  geom_boxplot() + theme_bw() + theme(legend.position = "none") + xlab("Truth Annotation") + 
  ggtitle(attr(fix1, "Transformed")) + ylab(colnames(fix2)[1])
ggsave(paste0("~/Desktop/HoldPlots/Sum/", colnames(fix2)[1], ".png"), plot)


# Create sum trelliscope
ScoreMeta %>%
  mutate(
    panel = map_plot(Scores, function(x) {
      ggdraw() + draw_image(paste0("~/Desktop/HoldPlots/Sum/", x, ".png"))
    })
  ) %>%
  trelliscope(name = "Sum", path = "~/Downloads/SS_Scores/")

## MAX plots

extractScoreMAX <- function(Score) {
  
  # Pull all the data associated with the score
  Data <- do.call(bind_rows, lapply(MAX, function(x) {
    Res <- x %>% dplyr::select(!!Score, Truth.Annotation)
    Res[[Score]] <- as.numeric(Res[[Score]])
    return(Res)
  }))
  Data$Truth.Annotation <- gsub(".", " ", Data$Truth.Annotation, fixed = T)
  Data$Truth.Annotation[Data$Truth.Annotation == ""] <- "Unknown"
  Data$Truth.Annotation <- factor(Data$Truth.Annotation, levels = c("True Positive", "True Negative", "Unknown"))
  
  attr(Data, "Transformed") <- "Not transformed"
  
  # If score is greater than 1000, than take 1/log10(Score)
  if (max(Data[[Score]], na.rm = T) > 250) {
    Data[[Score]] <- lapply(Data[[Score]], function(x) {
      if (is.na(x)) {return(NA)} else if (x == 0) {return(0)} else {log10(1/x)}
    }) %>% unlist()
    attr(Data, "Transformed") <- "Transformed"
  }
  
  return(Data)
}

ScoreDataMAX <- lapply(unlist(ScoreMeta$Scores), function(x) {message(x); extractScoreMAX(x)})

registerDoParallel(detectCores())
foreach(n = 1:length(ScoreDataMAX)) %dopar% {
  data <- ScoreDataMAX[[n]]
  plot <- ggplot(data, aes(x = `Truth.Annotation`, y = data[[1]], fill = `Truth.Annotation`)) +
    geom_boxplot() + theme_bw() + theme(legend.position = "none") + xlab("Truth Annotation") + 
    ggtitle(attr(data, "Transformed")) + ylab(colnames(data)[1])
  ggsave(paste0("~/Desktop/HoldPlots/Max/", colnames(data)[[1]], ".png"), plot)
}

fix3 <- ScoreDataMAX[[50]]
fix3[fix3$`Symmetric Chi Squared Distance` > 50, "Symmetric Chi Squared Distance"] <- log10(1 / (fix3[fix3$`Symmetric Chi Squared Distance` > 50, "Symmetric Chi Squared Distance"] %>% unlist()))
fix3 <- fix3[!is.infinite(fix3$`Symmetric Chi Squared Distance`),]

plot <- ggplot(fix3, aes(x = `Truth.Annotation`, y = `Symmetric Chi Squared Distance`, fill = `Truth.Annotation`)) +
  geom_boxplot() + theme_bw() + theme(legend.position = "none") + xlab("Truth Annotation") + 
  ggtitle(attr(fix1, "Transformed")) + ylab(colnames(fix3)[1])
ggsave(paste0("~/Desktop/HoldPlots/Max/", colnames(fix3)[1], ".png"), plot)


# Remove plots that aren't supposed to be included with max 
toRemove <- ScoreMeta[!ScoreMeta$IncludeMA, "Scores"] %>% unlist()
toRemove <- toRemove[toRemove != "Squared-chord Distance"]
lapply(toRemove, function(score) {
  ggsave(paste0("~/Desktop/HoldPlots/Max/", score, ".png"), ggplot() + ggtitle("Score Removed") + theme_bw())
})

ScoreMeta %>%
  unique() %>%
  mutate(
    panel = map_plot(Scores, function(x) {
      ggdraw() + draw_image(paste0("~/Desktop/HoldPlots/Max/", x, ".png"))
    })
  ) %>%
  trelliscope(name = "Max", path = "~/Downloads/SS_Scores/")


