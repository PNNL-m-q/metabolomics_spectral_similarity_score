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
ScoreData_Max <- lapply(unlist(ScoreMetadata$Score), function(x) {extractScore(MAX, x)})

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
    round(sum(overlap$Overlap, na.rm = T) / nrow(overlap), 8)
  )

}

# 7. Add overlap score to ScoreMetadata
OS <- lapply(ScoreData_Sum, function(x) {calc_overlap(x)}) %>% unlist()
OS2 <- lapply(ScoreData_Max, function(x) {calc_overlap(x)}) %>% unlist()
ScoreMetadata$Overlap_Sum <- OS
ScoreMetadata$Overlap_Max <- OS2
fwrite(ScoreMetadata, "~/Git_Repos/metabolomics_spectral_similarity_score/Metadata/Score_Metadata.txt", row.names = F, quote = F, sep = "\t")

# 8. SUM: Make and hold plots 
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

# 9. SUM: Make trelliscope display      
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
        `Sum TStatistic` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(TStatisticSum) %>% unlist() %>% as.numeric() %>% round(8), desc = "TStatistic_Sum"),
        `Sum Overlap` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(Overlap_Sum) %>% unlist() %>% as.numeric() %>% round(8), desc = "Overlap_Sum"),
        `Max TStatistic` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(TStatisticMax) %>% unlist() %>% as.numeric() %>% round(8), desc = "TStatistic_Max"),
        `Max Overlap` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(Overlap_Max) %>% unlist() %>% as.numeric() %>% round(8), desc = "Overlap_Max")
      )})
    ) %>%
  dplyr::select(Score, panel, cogs) %>%
  trelliscope(name = "Sum", path = "~/Git_Repos/metabolomics_spectral_similarity_score/Trelliscopes/SS_Scores/", thumb = T)

# 10. MAX: Make and hold plots 
lapply(1:length(ScoreData_Max), function(n) {
  message(ScoreMetadata$Score[n])
  data <- ScoreData_Max[[n]]
  colnames(data) <- c("Score", "Truth Annotation")
  plot <- ggplot(data, aes(x = `Truth Annotation`, fill = `Truth Annotation`, y = Score)) +
    geom_boxplot() + 
    theme_bw() + 
    ylab(ifelse(attr(data, "Transformed") == "Transformed", "Transformed Score", "Score")) +
    ggtitle(ScoreMetadata$Score[n])
  ggsave(file.path("~/Downloads/HoldPlots/MAX", paste0(ScoreMetadata$Score[n], ".png")), plot = plot)
})

# 11. SUM: Make trelliscope display      
ScoreMetadata %>%      
  mutate(
    panel = map_plot(Score, function(x) {
      ggdraw() + draw_image(paste0("~/Downloads/HoldPlots/MAX/", x, ".png"))
    }),
    cogs = map_cog(Score, function(x) {
      tibble(
        Family = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(Family) %>% unlist(), desc = "Family"),
        Type = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(Type) %>% unlist(), desc = "Type"),
        `Theoretical Bounds` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(TheoreticalBounds) %>% unlist(), desc = "Theoretical Bounds"),
        Cluster = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(SumCluster) %>% unlist() %>% as.factor(), desc = "Cluster"),
        `Sum TStatistic` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(TStatisticSum) %>% unlist() %>% as.numeric() %>% round(8), desc = "TStatistic_Sum"),
        `Sum Overlap` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(Overlap_Sum) %>% unlist() %>% as.numeric() %>% round(8), desc = "Overlap_Sum"),
        `Max TStatistic` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(TStatisticMax) %>% unlist() %>% as.numeric() %>% round(8), desc = "TStatistic_Max"),
        `Max Overlap` = cog(ScoreMetadata %>% filter(Score == x) %>% dplyr::select(Overlap_Max) %>% unlist() %>% as.numeric() %>% round(8), desc = "Overlap_Max")
      )})
  ) %>%
  dplyr::select(Score, panel, cogs) %>%
  trelliscope(name = "Max", path = "~/Git_Repos/metabolomics_spectral_similarity_score/Trelliscopes/SS_Scores/", thumb = T)









