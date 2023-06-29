library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(trelliscopejs)
library(cowplot)
library(doParallel)
library(foreach)

## First SS analysis 
setwd("~/Desktop/MQ/Similarity_Metrics_Paper/FilesFromLisa/")

# 0. Load names of scores
Name_Maps <- fread("cname_mapping.csv")
Name_Maps$Sum_Metric_cnames[46] <- "Squared-chord Distance"

# 1. Identify similarities (0-1) and distances 
MaxEx1 <- fread("/Users/degn400/Desktop/MQ/MQ_Data/Final_Dataset/MAX/Urine_01_Max_Abundance_TruthAnnotatedFeb2022.tsv")

## Define metric types
MetricTypes <- data.table(
  "Score" = colnames(MaxEx1)[13:85],
  "Type" = lapply(colnames(MaxEx1)[13:85], function(x) {
    if (grepl("Distance", x)) {return("Distance")}
    if (max(MaxEx1[[x]], na.rm = T) <= 1)  {"Similarity"} else{"Distance"}
  }) %>% unlist()
)
  
# 2. SUM: Load all data
SUM <- lapply(list.files("~/Desktop/MQ/MQ_Data/Final_Dataset/SUM/", full.names = T), function(path) {
  message(path)
  return(fread(path))
})
saveRDS(SUM, "~/Downloads/SUM.RDS")
SUM <- readRDS("~/Downloads/SUM.RDS")

# 3. MAX: Load all data 
MAX <- lapply(list.files("~/Desktop/MQ/MQ_Data/Final_Dataset/MAX/", full.names = T), function(path) {
  message(path)
  data <- fread(path)
  data$`Vicis Wave Hadges Distance` <- as.numeric(data$`Vicis Wave Hadges Distance`)
  return(data)
})
saveRDS(MAX, "~/Downloads/MAX.RDS")
MAX <- readRDS("~/Downloads/MAX.RDS")

# 4. RAW: Load all data 
RAW <- lapply(list.files("~/Desktop/MQ/MQ_Data/Final_Dataset/RAW/", full.names = T), function(path) {
  message(path)
  return(fread(path))
})
saveRDS(RAW, "~/Downloads/RAW.RDS")
RAW <- readRDS("~/Downloads/RAW.RDS")

## Double check counts of identifications
DC <- data.frame(
  "Names" = list.files("~/Desktop/MQ/MQ_Data/Final_Dataset/RAW/") %>% unlist(),
  "SUM" = lapply(SUM, nrow) %>% unlist(), 
  "MAX" = lapply(MAX, nrow) %>% unlist(),
  "RAW" = lapply(RAW, nrow) %>% unlist()
) %>% dplyr::mutate(
  "NoFixes" = SUM == MAX & MAX == RAW
)

extractScore <- function(Score) {
  
  # Pull all the data associated with the score
  Data <- do.call(bind_rows, lapply(SUM, function(x) {
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


# Start dataframe to hold all scores 
ScoreMeta <- data.table(
  Scores = Name_Maps$Sum_Metric_cnames %>% gsub(pattern = "_", replacement = " "),
  IncludeMA = Name_Maps$Sum_Metric_cnames %in% Name_Maps$Max_Metric_cnames
) 
ScoreData <- lapply(unlist(ScoreMeta$Scores), function(x) {message(x); extractScore(x)})

# SUM plots
registerDoParallel(detectCores())
foreach(n = 1:length(ScoreData)) %dopar% {
  data <- ScoreData[[n]]
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


