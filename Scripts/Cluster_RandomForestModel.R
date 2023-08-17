library(data.table)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(patchwork)
library(DT)
library(rsample)
library(randomForest)
library(recipes)
library(ranger)

set.seed(5373)

#' @returns A list of the following objects:
#' 1. The best fitted random forest
#' 2. The model fit to the testing dataset
#' 3. The model accuracy
#' 4. Variable importance
run_randomforest <- function(theDataset) {
  
  # Split data into testing and training 
  counts <- theDataset$Cluster %>% 
    table(dnn = "Cluster") %>% 
    data.frame() %>%
    mutate(
      TestSize = floor(Freq*0.25)
    )
  testRows <- c(sample(which(theDataset$Cluster == 1), counts$TestSize[counts$Cluster == 1]),
                sample(which(theDataset$Cluster == 2), counts$TestSize[counts$Cluster == 2]),
                sample(which(theDataset$Cluster == 3), counts$TestSize[counts$Cluster == 3]),
                sample(which(theDataset$Cluster == 4), counts$TestSize[counts$Cluster == 4]))
  test <- theDataset[testRows,]
  train <- theDataset[1:nrow(theDataset) %in% testRows == FALSE,]
  
  # Create a recipe 
  recipe <- recipes::recipe(Cluster ~ ., data = train) %>%
    recipes::step_normalize(Min, Median, Max)
  
  # Designate a resampling scheme 
  folds <- rsample::vfold_cv(train,
                             v = 10,
                             repeats = 5,
                             strata = Cluster,
                             breaks = 4)
  
  # Designate modeling engine
  rf_engine <- parsnip::rand_forest(trees = tune::tune(), mtry = tune::tune()) %>%
    parsnip::set_engine("ranger", importance = "impurity") %>%
    parsnip::set_mode("classification")
  
  # Designate test metrics 
  test_metrics <- yardstick::metric_set(yardstick::accuracy)
  
  # Specify number of trees based on a number of levels and m-try
  rf_grid <- dials::grid_regular(dials::trees(), dials::mtry(range = c(1, ncol(theDataset)-1)), levels = 5) 
  
  # Run the workflow - currently 1 random forest model
  wf <- workflowsets::workflow_set(
    preproc = list("nzv" = recipe),
    models = list(rf = rf_engine)
  )
  wfid_names <- wf$wflow_id
  train_wfs <- wf %>% 
    workflowsets::option_add(id = wfid_names[which(grepl("rf", wfid_names))], grid = rf_grid)
  
  # Designate half of the cores to run 
  ncores <- floor(parallel::detectCores(logical = TRUE)/2)
  
  # Set the number of cores
  cl <- parallel::makePSOCKcluster(ncores)
  
  # Register cores to use
  doParallel::registerDoParallel(cl)
  
  # Run and resolve workflows 
  wf_complete <- wf %>%
    workflowsets::workflow_map(resamples = folds,
                               verbose = TRUE, 
                               seed = 339,
                               metrics = test_metrics)
  
  # Stop parallel computing cores 
  parallel::stopCluster(cl)
  
  # extract best fit parameters
  best_rf_sum <- wf_complete %>%
    workflowsets::extract_workflow_set_result(id = "nzv_rf") %>%
    tune::select_best(metric = "accuracy")
  
  # fit to training
  fitted_randf_all <- wf_complete %>%
    workflowsets::extract_workflow(id = "nzv_rf") %>%
    tune::finalize_workflow(best_rf_sum) %>%
    parsnip::fit(data = train)
  
  # fit to testing
  cnmc_all <- predict(fitted_randf_all, test) %>%
    dplyr::bind_cols(predict(fitted_randf_all, test, type = "prob")) %>%
    dplyr::bind_cols(test %>% dplyr::select(Cluster)) %>%
    dplyr::rename(pred_class = .pred_class)
  
  # model accuracy
  cnmc_all_accuracy <- cnmc_all %>%
    test_metrics(truth = Cluster, estimate = pred_class)
  
  # Variable importance
  importance <- fitted_randf_all %>%
    workflows::extract_fit_parsnip() %>%
    vip::vi(scale = TRUE) %>%
    dplyr::arrange(dplyr::desc(Importance)) 
  
  return(list(
    best_rf_sum,
    cnmc_all,
    cnmc_all_accuracy,
    importance
    
  ))
  
}

ScoreMetadata <- fread("~/Git_Repos/metabolomics_spectral_similarity_score/Metadata/Score_Metadata.txt")

# Extract only Sum Relevant Metadata
SumMetadata <- ScoreMetadata %>% 
  dplyr::select(SumCluster, Family, Type, TheoreticalBounds, Sum_ScoreMin, Sum_ScoreMedian, Sum_ScoreMax, TStatisticSum, Overlap_Sum) %>%
  dplyr::mutate(
    SumCluster = factor(SumCluster, levels = 1:4),
    Family = as.factor(Family),
    Type = as.factor(Type), 
    TheoreticalBounds = as.factor(TheoreticalBounds)
  ) %>% 
  dplyr::rename(Min = Sum_ScoreMin, Median = Sum_ScoreMedian, Max = Sum_ScoreMax, Cluster = SumCluster)

Sum_RF_Results <- run_randomforest(SumMetadata)
saveRDS(Sum_RF_Results, "~/Git_Repos/metabolomics_spectral_similarity_score/Data/Sum_RF_Results.RDS")

# Extract only Max Relevant Metadata
MaxMetadata <- ScoreMetadata %>% 
  dplyr::select(MaxCluster, Family, Type, TheoreticalBounds, Max_ScoreMin, Max_ScoreMedian, Max_ScoreMax, TStatisticMax, Overlap_Max) %>%
  dplyr::mutate(
    MaxCluster = factor(MaxCluster, levels = 1:4),
    Family = as.factor(Family),
    Type = as.factor(Type), 
    TheoreticalBounds = as.factor(TheoreticalBounds)
  ) %>% 
  dplyr::rename(Min = Max_ScoreMin, Median = Max_ScoreMedian, Max = Max_ScoreMax, Cluster = MaxCluster) %>%
  dplyr::mutate(
    Min = as.numeric(Min),
    Median = as.numeric(Median),
    Max = as.numeric(Max)
  )

Max_RF_Results <- run_randomforest(MaxMetadata)
saveRDS(Max_RF_Results, "~/Git_Repos/metabolomics_spectral_similarity_score/Data/Max_RF_Results.RDS")

