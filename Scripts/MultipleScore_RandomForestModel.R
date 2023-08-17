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


seed_num <- 4656
set.seed(seed_num)


#' @returns A list of the following objects:
#' 1. The best fitted random forest
#' 2. The model fit to the testing dataset
#' 3. The model accuracy
#' 4. Variable importance
run_randomforest <- function(theDataset) {
  
  # Start with 1000 samples 
  n_samps <- 1000
  
  # NA filter data
  theDataset <- theDataset %>% na.omit()
  
  # Scale data (negative log transform of x / max(x))
  for (x in ScoreMetadata$Score) {
    theDataset[[x]] <- abs(theDataset[[x]] / max(theDataset[[x]]))
  }
  
  # Subset down to true positives 
  TPs <- theDataset %>% filter(Truth.Annotation == "True.Positive")
  TPs <- TPs[sample(1:nrow(TPs), n_samps)]
  
  # Subset down to true negatives and randomly sample 
  TNs <- theDataset %>% filter(Truth.Annotation == "True.Negative")
  TNs <- TNs[sample(1:nrow(TNs), n_samps)]
  
  # Create dataset
  toModel <- rbind(TPs, TNs) %>%
    dplyr::select(Truth.Annotation, ScoreMetadata$Score) %>%
    dplyr::mutate(Truth.Annotation = factor(ifelse(Truth.Annotation == "True.Positive", "1", "0"), levels = c("0", "1"))) %>%
    data.frame()
  
  # Split data into testing and training 
  counts <- toModel %>%
    count(Truth.Annotation) %>%
    mutate(
      TestSize = floor(n * 0.25)
    )
  testRows <- c(sample(which(toModel$Truth.Annotation == "0"), counts$TestSize[counts$Truth.Annotation == "0"]),
                sample(which(toModel$Truth.Annotation == "1"), counts$TestSize[counts$Truth.Annotation == "1"]))
  test <- toModel[testRows,]
  train <- toModel[1:nrow(toModel) %in% testRows == FALSE,]
  
  # Create a recipe 
  recipe <- recipes::recipe(Truth.Annotation ~ ., data = train)
  
  # Designate a resampling scheme 
  folds <- rsample::vfold_cv(train,
                             v = 10,
                             repeats = 5,
                             strata = Truth.Annotation,
                             breaks = 4)
  
  # Designate modeling engine
  rf_engine <- parsnip::rand_forest(trees = tune::tune(), mtry = tune::tune()) %>%
    parsnip::set_engine("ranger", importance = "impurity") %>%
    parsnip::set_mode("classification")
  
  # Designate test metrics 
  test_metrics <- yardstick::metric_set(yardstick::accuracy)
  
  # Specify number of trees based on a number of levels and m-try
  rf_grid <- dials::grid_regular(dials::trees(), dials::mtry(range = c(1, ncol(toModel)-1)), levels = 3) 
  
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
                               seed = seed_num,
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
    dplyr::bind_cols(test %>% dplyr::select(Truth.Annotation)) %>%
    dplyr::rename(pred_class = .pred_class)
  
  # model accuracy
  cnmc_all_accuracy <- cnmc_all %>%
    test_metrics(truth = Truth.Annotation, estimate = pred_class)
  
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

## Sum ##

Sum <- fread("~/Downloads/SumData.tsv")

# Make sure all columns are numeric
data.frame(
  Score = ScoreMetadata$Score,
  Class = lapply(ScoreMetadata$Score, function(x) {class(Sum[[x]])}) %>% unlist()
) %>%
  filter(Class != "numeric")
Sum$`Kumar Johnson Distance` <- as.numeric(Sum$`Kumar Johnson Distance`)
Sum$`Kumar Johnson Divergence` <- as.numeric(Sum$`Kumar Johnson Divergence`)

saveRDS(run_randomforest(Sum), "~/Downloads/Sum_MultipleScoreRF.RDS")

Results <- readRDS("~/Downloads/Sum_MultipleScoreRF.RDS")

Results[[2]] %>% dplyr::select(pred_class, Truth.Annotation) %>%
  table()




