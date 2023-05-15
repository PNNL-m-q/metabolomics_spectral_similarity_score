library(dplyr)
library(magrittr)
library(ggplot2)
library(rules)
library(baguette)
library(doParallel)

## Data Import/Processing -------------------------------------------------------------

metfreqs <- readRDS(file = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\NMR\\FY23_Recommender\\Data\\EMSL_metfreqs_nodupes.rds")
fy22_top100 <- readRDS(file = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\NMR\\FY23_Recommender\\Data\\EMSL_fy22_top100.rds")
freqthresh <- as.numeric(summary(metfreqs$Freq)["Median"])
metfreqs <- metfreqs %>% dplyr::filter(Freq >= freqthresh) %>%
  dplyr::filter(!(rownames(.) %in% c("DSS-d6 (Chemical Shape Indicator)", "DSS (Chemical Shape Indicator)"))) # Retain only compounds with at least as many cases as the median amount and remove DSS

metfreqs %>% dplyr::filter(rownames(.) %in% fy22_top100) %>% .$Freq
fy22_top100 %in% rownames(metfreqs)

# Import concentration data
fconc <- readRDS(file = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\NMR\\FY23_Recommender\\Data\\FilteredConcentrations_20221207_nodupes.rds")

# 0.1 bw ------------------------------------------------------------------

# Import chemshift data
fbinppm1 <- data.table::fread("//PNL/Projects/NMRsoftware/FY23_Project_Files/Data/WK_datamine/Binned_Data_Feb2023/StandardisedArea/Filtered_BinnedData_0.1-bw.csv.gz", header = TRUE)
names(fbinppm1)[1] <- "CNX_id"

# Note: the processed fconc data used here contains data for the same number of samples as in fbinppm1 
setdiff(fbinppm1$CNX_id, fconc$CNX_id) # all fbinppm1 samples are in fconc
setdiff(fconc$CNX_id, fbinppm1$CNX_id) # all fconc samples are in fbinppm1

# -------------------------------------------------------------------------


# Model Fitting -----------------------------------------------------------


# Models ------------------------------------------------------------------

test_metrics <- yardstick::metric_set(yardstick::accuracy, 
                                      #yardstick::roc_auc,
                                      #yardstick::pr_auc,
                                      yardstick::f_meas,
                                      yardstick::recall,
                                      yardstick::precision,
                                      yardstick::sensitivity,
                                      yardstick::specificity)


# Normalizing by predictors - min, median, and max values 
mod_rf <- parsnip::rand_forest(trees = tune::tune()) %>%
  parsnip::set_engine("ranger") %>%
  parsnip::set_mode("classification")
rf_grid <- dials::grid_regular(dials::trees(),
                               levels = 3)


# Tuning / Fitting / Evaluation -------------------------------------------

i <- 1
met <- rownames(metfreqs)[i]

# Create outcome by dichotomizing measured concentration for specified metabolite
fconc_met <- fconc %>%
  dplyr::select(CNX_id, dplyr::all_of(met)) %>%
  dplyr::mutate(met_bin = factor(ifelse(is.na(.data[[met]]), "Absent", "Present"), levels = c("Present", "Absent"))) %>%
  dplyr::select(-dplyr::all_of(met))

# Create dataset by joining dichotomous outcome to spectral data
met_dat <- dplyr::full_join(fconc_met, fbinppm1, by = "CNX_id") %>%
  select(-CNX_id) %>%
  setNames(c("met_bin", paste0("X", names(.)[-1])))

dat_split <- rsample::initial_split(data = met_dat, prop = 3/4, strata = met_bin)
train_dat <- rsample::training(dat_split)
test_dat  <- rsample::testing(dat_split)

# Recipe creation
dat_rec <- recipes::recipe(met_bin ~ ., data = train_dat) %>%
  recipes::step_nzv(recipes::all_predictors())

# Define resampling scheme
set.seed(1)
dat_folds <- rsample::vfold_cv(train_dat,
                               v = 10,
                               repeats = 5,
                               strata = met_bin)

# Define workflow set
all_wfs <- workflowsets::workflow_set(
  preproc = list("nzv" = dat_rec),
  models = list(plogr = mod_plogr, rf = mod_rf, #btree = mod_btree, 
                nn = mod_nn, svml = mod_svml, svmp = mod_svmp, svmr = mod_svmr)
)
wfid_names <- all_wfs$wflow_id

all_wfs <- all_wfs %>%
  workflowsets::option_add(id = wfid_names[which(grepl("rf", wfid_names))], grid = rf_grid)
 
ncores <- floor(parallel::detectCores(logical = TRUE)/2)
cl <- parallel::makePSOCKcluster(ncores)
registerDoParallel(cl)
allwf_res <- all_wfs %>%
  workflowsets::workflow_map(resamples = dat_folds,
                             control = tune::control_grid(parallel_over = "everything"),
                             verbose = TRUE, 
                             seed = 1,
                             metrics = test_metrics)
parallel::stopCluster(cl)

# Extract best training results
train_perf_allmods <- allwf_res %>%
  workflowsets::rank_results(select_best = TRUE, rank_metric = "pr_auc")

# Extract best of each model
best_rf <- allwf_res %>%
  workflowsets::extract_workflow_set_result(id = "nzv_rf") %>%
  tune::select_best(metric = "pr_auc")

# Fit best of each model
fitted_rf <- allwf_res %>%
  workflowsets::extract_workflow(id = "nzv_rf") %>%
  tune::finalize_workflow(best_rf) %>%
  parsnip::fit(data = train_dat)

# Predictive performance
predperf_rf <- predict(fitted_rf, test_dat) %>%
  dplyr::bind_cols(predict(fitted_rf, test_dat, type = "prob")) %>%
  dplyr::bind_cols(test_dat %>% dplyr::select(met_bin)) %>%
  dplyr::rename(pred_class = .pred_class,
                pred_Present = .pred_Present,
                pred_Absent = `.pred_Absent`) %>%
  test_metrics(truth = met_bin, `pred_Present`, estimate = pred_class) %>%
  mutate(Metabolite = met,
         present_count = sum(test_dat$met_bin == "Present"),
         present_perc = present_count/nrow(test_dat),
         model = "rand_forest")


# saveRDS(met_fits, file = "C:\\Users\\flor829\\OneDrive - PNNL\\Projects\\NMR\\FY23_Recommender\\Results\\metmods_nodupes.rds")
# rm(met_fits, pb, i, fconc1, j)

