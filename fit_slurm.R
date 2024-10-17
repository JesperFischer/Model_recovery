

fit_model_recovery = function(){
  
  pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr,bridgesampling, rstan, brms, faux,LRO.utilities,reticulate)
  message("Number of CPU cores in R: ", parallelly::availableCores())
  
  source(here::here("scripts","Fitmodelrecovery_matlab.R"))


  data = c("gaussian","gumbel","hyper", "logit")
  id = 1:100
  
  parameters = expand.grid(data = data,
                           id = id) %>% mutate(ids = 1:n())
  
  data_list <- split(parameters, parameters$ids)
  
  cores = parallelly::availableCores() / 2
  
  plan(multisession, workers = cores)
  
  possfit_model = possibly(.f = Model_recovery, otherwise = "Error")
  
  test = possfit_model(data_list[[1]])

  results <- future_map(data_list, ~possfit_model(.x), .options = furrr_options(seed = TRUE),.progress = TRUE) 
  
  saveRDS(results,here::here("model_recovery.rds"))

  
  #map_dfr(results,bind_rows) %>% filter(elpd_diff == 0) %>% group_by(dataset, models) %>% summarize(n()) 
  
  
}


fit_model_recovery()