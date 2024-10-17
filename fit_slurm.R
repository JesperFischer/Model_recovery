

fit_model_recovery = function(){
  
  pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr,bridgesampling, rstan, brms, faux,LRO.utilities,reticulate)
  message("Number of CPU cores in R: ", parallelly::availableCores())
  
  source(here::here("Model recovery", "scripts","Fit model recovery.R"))


  subjects = 20
  trials = 50
  model = c("gaussian","gumbel","hyper", "logit")
  replicate = 1:100
  
  print(subjects)
  
  print(trials)
  parameters = expand.grid(subjects = subjects,
                           trials = trials,
                           model = model,
                           replicate = replicate) %>% mutate(ids = 1:n()) %>% 
    group_by(model) %>% 
    mutate(id = 1:n())
  
  data_list <- split(parameters, parameters$ids)
  
  cores = parallelly::availableCores() / 2
  # cores = 40 / 2  
#  q = Model_recovery_rstan(data_list[[1]])
 # qq = Model_recovery(data_list[[1]])
  
  plan(multisession, workers = cores)
  
  possfit_model = possibly(.f = Model_recovery, otherwise = "Error")
  
  # results <- future_map(data_list, ~possfit_model(.x), .options = furrr_options(seed = TRUE),.progress = TRUE)

  results <- future_map(data_list, ~possfit_model(.x), .options = furrr_options(seed = TRUE)) 
  
  saveRDS(results,here::here("Model recovery","results","moment_match",paste0("subjects = ", subjects, " trials = ",trials," results.rds")))

  
  #map_dfr(results,bind_rows) %>% filter(elpd_diff == 0) %>% group_by(dataset, models) %>% summarize(n()) 
  
  
}


fit_model_recovery()