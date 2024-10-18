

fit_model_recovery = function(){
  
  pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr,bridgesampling, rstan, brms, faux,LRO.utilities,reticulate)
  message("Number of CPU cores in R: ", parallelly::availableCores())
  
  source(here::here("Model recovery","scripts", "make_data_scripts.R"))
  
  subjects = 50
  trials = 150
  model = c("normal","logistic","hyperbolic","gumbel")
  replicate = 1:100
  
  
  parameters = expand.grid(subjects = subjects,
                           trials = trials,
                           model = model,
                           replicate = replicate) %>% 
    mutate(id = 1:nrow(.))
  
  data_list <- split(parameters, parameters$id)
  
  cores = 15
  print("running")
  plan(multisession, workers =  parallelly::availableCores())
  
  possfit_model = possibly(.f = make_datasets, otherwise = "Error")
  
  results <- future_map(data_list, ~possfit_model(.x), .options = furrr_options(seed = TRUE))
}

fit_model_recovery()