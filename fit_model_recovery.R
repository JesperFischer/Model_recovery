
############ setup

pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr, rstan, brms, faux)

source(here::here("scripts","Fitmodelrecovery_matlab.R"))

data = c("gaussian","gumbel","hyper", "logit")
id = 1:100

parameters = expand.grid(data = data,
                         id = id) %>% mutate(ids = 1:n())

data_list <- split(parameters, parameters$ids)

cores = parallelly::availableCores() / 4

plan(multisession, workers = cores)

possfit_model = possibly(.f = Model_recovery, otherwise = "Error")

#test it works:

test = possfit_model(data_list[[1]])

test

#full run

results <- future_map(data_list, ~possfit_model(.x), .options = furrr_options(seed = TRUE),.progress = TRUE) 

#save the results:
saveRDS(results,here::here("model_recovery.rds"))
