

Model_recovery = function(parameters){
  
  source(here::here("Model recovery","scripts", "Fit model recovery.R"))

  model = paste0("model = ",
                   as.character(parameters$model[1]))
  
  trials_n_subs = paste0("trials = ",
                         as.character(parameters$trials[1]),
                         " subjects = ",
                         as.character(parameters$subjects[1]))
  
  
  # Get a list of all files in the directory
  files <- sort(list.files(path = here::here("Model recovery","Data",as.character(parameters$model[1]),trials_n_subs), full.names = TRUE))
  
  data = read.csv(files[parameters$id])
  
  #rstanmod = rstan::stan_model(here::here("Model recovery","Stanmodels","mixed parameterization","Normal.stan"))
  
  mod_norm_mixed = cmdstanr::cmdstan_model(here::here("Model recovery","Stanmodels","Loo_trial","mixed","Mixed_Normal.stan"),   force_recompile = TRUE)
  mod_Hyperbolic_mixed = cmdstanr::cmdstan_model(here::here("Model recovery","Stanmodels","Loo_trial","mixed","Mixed_Hyperbolic.stan"))
  mod_Gumbel_mixed = cmdstanr::cmdstan_model(here::here("Model recovery","Stanmodels","Loo_trial","mixed","Mixed_Gumbel.stan"))
  mod_Logistic_mixed = cmdstanr::cmdstan_model(here::here("Model recovery","Stanmodels","Loo_trial","mixed","Mixed_Logistic.stan"))
  
  #rstanmod = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","Normal.stan"))
  #rstanmod = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","mixed","Mixed_Normal.stan"))
  
  
  mod_norm = cmdstanr::cmdstan_model(here::here("Model recovery","Stanmodels","Loo_trial","Normal.stan"))
  mod_Hyperbolic = cmdstanr::cmdstan_model(here::here("Model recovery","Stanmodels","Loo_trial","Hyperbolic.stan"))
  mod_Gumbel = cmdstanr::cmdstan_model(here::here("Model recovery","Stanmodels","Loo_trial","Gumbel.stan"))
  mod_Logistic = cmdstanr::cmdstan_model(here::here("Model recovery","Stanmodels","Loo_trial","Logistic.stan"))
  
  #mod_norm = cmdstanr::cmdstan_model(here::here("Model recovery","Stanmodels","Loo_subject","Normal_loo_subject.stan"),stanc_options = list("O1"))
  
  
  normal_fit1 = fit_model(data,mod_norm, model = "normal")
  
  if(normal_fit1[[2]]$mean_div != 0){
    normal_fit2 = fit_model(data,mod_norm_mixed, model = "normal")
    if(normal_fit1[[2]]$mean_div < normal_fit2[[2]]$mean_div){
      normal_fit = normal_fit1
    }else{
      normal_fit = normal_fit2
    }
  }else{
    normal_fit = normal_fit1
  }

  hyper_fit1 = fit_model(data,mod_Hyperbolic, model = "hyperbolic")

  if(hyper_fit1[[2]]$mean_div != 0){
    hyper_fit2 = fit_model(data,mod_Hyperbolic_mixed, model = "hyperbolic")
    
    if(hyper_fit1[[2]]$mean_div < hyper_fit2[[2]]$mean_div){
      hyper_fit = hyper_fit1
    }else{
      hyper_fit = hyper_fit2
    }
  }else{
    hyper_fit = hyper_fit1
  }
  
    
  gumbel_fit1 = fit_model(data,mod_Gumbel_mixed, model = "gumbel")
  
  
  if(gumbel_fit1[[2]]$mean_div != 0){
    gumbel_fit2 = fit_model(data,mod_Gumbel, model = "gumbel")
    if(gumbel_fit1[[2]]$mean_div < gumbel_fit2[[2]]$mean_div){
      gumbel_fit = gumbel_fit1
    }else{
      gumbel_fit = gumbel_fit2
    }
  }else{
    gumbel_fit = gumbel_fit1
  }
  
  logistic_fit1 = fit_model(data,mod_Logistic_mixed, model = "logistic")
  
  
  if(logistic_fit1[[2]]$mean_div != 0){
    logistic_fit2 = fit_model(data,mod_Logistic, model = "logistic")

    if(logistic_fit1[[2]]$mean_div < logistic_fit2[[2]]$mean_div){
      logistic_fit = logistic_fit1
    }else{
      logistic_fit = logistic_fit2
    }
  }else{
    logistic_fit = logistic_fit1
  }
  
  
  
  
  
  loo = loo::loo_compare(list(normal = normal_fit[[1]],
                   hyperbolic = hyper_fit[[1]],
                   gumbel = gumbel_fit[[1]],
                   logistic = logistic_fit[[1]]))
  
  
  fit_diagnostics = rbind(normal_fit[[2]],
                      hyper_fit[[2]],
                      gumbel_fit[[2]],
                      logistic_fit[[2]])
  
  
  result = inner_join(data.frame(loo) %>% rownames_to_column(var = "models") %>% 
                        dplyr::select(models, elpd_diff,se_diff) %>% mutate(elpd_ratio = elpd_diff/se_diff), fit_diagnostics) %>% 
    mutate(replicate = parameters$replicate,
           subjects = parameters$subjects,
           trials = parameters$trials,
           dataset = parameters$model)
  
  return(result)
  
}



fit_model = function(data,mod, model){
  
  
  # if(max(data$sessions) != 1){
  # 
  #   df = data%>% filter(sessions == 1) %>% group_by(participant_id,X) %>% summarize(yn = sum(resp), n = n())
  # }else{
  #   df = data %>% group_by(participant_id,X) %>% summarize(yn = sum(resp), n = n())
  # }

  datastan = list(Y = data$resp,
                  N = nrow(data),
                  S = length(unique(data$participant_id)),
                  S_id = data$participant_id,
                  X = matrix(c(rep(1,nrow(data)), data$X), ncol = 2, nrow = nrow(data)))
  
    
  # seems to be the case that the non-centered is just better.
  
  # fit = sampling(mod, data = datastan,
  #                iter = 3000,
  #                chains = 4,
  #                init = 0,
  #                cores = 4,
  #                refresh = 500)
    #fitting
    fit <- mod$sample(
      data = datastan,
      iter_sampling = 1500,
      iter_warmup = 1500,
      chains = 4,
      init = 0,
      parallel_chains = 4,
      refresh = 500,
      adapt_delta = 0.90,
      max_treedepth = 10
    )
    
      
      diags_norm = data.frame(fit$diagnostic_summary()) %>% 
        mutate(models = model) %>% group_by(models) %>% summarize(mean_div = mean(num_divergent),
                                                                   mean_treedepth = mean(num_max_treedepth))
      
      rhat_norm = data.frame(fit$summary(c("gm[1]","gm[2]","gm[3]",
                                                  "tau_u[1]","tau_u[2]","tau_u[3]"))) %>% 
        mutate(models = model) %>% group_by(models) %>% summarize(meanrhat = mean(rhat))
      

    
      options(mc.cores = 4)
      
      fit_loo = loo(fit$draws("log_lik"),moment_match=T)
      
      # fit_loomm = fit$loo(moment_match = T)
      # 
      # fit_loo = fit$loo()
      
      normal_loo_diag = sum(fit_loo$diagnostics$pareto_k > 0.7)
      
      diags = inner_join(diags_norm,rhat_norm) %>% mutate(pareto_k_over_0.7 = normal_loo_diag)
            
      return(list(fit_loo, diags,fit))
}


Model_recovery_rstan = function(parameters){
  
  source(here::here("Model recovery","scripts", "Fit model recovery.R"))
  
  model = paste0("model = ",
                 as.character(parameters$model[1]))
  
  trials_n_subs = paste0("trials = ",
                         as.character(parameters$trials[1]),
                         " subjects = ",
                         as.character(parameters$subjects[1]))
  
  
  # Get a list of all files in the directory
  files <- sort(list.files(path = here::here("Model recovery","Data",as.character(parameters$model[1]),trials_n_subs), full.names = TRUE))
  
  data = read.csv(files[parameters$id])
  

  mod_norm = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","Normal.stan"))
  mod_norm_mixed = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","mixed","Mixed_Normal.stan"))
  
  
  
  normal_fit1 = fit_model_rstan(data,mod_norm, model = "normal")
  
  if(normal_fit1[[2]]$sum_div != 0){
    normal_fit2 = fit_model_rstan(data,mod_norm_mixed, model = "normal")
    if(normal_fit1[[2]]$sum_div < normal_fit2[[2]]$sum_div){
      normal_fit = normal_fit1
    }else{
      normal_fit = normal_fit2
    }
  }else{
    normal_fit = normal_fit1
  }
  
  
  mod_Hyperbolic = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","Hyperbolic.stan"))
  mod_Hyperbolic_mixed = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","mixed","Mixed_Hyperbolic.stan"))
  
  
  hyper_fit1 = fit_model_rstan(data,mod_Hyperbolic, model = "hyperbolic")
  
  if(hyper_fit1[[2]]$sum_div != 0){
    hyper_fit2 = fit_model_rstan(data,mod_Hyperbolic_mixed, model = "hyperbolic")
    
    if(hyper_fit1[[2]]$sum_div < hyper_fit2[[2]]$sum_div){
      hyper_fit = hyper_fit1
    }else{
      hyper_fit = hyper_fit2
    }
  }else{
    hyper_fit = hyper_fit1
  }
  
  
  mod_Gumbel = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","Gumbel.stan"))
  mod_Gumbel_mixed = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","mixed","Mixed_Gumbel.stan"))
  
  
  gumbel_fit1 = fit_model_rstan(data,mod_Gumbel, model = "gumbel")
  
  if(gumbel_fit1[[2]]$sum_div != 0){
    gumbel_fit2 = fit_model_rstan(data,mod_Gumbel_mixed, model = "gumbel")
    if(gumbel_fit1[[2]]$sum_div < gumbel_fit2[[2]]$sum_div){
      gumbel_fit = gumbel_fit1
    }else{
      gumbel_fit = gumbel_fit2
    }
  }else{
    gumbel_fit = gumbel_fit1
  }
  
  
  mod_Logistic = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","Logistic.stan"))
  mod_Logistic_mixed = rstan::stan_model(here::here("Model recovery","Stanmodels","Loo_trial","mixed","Mixed_Logistic.stan"))
  
  
  logistic_fit1 = fit_model_rstan(data,mod_Logistic, model = "logistic")
  
  
  if(logistic_fit1[[2]]$sum_div != 0){
    logistic_fit2 = fit_model_rstan(data,mod_Logistic, model = "logistic")
    
    if(logistic_fit1[[2]]$sum_div < logistic_fit2[[2]]$sum_div){
      logistic_fit = logistic_fit1
    }else{
      logistic_fit = logistic_fit2
    }
  }else{
    logistic_fit = logistic_fit1
  }
  
  
  loo = loo::loo_compare(list(normal = normal_fit[[1]],
                              hyperbolic = hyper_fit[[1]],
                              gumbel = gumbel_fit[[1]],
                              logistic = logistic_fit[[1]]))
  
  
  fit_diagnostics = rbind(normal_fit[[2]],
                          hyper_fit[[2]],
                          gumbel_fit[[2]],
                          logistic_fit[[2]])
  
  
  result = inner_join(data.frame(loo) %>% rownames_to_column(var = "models") %>% 
                        dplyr::select(models, elpd_diff,se_diff) %>% mutate(elpd_ratio = elpd_diff/se_diff), fit_diagnostics) %>% 
    mutate(replicate = parameters$replicate,
           subjects = parameters$subjects,
           trials = parameters$trials,
           dataset = parameters$model)
  
  return(result)
  
}


fit_model_rstan = function(data,mod, model){
  

  
  datastan = list(Y = data$resp,
                  N = nrow(data),
                  S = length(unique(data$participant_id)),
                  S_id = data$participant_id,
                  X = matrix(c(rep(1,nrow(data)), data$X), ncol = 2, nrow = nrow(data)))
  
  
  # seems to be the case that the non-centered is just better.
  
  fit = sampling(mod, data = datastan,
                 iter = 2000,
                 chains = 4,
                 init = 0,
                 cores = 4,
                 refresh = 500)
  
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  num_divergent <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
  total_divergent <- sum(num_divergent)
  
  treedepth <- sampler_params[[1]][,3]
  total_treedepth <- sum(treedepth == 10)
  

  diags_norm = data.frame(models = model) %>% mutate(sum_div = total_divergent,
                                                     sum_tree = total_treedepth)
  
  
  rhat_norm = data.frame(models = model,rhat = sum(summary(fit)$summary[1:6,10] > 1.03))
  
  fit_loo = loo(fit,moment_match=T, cores = 4)
  
  normal_loo_diag = pareto_over_0.7 = sum(fit_loo$diagnostics$pareto_k >0.7)
  
  diags = inner_join(diags_norm,rhat_norm) %>% mutate(pareto_k_over_0.7 = normal_loo_diag)
  
  return(list(fit_loo, diags,fit))
}
