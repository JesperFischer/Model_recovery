


make_datasets = function(parameters){
  
  source(here::here("Model recovery","scripts", "make_data_scripts.R"))
  
  
  
  a = get_things_2con(parameters = data.frame(subjects = parameters$subjects,
                                              mu_alpha = rnorm(1,0,5),
                                              sd_alpha = rnorm(1,5,2),
                                              mu_beta = rnorm(1,-2, 1),
                                              sd_beta = rnorm(1, 1, 0.5),
                                              mu_lambda = rnorm(1,-4,0.5),
                                              sd_lambda = rnorm(1,2,0.5)
  )) %>% mutate(participant_id = 1:parameters$subjects)
  

  
  #getting Psi stimuli
  df = a %>% arrange(participant_id) %>% mutate(trials = parameters$trials)
  
  df$model = parameters$model
  
  print("Timer for PSI")
  tictoc::tic()
  ########################################################################################################### GETTING PSI!
  data = get_psi_stim(df)
  ###########################################################################################################
  tictoc::toc()
  
  #wrangling the data
  
  data = inner_join(df,data %>% dplyr::select(resp,X,participant_id,sessions,Estimatedthreshold, Estimatedslope,q5_threshold,q95_threshold,q5_slope,q95_slope), by = c("participant_id"))
  
  
  
  data = data %>% mutate(iter = parameters$id)%>% 
    mutate(real_effectsize_alpha = parameters$effect_size_alpha,
           real_effectsize_beta = parameters$effect_size_beta,
           subjects = parameters$subjects,
           trials = parameters$trials,
           model = parameters$model)
  
  directory = paste0("trials = ", df$trials[1],
                     " subjects = ", length(unique(df$participant_id))[1],
                     " iter = ", parameters$id[1],
                     " random id = ", round(rnorm(1,0,1000),2),".csv")
  
  model = as.character(parameters$model[[1]])
  
  dir1 = paste0("trials = ", df$trials[1],
                " subjects = ", length(unique(df$participant_id))[1])
  
  if(!dir.exists(here::here("Model recovery","Data",model))){
    dir.create(here::here("Model recovery","Data",model))
  }
  
  if(!dir.exists(here::here("Model recovery","Data",model,dir1))){
    dir.create(here::here("Model recovery","Data",model,dir1))
  }
  data = data %>% unnest()
  write.csv(data,here::here("Model recovery","Data",model,dir1,directory))
  
}




erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1


get_things_2con = function(parameters){
  
  alpha = rnorm(parameters$subjects, parameters$mu_alpha, parameters$sd_alpha)
  
  beta = exp(rnorm(parameters$subjects, parameters$mu_beta, parameters$sd_beta))
  
  lambda = brms::inv_logit_scaled(rnorm(parameters$subjects, parameters$mu_lambda, parameters$sd_lambda) ) / 2
  
  parameters2 = data.frame(alpha = alpha, beta = beta, lambda = lambda)

  return(parameters2)
}

get_psi_stim = function(parameters){
  
  if(parameters$model[1] == "normal"){
    python_script <- here::here("Model recovery","python","PSI_normal.py")
  }else if(parameters$model[1] == "logistic"){
    python_script <- here::here("Model recovery","python","PSI_logistic.py")
  }else if(parameters$model[1] == "gumbel"){
  python_script <- here::here("Model recovery","python","PSI_gumbel.py")
  }else if(parameters$model[1] == "hyperbolic"){
    python_script <- here::here("Model recovery","python","PSI_hyperbolic.py")
  }

  alpha = parameters$alpha
  
  beta = parameters$beta
  
  lapse = parameters$lambda
  
  trials = as.integer(parameters$trials)
  
  ids = as.integer(parameters$participant_id)
  
  subjects = max(as.integer(unique(parameters$participant_id)))
  
  sessions = as.integer(1)
  
  library(reticulate)
  
  # Use reticulate to run the Python script with arguments
  
  source_python(python_script, convert = FALSE)
  
  d = get_stim(lapse, alpha, beta, ids, trials, subjects, sessions)
  
  #Produces a warning about will be removed in future versions (but seems to be a thing with reticulate and not my code)
  dd = reticulate::py_to_r(d)
  dd = reticulate::py_to_r(d)
  
  d = data.frame(lapse = dd$lapse, alpha = dd$alpha, beta = dd$beta, participant_id = unlist(dd$participant_id),
                 trials = unlist(dd$trials), subs = unlist(dd$subs), X = unlist(dd$X), resp = unlist(dd$resp), sessions = unlist(dd$sessions),
                 Estimatedthreshold = unlist(dd$Estimatedthreshold), Estimatedslope = unlist(dd$Estimatedslope), q5_threshold =  unlist(dd$q5_threshold),
                 q95_threshold =  unlist(dd$q95_threshold), q5_slope =  unlist(dd$q5_slope), q95_slope =  unlist(dd$q95_slope))
  
  
  #plot for it to make sense:
  # d  %>%  group_by(participant_id,sessions) %>% mutate(trial = 1:n()) %>%
  #   ggplot(aes(x = trial, y = Estimatedslope, col = sessions))+geom_pointrange(aes(ymin = q5_slope, ymax = q95_slope))+
  #   facet_wrap(~participant_id)+geom_hline(aes(yintercept = beta, col = sessions))
  # 
  # 
  # d  %>%  group_by(participant_id,sessions) %>% mutate(trial = 1:n()) %>%
  #   ggplot(aes(x = trial, y = Estimatedthreshold, col = sessions))+geom_pointrange(aes(ymin = q5_threshold, ymax = q95_threshold))+
  #   facet_wrap(~participant_id)+geom_hline(aes(yintercept = alpha, col = sessions))
  # 
  
  return(d)
  
}





#function to get confidence intervals used in plot_intervals.
ci = function(x){
  list = list(which(cumsum(x)/sum(x) > 0.025)[1],
              last(which(cumsum(x)/sum(x) < 0.975)))
  return(list)
  
}

#function to get the line intervals from the exteropost and interopost dataframe
get_line_intervals = function(data, parameter){
  
  
  if(parameter == "alpha"){
    confidence = ci(rowMeans(data))
  }else if(parameter == "beta"){
    confidence = ci(colMeans(data))
  }
  rg = seq(-50.5,50.5, by = 1)
  rg = seq(0.1, 25, by = 0.1)
  
  
  upper = rg[confidence[[1]]]
  lower = rg[confidence[[2]]]
  
  
  data = data.frame(upper = upper,lower = lower)
  
  return(data)
}
