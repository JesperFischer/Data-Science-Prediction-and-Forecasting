#model recovery


model_fitter_rwdata = function(parameters_rw){
  
  
  df_rw = our_hier_rw_agent(parameters_rw)[[1]]
  
  #fitting hierachically
  data_rw = df_rw %>% filter(id %in% unique(id)[1:nsubs])
  
  data_rw = data_rw %>% mutate(across(percept, ~case_when(
    . < 0.001 ~ 0.001,
    . > 0.999 ~ 0.999,
    TRUE ~ .
  )))
  
  
  data1_rw = list(nsubs = length(unique(data_rw$id)),
               ntrials = nrow(data_rw %>% filter(id %in% unique(id)[1])),
               percept = as.matrix(data_rw %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL)),
               expectPain = as.matrix(data_rw %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL)),
               percept_bin = as.matrix(data_rw %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL)),
               stim = as.matrix(data_rw %>% pivot_wider(id_cols = trial, names_from = id, values_from = stim)%>% mutate(trial= NULL)),
               cues = as.matrix(data_rw %>% pivot_wider(id_cols = trial, names_from = id, values_from = cue)%>% mutate(trial= NULL)))
  
  id = rnorm(1,0,1)
  
  
  
  rw_model = here::here("Stan","myRW_real.stan")
  wb_model = here::here("Stan","hier_weighted_bayes.stan")
  logs = here::here("Stan","rw_vs_rw_hier_complete.stan")
  
  data1_rw_wb = data1_rw
  data1_rw_logs = data1_rw
  
  names(data1_rw_wb) = c("nsubs","ntrials","percept","pred","percept_bin","stim","cue")
  names(data1_rw_logs) = c("nsubs","ntrials","percept","pred","percept_bin","stim","cues")
  data1_rw_wb$stim = ifelse(data1_rw_wb$stim == 1, 0.999, 0.001)
  
  data1_rw_logs$u = as.matrix(data_rw %>% pivot_wider(id_cols = trial, names_from = id, values_from = u)%>% mutate(trial= NULL))
  
  
  loo_rw_rw = get_loo(rw_model, data1_rw)
  loo_rw_wb = get_loo(wb_model, data1_rw_wb)
  loo_rw_logs = get_loo(logs, data1_rw_logs)
  
  
  comparison = data.frame(loo_compare(list(loo_rw_rw[[3]],loo_rw_wb[[3]],loo_rw_logs[[3]])))
  comparison$id = id

  return(list(data.frame(names = c("rw",
                                   "wb",
                                   "logs"),
                         means = c(loo_rw_rw[[1]],loo_rw_wb[[1]],loo_rw_logs[[1]]),
                         sds = c(loo_rw_rw[[2]],loo_rw_wb[[2]],loo_rw_logs[[2]]),
                         id = id),
              comparison))
}


model_fitter_wbdata = function(parameters_wb){
  
  
  df_wb = hier_Bayes_f(parameters_wb)[[1]]
  
  #fitting hierachically
  data_wb = df_wb %>% filter(id %in% unique(id)[1:nsubs])
  
  data_wb = data_wb %>% mutate(across(percept, ~case_when(
    . < 0.001 ~ 0.001,
    . > 0.999 ~ 0.999,
    TRUE ~ .
  )))
  
  data_wb$stim = ifelse(data_wb$stim > 0.5, 1, 0)
  
  data1_wb = list(nsubs = length(unique(data_wb$id)),
                  ntrials = nrow(data_wb %>% filter(id %in% unique(id)[1])),
                  percept = as.matrix(data_wb %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL)),
                  pred = as.matrix(data_wb %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL)),
                  percept_bin = as.matrix(data_wb %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL)),
                  stim = as.matrix(data_wb %>% pivot_wider(id_cols = trial, names_from = id, values_from = stim)%>% mutate(trial= NULL)),
                  cue = as.matrix(data_wb %>% pivot_wider(id_cols = trial, names_from = id, values_from = cue)%>% mutate(trial= NULL)))
  
  
  data1_wb_rw = data1_wb
  data1_rw_logs = data1_wb
  
  names(data1_wb_rw) = c("nsubs","ntrials","percept","expectPain","percept_bin","stim","cues")
  names(data1_rw_logs) = c("nsubs","ntrials","percept","pred","percept_bin","stim","cues")
  
  data1_rw_logs$u = as.matrix(data_wb %>% pivot_wider(id_cols = trial, names_from = id, values_from = u)%>% mutate(trial= NULL))
  
  
  data1_wb$stim = ifelse(data1_wb$stim == 1, 0.999, 0.001)
  
  rw_model = here::here("Stan","myRW_real.stan")
  wb_model = here::here("Stan","hier_weighted_bayes.stan")
  logs = here::here("Stan","rw_vs_rw_hier_complete.stan")
  
  id = rnorm(1,0,1)
  
  loo_wb_rw = get_loo(rw_model, data1_wb_rw)
  loo_wb_wb = get_loo(wb_model, data1_wb)
  loo_rw_logs = get_loo(logs, data1_rw_logs)
  
  comparison = data.frame(loo_compare(list(loo_wb_rw[[3]],loo_wb_wb[[3]],loo_rw_logs[[3]])))
  comparison$id = id
  
  return(list(data.frame(names = c("rw",
                                   "wb",
                                   "kalman"),
                         means = c(loo_wb_rw[[1]],loo_wb_wb[[1]],loo_rw_logs[[1]]),
                         sds = c(loo_wb_rw[[2]],loo_wb_wb[[2]],loo_rw_logs[[2]]),
                         id = id),
              comparison))

  
}



model_fitter_kalmandata = function(parameters_kalman){
  
  
  df_kalman = hier_rw_logistic(parameters_kalman)[[1]]
  
  #fitting hierachically
  data_kalman = df_kalman %>% filter(id %in% unique(id)[1:nsubs])
  
  data_kalman = data_kalman %>% mutate(across(percept, ~case_when(
    . < 0.001 ~ 0.001,
    . > 0.999 ~ 0.999,
    TRUE ~ .
  )))
  
  
  
  
  data1_kalman = list(nsubs = length(unique(data_kalman$id)),
                      ntrials = nrow(data_kalman %>% filter(id %in% unique(id)[1])),
                      percept = as.matrix(data_kalman %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL)),
                      pred = as.matrix(data_kalman %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL)),
                      percept_bin = as.matrix(data_kalman %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL)),
                      stim = as.matrix(data_kalman %>% pivot_wider(id_cols = trial, names_from = id, values_from = stim)%>% mutate(trial= NULL)),
                      cues = as.matrix(data_kalman %>% pivot_wider(id_cols = trial, names_from = id, values_from = cue)%>% mutate(trial= NULL)),
                      u = as.matrix(data_kalman %>% pivot_wider(id_cols = trial, names_from = id, values_from = u)%>% mutate(trial= NULL))
  )
  
  
  rw_model = here::here("Stan","myRW_real.stan")
  wb_model = here::here("Stan","hier_weighted_bayes2.stan")
  logs = here::here("Stan","rw_vs_rw_hier_complete.stan")
  
  
  
  data1_kalman_rw = data1_kalman
  data1_kalman_wb = data1_kalman
  
  data1_kalman_wb$u = NULL
  data1_kalman_rw$u = NULL
  names(data1_kalman_rw) = c("nsubs","ntrials","percept","expectPain","percept_bin","stim","cues")
  names(data1_kalman_wb) = c("nsubs","ntrials","percept","pred","percept_bin","stim","cue")
  data1_kalman_wb$stim = ifelse(data1_kalman_wb$stim == 1, 0.999, 0.001)
  
  
  
  loo_kalman_rw = get_loo(rw_model, data1_kalman_rw)
  loo_kalman_wb = get_loo(wb_model, data1_kalman_wb)
  loo_kalman_kalman = get_loo(logs, data1_kalman)

  id = rnorm(1,0,1)
  
  comparison = data.frame(loo_compare(list(loo_kalman_rw[[3]],loo_kalman_wb[[3]],loo_kalman_kalman[[3]])))
  comparison$id = id
  
  return(list(data.frame(names = c("rw",
                                   "wb",
                                   "kalman"),
                         means = c(loo_kalman_rw[[1]],loo_kalman_wb[[1]],loo_kalman_kalman[[1]]),
                         sds = c(loo_kalman_rw[[2]],loo_kalman_wb[[2]],loo_kalman_kalman[[2]]),
                         id = id),
              comparison))
  
  
  
    
  
}



get_loo = function(file, data1){
  mod = cmdstan_model(file)
  
  fit <- mod$sample(
    data = data1, 
    seed = 123, 
    chains = 4, 
    parallel_chains = 4,
    refresh = 100
  )
  
  loglik = as_draws_df(fit$draws(variables = "log_lik")) %>% drop_na()
  q = loo::loo(as.matrix(loglik))
  
  
  loo1 = q$estimates[3,1]
  loo1sd = q$estimates[3,2]
  
  return(list(loo1,loo1sd,q))
}

