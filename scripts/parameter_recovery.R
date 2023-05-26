#parameter_recovery


parameter_recovery_rw = function(parameters){
  
  df = our_hier_rw_agent(parameters)[[1]]
  hier = our_hier_rw_agent(parameters)[[2]]
  
  #fitting hierachically
  data = df %>% filter(id %in% unique(id)[1:nsubs])
  
  data = data %>% mutate(across(percept, ~case_when(
    . < 0.001 ~ 0.001,
    . > 0.999 ~ 0.999,
    TRUE ~ .
  )))
  
  
  
  data1 = list(nsubs = length(unique(data$id)),
               ntrials = nrow(data %>% filter(id %in% unique(id)[1])),
               percept = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL)),
               expectPain = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL)),
               percept_bin = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL)),
               stim = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = stim)%>% mutate(trial= NULL)),
               cues = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = cue)%>% mutate(trial= NULL)))
  
  
  mod = cmdstan_model(here::here("Stan","no generated quantities","myRW_real_nogen.stan"))
  

  
  fitter = function(){
    fit <- mod$sample(
      data = data1,
      chains = 4, 
      parallel_chains = 4,
      refresh = 100,
      adapt_delta = 0.90,
      max_treedepth = 12
    )
    return(fit)
  }
  
  
  fit <- tryCatch({
    fit = withTimeout(fitter(), timeout = 1000, cpu = 1000)
  }, error = function(e) {
    print("An error occurred!")
    return(NA)
  })
  
  
  
  if(is.environment(fit)){
    
    index = rnorm(1,0,1)
    
    hier_parameters = c("kappa_w1","mu_w1","kappa_alpha","mu_alpha","sd_beta","sd_precision_percept")
    
    sub_parameters = c("w1","alpha","beta","precision_percept")
    return(list(hier = data.frame(fit$summary(hier_parameters),
                                  reals = hier %>% dplyr::select(all_of(hier_parameters)) %>% pivot_longer(everything()) %>% rename(reals = value),
                                  div = sum(fit$diagnostic_summary()$num_divergent),index = index),
                sub = data.frame(fit$summary(sub_parameters) %>% arrange(variable),
                                 reals = df %>% filter(trial == 1) %>% mutate(ids = 1:nrow(.)) %>% dplyr::select(all_of(sub_parameters),ids) %>% pivot_longer(cols = -ids) %>% rename(reals =value) %>% mutate(names = paste0(name,"[",ids,"]")) %>% arrange(names) %>% select(names,reals),
                                 div = sum(fit$diagnostic_summary()$num_divergent),index = index)))}else{
                                   return(list(NA,NA))
                                 }
}


parameter_recovery_wb = function(parameters){
  
  df = hier_Bayes_f(parameters)[[1]]
  hier = hier_Bayes_f(parameters)[[2]]
  
  #fitting hierachically
  data = df %>% filter(id %in% unique(id)[1:nsubs])
  
  data = data %>% mutate(across(percept, ~case_when(
    . < 0.001 ~ 0.001,
    . > 0.999 ~ 0.999,
    TRUE ~ .
  )))
  
  
  
  data1 = list(nsubs = length(unique(data$id)),
               ntrials = nrow(data %>% filter(id %in% unique(id)[1])),
               percept = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL)),
               pred = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL)),
               percept_bin = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL)),
               stim = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = stim)%>% mutate(trial= NULL)),
               cue = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = cue)%>% mutate(trial= NULL)))
  
  
  
  

  mod = cmdstan_model(here::here("Stan","no generated quantities","hier_weighted_bayes_nogen.stan"))
  
  
  
  fitter = function(){
    fit <- mod$sample(
      data = data1,
      chains = 4, 
      parallel_chains = 4,
      refresh = 100,
      adapt_delta = 0.90,
      max_treedepth = 12
    )
    return(fit)
  }
  
  
  fit <- tryCatch({
    fit = withTimeout(fitter(), timeout = 1500)
  }, error = function(e) {
    print("An error occurred!")
    message(e)
    return(NA)
  }, warning = FALSE)
  
  
  
  if(is.environment(fit)){
    
    index = rnorm(1,0,1)
    
    hier_parameters = c("kappa_alpha","mu_alpha","mu_w1","kappa_w1","mu_w2","kappa_w2","sd_beta","sd_percept_precision")
    
    sub_parameters = c("alpha","w1","w2","beta","percept_precision")
    return(list(hier = data.frame(fit$summary(hier_parameters),
                                  reals = hier %>% dplyr::select(all_of(hier_parameters)) %>% pivot_longer(everything()) %>% rename(reals = value),
                                  div = sum(fit$diagnostic_summary()$num_divergent),index = index),
                sub = data.frame(fit$summary(sub_parameters) %>% arrange(variable),
                                 reals = df %>% filter(trial == 1) %>% mutate(ids = 1:nrow(.)) %>% dplyr::select(all_of(sub_parameters),ids) %>% pivot_longer(cols = -ids) %>% rename(reals =value) %>% mutate(names = paste0(name,"[",ids,"]")) %>% arrange(names) %>% select(names,reals),
                                 div = sum(fit$diagnostic_summary()$num_divergent),index = index)))}else{
                                   return(list(NA,NA))
                                 }
}



parameter_recovery_kalman = function(parameters){
  
  df = our_hier_kalman_agent(parameters)[[1]]
  hier = our_hier_kalman_agent(parameters)[[2]]
  
  #fitting hierachically
  data = df %>% filter(id %in% unique(id)[1:nsubs])
  
  data = data %>% mutate(across(percept, ~case_when(
    . < 0.001 ~ 0.001,
    . > 0.999 ~ 0.999,
    TRUE ~ .
  )))
  
  
  
  data1 = list(nsubs = length(unique(data$id)),
               ntrials = nrow(data %>% filter(id %in% unique(id)[1])),
               percept = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL)),
               expectPain = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL)),
               percept_bin = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL)),
               stim = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = stim)%>% mutate(trial= NULL)),
               cues = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = cue)%>% mutate(trial= NULL)),
               u = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = u)%>% mutate(trial= NULL))
  )
  
  mod = cmdstan_model(here::here("Stan","no generated quantities","myKalmanfilter_v2_nogen.stan"))
  
  
  fitter = function(){
    fit <- mod$sample(
      data = data1,
      chains = 4, 
      parallel_chains = 4,
      refresh = 100,
      adapt_delta = 0.90,
      max_treedepth = 12
    )
    
    return(fit)
  }
  
  
  fit <- fitter()
  
  
  
  if(is.environment(fit)){
    
    index = rnorm(1,0,1)
    
    hier_parameters = c("mu_sigmaEpsilon","mu_sigmaEta","mu_sigmaPsi","sd_sigmaEpsilon","sd_sigmaEta","sd_sigmaPsi")
    
    sub_parameters = c("sigmaEpsilon","sigmaEta","sigmaPsi")
    return(list(hier = data.frame(fit$summary(hier_parameters),
                                  reals = hier %>% dplyr::select(all_of(hier_parameters)) %>% pivot_longer(everything()) %>% rename(reals = value),
                                  div = sum(fit$diagnostic_summary()$num_divergent),
                                  index = index),
                sub = data.frame(fit$summary(sub_parameters) %>% arrange(variable),
                                 reals = df %>% filter(trial == 1) %>% mutate(ids = 1:nrow(.)) %>% dplyr::select(all_of(sub_parameters),ids) %>% pivot_longer(cols = -ids) %>% rename(reals =value) %>% mutate(names = paste0(name,"[",ids,"]")) %>% arrange(names) %>% select(names,reals),
                                 div = sum(fit$diagnostic_summary()$num_divergent),index = index)))}else{
                                   return(list(NA,NA))
                                 }
}




parameter_recovery_logistic_rw = function(parameters){
  
  df = hier_rw_logistic(parameters)[[1]]
  hier = hier_rw_logistic(parameters)[[2]]
  
  #fitting hierachically
  data = df %>% filter(id %in% unique(id)[1:nsubs])
  
  data = data %>% mutate(across(percept, ~case_when(
    . < 0.001 ~ 0.001,
    . > 0.999 ~ 0.999,
    TRUE ~ .
  )))
  
  
  
  data1 = list(nsubs = length(unique(data$id)),
               ntrials = nrow(data %>% filter(id %in% unique(id)[1])),
               pred = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = pred)%>% mutate(trial= NULL)),
               percept_bin = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept_bin)%>% mutate(trial= NULL)),
               percept = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = percept)%>% mutate(trial= NULL)),
               stim = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = stim)%>% mutate(trial= NULL)),
               u = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = u)%>% mutate(trial= NULL)),
               cues = as.matrix(data %>% pivot_wider(id_cols = trial, names_from = id, values_from = cue)%>% mutate(trial= NULL))
  )
  
  mod = cmdstan_model(here::here("Stan","rw_vs_rw_hier_complete.stan"))
  
  
  
  fitter = function(){
    fit <- mod$sample(
      data = data1,
      chains = 4, 
      parallel_chains = 4,
      refresh = 100,
      adapt_delta = 0.90,
      max_treedepth = 12
    )
    return(fit)
  }
  
  
  fit <- tryCatch({
    fit = withTimeout(fitter(), timeout = 1000, cpu = 1000)
  }, error = function(e) {
    print("An error occurred!")
    return(NA)
  })
  
  
  
  if(is.environment(fit)){
    
    index = rnorm(1,0,1)
    
    hier_parameters = c("mu_b0","mu_b1","mu_b2","sd_b0","sd_b1","sd_b2","mu_b0_b","mu_b1_b","mu_b2_b","sd_b0_b","sd_b1_b","sd_b2_b","kappa_alpha","mu_alpha","sd_percept_precision")
    
    sub_parameters = c("b0","b1","b2","b0_b","b1_b","b2_b","alpha","percept_precision")
    
    return(list(hier = data.frame(fit$summary(hier_parameters),
                                  reals = hier %>% dplyr::select(all_of(hier_parameters)) %>% pivot_longer(everything()) %>% rename(reals = value),
                                  div = sum(fit$diagnostic_summary()$num_divergent),index = index),
                sub = data.frame(fit$summary(sub_parameters) %>% arrange(variable),
                                 reals = df %>% filter(trial == 1) %>% mutate(ids = 1:nrow(.)) %>% dplyr::select(all_of(sub_parameters),ids) %>% pivot_longer(cols = -ids) %>% rename(reals =value) %>% mutate(names = paste0(name,"[",ids,"]")) %>% arrange(names) %>% select(names,reals),
                                 div = sum(fit$diagnostic_summary()$num_divergent),index = index)))}else{
                                   return(list(NA,NA))
                                 }
}
