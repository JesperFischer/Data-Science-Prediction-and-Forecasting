#parameter_recovery


parameter_recovery = function(parameters){
  
  df = our_hier_rw_agent(parameters)[[1]]
  hier = our_hier_rw_agent(parameters)[[2]]
  
  
  source(here::here("stan_functions.R"))
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
  
  mod = cmdstan_model(here::here("myRW_real.stan"))
  
  fit <- mod$sample(
    data = data1,
    chains = 4, 
    parallel_chains = 4,
    refresh = 100,
    adapt_delta = 0.80,
    max_treedepth = 10
  )
  
  index = rnorm(1,0,1)
  
  hier_parameters = c("kappa_w1","mu_w1","kappa_alpha","mu_alpha","sd_beta","sd_precision_percept")
  
  sub_parameters = c("w1","alpha","beta","precision_percept")
  return(list(hier = data.frame(fit$summary(hier_parameters),
                    reals = hier %>% dplyr::select(all_of(hier_parameters)) %>% pivot_longer(everything()) %>% rename(reals = value),
                    div = sum(fit$diagnostic_summary()$num_divergent),index = index),
              sub = data.frame(fit$summary(sub_parameters) %>% arrange(variable),
                               reals = df %>% filter(trial == 1) %>% mutate(ids = 1:nrow(.)) %>% dplyr::select(all_of(sub_parameters),ids) %>% pivot_longer(cols = -ids) %>% rename(reals =value) %>% mutate(names = paste0(name,"[",ids,"]")) %>% arrange(names) %>% select(names,reals),
                               div = sum(fit$diagnostic_summary()$num_divergent),index = index)))
}



