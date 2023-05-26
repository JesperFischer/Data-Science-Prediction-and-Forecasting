#visualize parameters


RW = function(dd){
  
  
  scripts = c("agent.R","parameter_recovery.R","stan_functions.R")
  q = lapply(here("scripts", scripts), source)
  alpha = seq(0.1,0.5,length.out = 3)
  w1 = seq(0.1,0.5,length.out = 3)
  #precision_percept = seq(10,100,50)
  #beta = seq(1,50,20)
  
  beta = 50
  precision_percept = seq(10,300,length.out = 3)
  
  parameters = expand.grid(alpha = alpha,w1 =  w1, precision_percept = precision_percept, beta = beta)
  
  
  
  parameters$id = 1:nrow(parameters)
  
  data_list <- split(parameters, parameters$id)
  
  plan(multisession, workers = 12)
  
  
  dd = furrr::future_map_dfr(data_list, ~our_rw_agent(.x, dd), .progress = TRUE, .options = furrr_options(seed = 123))
  
  dd$stim = as.factor(dd$stim)
  dd$precision_percept = as.factor(dd$precision_percept)
  
  
  levels(dd$stim) = c("Warm","Hot")
  
  
  plot1 = dd %>% 
    ggplot()+
    geom_line(aes(x = trial, y = association))+
    geom_point(aes(x = trial, y = u), col = "orange")+
    geom_line(aes(x = trial, y = desired))+
    facet_grid(alpha~w1, labeller = label_both)+
    theme_classic()+
    scale_y_continuous("Cue-stimulus association", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("Trial", breaks = scales::pretty_breaks(n = 3))    
  
  
  plot2 = dd %>% ggplot(aes(exp, y = percept,col = precision_percept, group = interaction(stim, precision_percept)))+
    geom_point(aes(shape = stim))+
    facet_grid(alpha~w1, labeller = label_both)+
    geom_smooth(method = "lm")+
    theme_classic()+
    scale_y_continuous("Mean of percept", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("% Expectation of Stimulus being Hot", breaks = scales::pretty_breaks(n = 3))+ theme(legend.position = "top")+ theme(text = element_text(size = 14))    
  
 
  plot1+plot2
   
}


KALMAN = function(){
  
  scripts = c("agent.R","parameter_recovery.R","stan_functions.R")
  q = lapply(here("scripts", scripts), source)
  #kalman filter:
  
  
  sigmaEta = seq(0.01,0.1,length.out = 5)
  sigmaPsi = seq(0.01,2,length.out = 5)
  sigmaEpsilon = seq(0.1,1,length.out = 5)
  
  
  parameters = expand.grid(sigmaEta = sigmaEta,sigmaPsi =  sigmaPsi, sigmaEpsilon = sigmaEpsilon)
  
  
  
  parameters$id = 1:nrow(parameters)
  
  data_list <- split(parameters, parameters$id)
  
  plan(multisession, workers = 12)
  
  dd = get_experiment()
  
  dd = furrr::future_map_dfr(data_list, ~our_kalman_agent(.x,dd), .progress = TRUE, .options = furrr_options(seed = 123))
  
  dd$sigmaEpsilon = as.factor(dd$sigmaEpsilon)
  dd$stim = as.factor(dd$stim)
  
  plot1 = dd %>% 
    ggplot()+
    geom_line(aes(x = trial, y = association, col = as.factor(sigmaEpsilon)))+
    geom_point(aes(x = trial, y = u), col = "blue")+
    geom_line(aes(x = trial, y = desired))+
    facet_grid(sigmaEta~sigmaPsi, labeller = label_both)+
    theme_classic()+
    scale_y_continuous("Cue-stimulus association", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("Trial", breaks = scales::pretty_breaks(n = 3))+ theme(legend.position = "none")    
  
  
  
  plot2 = dd %>% ggplot(aes(exp_mu, y = percept, color = sigmaEpsilon))+
    geom_point(aes(shape = stim))+
    geom_smooth(method = "lm")+
    facet_grid(sigmaEta~sigmaPsi, labeller = label_both)+
    theme_classic()+
    scale_y_continuous("Mean of percept", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("% belief of Stimulus being Hot", breaks = scales::pretty_breaks(n = 3))    
  
  
  plot1+plot2+plot_layout(guides = 'collect')
  
  
}

KALMAN_v2 = function(dd){

  scripts = c("agent.R","parameter_recovery.R","stan_functions.R")
  q = lapply(here("scripts", scripts), source)
  #kalman filter:
  
  
  sigmaEta = seq(0.01,1,length.out = 3)
  sigmaPsi = seq(0.01,5,length.out = 3)
  sigmaEpsilon = seq(0.1,5,length.out = 3)
  beta = 50
  precision_percept = 155
  
  
  parameters = expand.grid(sigmaEta = sigmaEta,sigmaPsi =  sigmaPsi, sigmaEpsilon = sigmaEpsilon, beta = beta, precision_percept = precision_percept)
  
  
  
  parameters$id = 1:nrow(parameters)
  
  data_list <- split(parameters, parameters$id)
  
  plan(multisession, workers = 12)

  
  dd = furrr::future_map_dfr(data_list, ~our_kalman_agent22(.x,dd), .progress = TRUE, .options = furrr_options(seed = 123))
  
  dd$sigmaEpsilon = as.factor(dd$sigmaEpsilon)
  dd$stim = as.factor(dd$stim)
  
  levels(dd$stim) = c("Warm","Hot")
  
  plot1 = dd %>% 
    ggplot()+
    geom_line(aes(x = trial, y = association, col = sigmaEpsilon))+
    geom_point(aes(x = trial, y = u), col = "orange")+
    geom_line(aes(x = trial, y = desired))+
    facet_grid(sigmaEta~sigmaPsi, labeller = label_both)+
    theme_classic()+
    scale_y_continuous("Cue-stimulus association", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("Trial", breaks = scales::pretty_breaks(n = 3))+theme(legend.position = "top")    
  
  
  
  plot2 = dd %>% 
    ggplot(aes(exp_mu, y = percept, color = sigmaEpsilon, group = interaction(stim, sigmaEpsilon)))+
    geom_point(aes(shape = stim))+
    geom_smooth(method = "lm")+
    facet_grid(sigmaEta~sigmaPsi, labeller = label_both)+
    theme_classic()+
    scale_y_continuous("Mean of percept", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("% Expectation of Stimulus being Hot", breaks = scales::pretty_breaks(n = 3))+theme(legend.position = "top")+
    ggtitle("precision percept = 155")    
  
  
  plot1+plot2
  
  
}

WB = function(dd){
  
  
  scripts = c("agent.R","parameter_recovery.R","stan_functions.R")
  q = lapply(here("scripts", scripts), source)
  alpha = seq(0.1,0.5,length.out = 3)
  w1 = seq(0.1,0.5,length.out = 3)
  w2 = seq(0.1,0.5,length.out = 3)
  percept_precision = seq(10,300,length.out = 3)
  beta = 50
  
  
  
  parameters = expand.grid(alpha = alpha, w1 =  w1, w2 = w2, percept_precision = percept_precision, beta = beta)
  
  
  
  parameters$id = 1:nrow(parameters)
  
  data_list <- split(parameters, parameters$id)
  
  plan(multisession, workers = 6)
  
  
  dd = furrr::future_map_dfr(data_list, ~weighted_Bayes_f(.x,dd), .progress = TRUE, .options = furrr_options(seed = 123))
  
  
  dd$stim = ifelse(dd$stim > 0.5, 1, 0)
  
  dd$alpha = as.factor(dd$alpha)
  dd$stim = as.factor(dd$stim)
  
  levels(dd$stim) = c("Warm","Hot")
  
  plot1 = dd %>% filter(percept_precision == 10) %>% 
    ggplot()+
    geom_line(aes(x = trial, y = association, col = alpha))+
    geom_point(aes(x = trial, y = u), col = "orange")+
    geom_line(aes(x = trial, y = desired))+
    facet_grid(w1~w2, labeller = label_both)+
    theme_classic()+
    scale_y_continuous("Cue-stimulus association", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("Trial", breaks = scales::pretty_breaks(n = 3))+theme(legend.position = "top")    
  
  
  dd$percept_precision = as.factor(dd$percept_precision)
  
  plot2 = dd %>% filter(alpha == 0.3) %>%  ggplot(aes(expect, y = percept, color = percept_precision, group = interaction(stim, percept_precision)))+
    geom_point(aes(shape = stim))+
    geom_smooth(method = "lm")+
    facet_grid(w1~w2, labeller = label_both)+
    theme_classic()+
    scale_y_continuous("Mean of percept", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("% Expectation of Stimulus being Hot", breaks = scales::pretty_breaks(n = 3))+theme(legend.position = "top")+
    ggtitle("alpha = 0.3")    
  
  
  
  plot1+plot2
  
  
}

RW_logistic_vis = function(dd){
  
  
  scripts = c("agent.R","parameter_recovery.R","stan_functions.R")
  q = lapply(here("scripts", scripts), source)
  
  alpha = seq(0.1,0.5,length.out = 3)
  b0 = -1.5
  b1 = seq(0.1,2,length.out = 3)
  b2 = seq(0.1,2,length.out = 3)
  
  b0_b = -4
  b1_b = 7.8
  b2_b = 0.7
  
  
  #precision_percept = seq(10,100,50)
  #beta = seq(1,50,20)
  
  percept_precision = seq(10,300,length.out = 3)
  
  parameters = expand.grid(alpha = alpha,b0=b0,b1=b1,b2=b2,b0_b=b0_b,b1_b=b1_b,b2_b=b2_b, percept_precision = percept_precision)
  
  
  
  parameters$id = 1:nrow(parameters)
  
  data_list <- split(parameters, parameters$id)
  
  plan(multisession, workers = 12)
  
  dd = get_experiment()
  
  dd = furrr::future_map_dfr(data_list, ~rw_logistic_v2(.x, dd), .progress = TRUE, .options = furrr_options(seed = 123))
  
  dd$stim = as.factor(dd$stim)
  dd$percept_precision = as.factor(dd$percept_precision)
  
  
  levels(dd$stim) = c("Warm","Hot")
  
  
  plot1 = dd %>% 
    ggplot()+
    geom_line(aes(x = trial, y = belief))+
    geom_point(aes(x = trial, y = u), col = "orange")+
    geom_line(aes(x = trial, y = desired))+
    facet_grid(alpha~b2, labeller = label_both)+
    theme_classic()+
    scale_y_continuous("Cue-stimulus association", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("Trial", breaks = scales::pretty_breaks(n = 3))    
  
  
  plot2 = dd %>% filter(alpha == 0.1) %>% ggplot(aes(expect, y = percept,col = percept_precision, group = interaction(stim, percept_precision)))+
    geom_point(aes(shape = stim))+
    facet_grid(b1~b2, labeller = label_both)+
    geom_smooth(method = "lm")+
    theme_classic()+
    scale_y_continuous("Mean of percept", breaks = scales::pretty_breaks(n = 3))+
    scale_x_continuous("% Expectation of Stimulus being Hot", breaks = scales::pretty_breaks(n = 3))+ theme(legend.position = "top")+
    ggtitle("Alpha = 0.1")    
  
  
  plot1+plot2
  
}



plot_parameterrecovery = function(pp, div, subs){
  
  if(div){
  hier_pp <- do.call(rbind, lapply(pp, "[[", 1)) %>% drop_na()
  sub_pp <- do.call(rbind, lapply(pp, "[[", 2))%>% drop_na()
  }else{
    hier_pp <- do.call(rbind, lapply(pp, "[[", 1)) %>% filter(div < 1)%>% drop_na()
    sub_pp <- do.call(rbind, lapply(pp, "[[", 2)) %>% filter(div < 1)%>% drop_na()
  }
  
  hier = hier_pp %>% mutate(divergences = ifelse(div > 0,"yes","no")) %>% 
    ggplot(aes(x = mean, y = reals.reals))+
    geom_point()+
    {if(div)geom_point(aes(x = mean, y = reals.reals, col = as.factor(divergences)))}+
    facet_wrap(~variable, scales = "free")+
    theme_classic()+
    geom_abline(slope = 1, intercept = 0)+
    scale_color_manual("divergence",values = c("black","red"))+
    scale_y_continuous("simulated value",breaks = scales::pretty_breaks(n = 4))+
    scale_x_continuous("Recovered value", breaks = scales::pretty_breaks(n = 4))
  

  # Generate the regular expression pattern dynamically
  pattern <- paste0("\\[", 1:subs, "\\]")
  
  # Filter the column
  sub_pp <- sub_pp[grepl(paste(pattern, collapse = "|"), sub_pp$variable), ]
  
  
  
  sub = sub_pp %>% mutate(divergences = ifelse(div > 0,"yes","no")) %>% 
    ggplot(aes(x = mean, y = reals.reals))+
    geom_point()+
    {if(div)geom_point(aes(col = as.factor(divergences)))}+
    facet_wrap(~variable, scales = "free", ncol = subs)+
    theme_classic()+
    geom_abline(slope = 1, intercept = 0)+
    scale_color_manual("divergence",values = c("black","red"))+
    scale_y_continuous("simulated value",breaks = scales::pretty_breaks(n = 4))+
    scale_x_continuous("Recovered value", breaks = scales::pretty_breaks(n = 4))
  
  
  
  div = hier_pp %>% group_by(index) %>%
    slice(1) %>%
    ungroup() %>% ggplot(aes(x = div))+geom_histogram()+theme_classic()+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
    scale_x_continuous("Number of divergences", breaks = scales::pretty_breaks(n = 4))
  
  
  
  
  return(list(hier = hier,sub = sub, div = div))
  
}


get_corplot = function(fit1,fit2,variable){


  df1 = (fit1$summary(variables = variable))
  df2 = (fit2$summary(variables = variable))
  
  df1$sess = df2$mean
  
  df1$sess_sd = df2$sd
  
  
  #fit <- lsfit(x = df1$sess, y = df1$mean, weights = list(x = df1$sess_sd, y = df1$sd))
  
  # Extract the correlation coefficient
  #correlation <- cor(fit$residuals$x, fit$residuals$y)
  

  plot = df1 %>% ggplot(aes(x = mean, y = sess)) +
    geom_point() +
    geom_errorbar(aes(ymin = sess - sess_sd, ymax = sess + sess_sd), width = 0.005,linetype = "dashed") +
    geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), height = 0.005,linetype = "dashed") +
    labs(x = "Session 1", y = "Session 2") +
    ggtext::geom_richtext(aes(x = min(mean)+sd(mean), y = max(sess),
                              label = paste0("r = ",round(cor.test(df1$mean, df1$sess)$estimate[[1]],2), " [", round(cor.test(df1$mean, df1$sess)$conf.int[[1]],2)," ; ", round(cor.test(df1$mean, df1$sess)$conf.int[[2]],2),"]")))+ggtitle(paste0("test retest correlation of ",variable))+
    theme_classic()
  
  return(plot)
  
}