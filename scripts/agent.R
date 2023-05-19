
get_experiment = function(){
  trials = c(40,40,40,40)
  bias = rep(c(0.1,0.9,0.2,0.8),trials)
  
  
  cue = rbinom(sum(trials),1,0.5)
  
  stim = array(sum(trials, NA))
  for(i in 1:sum(trials)){
    stim[i] = ifelse(cue[i] == 1, rbinom(1,1,bias[i]), rbinom(1,1,(1-bias[i])))
    
  }
  
  u = ifelse(cue == stim, 1,0)
  dd = data.frame(stim = stim, cue = cue, u = u, bias = bias)
  
  return(dd)
  
  
}




our_rw_agent = function(parameters){
  

  dd = get_experiment()
  
  cue = dd$cue
  stim = dd$stim
  u = dd$u
  bias = dd$bias
  
  
  
  ntrials = nrow(dd)
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  exp = array(NA, ntrials)
  
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  pe = array(NA,ntrials)
  percept_bin = array(NA,ntrials)
  
  
  w1 = parameters$w1
  alpha = parameters$alpha
  precision_percept = parameters$precision_percept
  beta = parameters$beta
  
  
  association[1] = 0.5
  

  for (i in seq(ntrials)){
    
    #get the expectation of the first trial i.e. (how likeli is hot?)
    exp[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    #the prediction
    pred[i] = rbinom(1,1,(exp[i]^beta)/((exp[i]^beta)+(1-exp[i])^(beta)))
    
    #the mean percept given from the stimulus received and the expectation
    perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*exp[i] == 0, 0.01,
                          ifelse(w1*stim[i]+(1-w1)*exp[i] == 1, 0.99, w1*stim[i]+(1-w1)*exp[i]))
    
    #generating the percept
    percept[i] = extraDistr::rprop(1,precision_percept,perceptmu[i])
    
    #now we get the binary percept (did you feel hot or warm?)
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    #prediction error for the association
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -exp[i]),(-(perceptmu[i] -exp[i])))  
    
    #association update
    association[i+1] = association[i]+alpha*pe[i]
    
  }
  
  df = data.frame(u = u[1:ntrials], stim = stim[1:ntrials], percept = percept[1:ntrials], pred = pred[1:ntrials],association = association[1:ntrials],
                  exp = exp[1:ntrials],pe = pe[1:ntrials],cue = cue[1:ntrials], alpha = alpha,
                  percept_bin = percept_bin[1:ntrials],
                  w1 = w1,precision_percept = precision_percept,
                  trial = 1:ntrials, perceptmu = perceptmu[1:ntrials], desired = rep(bias,1), id = rnorm(1,0,1))
  
  return(df)

  
  
}

our_rw_agent_v2 = function(w1,alpha,precision_percept,beta){
  
  
  dd = get_experiment()
  
  cue = dd$cue
  stim = dd$stim
  u = dd$u
  bias = dd$bias
  
  
  
  ntrials = nrow(dd)
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  exp = array(NA, ntrials)
  
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  pe = array(NA,ntrials)
  percept_bin = array(NA,ntrials)
  
  
  association[1] = 0.5
  
  
  for (i in seq(ntrials)){
    
    #get the expectation of the first trial i.e. (how likeli is hot?)
    exp[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    #the prediction
    pred[i] = rbinom(1,1,(exp[i]^beta)/((exp[i]^beta)+(1-exp[i])^(beta)))
    
    #the mean percept given from the stimulus received and the expectation
    perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*exp[i] == 0, 0.01,
                          ifelse(w1*stim[i]+(1-w1)*exp[i] == 1, 0.99, w1*stim[i]+(1-w1)*exp[i]))
    
    #generating the percept
    percept[i] = extraDistr::rprop(1,precision_percept,perceptmu[i])
    
    #now we get the binary percept (did you feel hot or warm?)
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    #prediction error for the association
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -exp[i]),(-(perceptmu[i] -exp[i])))  
    
    #association update
    association[i+1] = association[i]+alpha*pe[i]
    
  }
  
  df = data.frame(u = u[1:ntrials], stim = stim[1:ntrials], percept = percept[1:ntrials], pred = pred[1:ntrials],association = association[1:ntrials],
                  exp = exp[1:ntrials],pe = pe[1:ntrials],cue = cue[1:ntrials], alpha = alpha,
                  percept_bin = percept_bin[1:ntrials],
                  w1 = w1,
                  precision_percept = precision_percept,
                  beta = beta,
                  trial = 1:ntrials, perceptmu = perceptmu[1:ntrials], desired = rep(bias,1), id = rnorm(1,0,1))
  
  return(df)
  
  
  
}


our_hier_rw_agent = function(parameters){
  
  nsubs = parameters$nsubs
  
  df = data.frame()
  for (s in 1:nsubs){
    
    df1 = our_rw_agent_v2(w1 = extraDistr::rprop(1,parameters$kappa_w1,parameters$mu_w1),
                 alpha = extraDistr::rprop(1,parameters$kappa_alpha,parameters$mu_alpha),
                 precision_percept = rexp(1,parameters$sd_precision_percept),
                 beta = rexp(1,parameters$sd_beta)
                 )
    
    df = rbind(df,df1)
    
  }
  
  
  return(list(df, parameters))
  
  
  
}






our_kalman_agent = function(parameters){
  
  dd = get_experiment()
  
  cue = dd$cue
  stim = dd$stim
  u = dd$u
  bias = dd$bias
  
  
  
  ntrials = nrow(dd)
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  perceptvar = array(NA, ntrials)
  perceptprec = array(NA,ntrials)
  
  
  
  exp_mu = array(NA, ntrials)
  exp_var = array(NA, ntrials)
  exp_prec = array(NA, ntrials)  
  
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  pe = array(NA,ntrials)
  percept_bin = array(NA,ntrials)
  
  association = array(NA,ntrials)
  
  association[1] = 0.5
  exp_var[1] = 2
  
  sigmaEta = parameters$sigmaEta
  sigmaEpsilon = parameters$sigmaEpsilon
  sigmaPsi = parameters$sigmaPsi
  
  for(t in 1:ntrials){
    
    exp_mu[t] = ifelse(cue[t] == 1, association[t], 1-association[t])
    
    exp_prec[t] = 1/exp_var[t]
    
    exp_mu[t] = ifelse(exp_mu[t] > 0.999, 0.999, ifelse(exp_mu[t] < 0.001, 0.001, exp_mu[t]))
    
    pred[t] = rbinom(1,1,(exp_mu[t]^exp_prec[t])/((exp_mu[t]^exp_prec[t])+(1-exp_mu[t])^(exp_prec[t])))
    
    
    perceptmu[t] =  (sigmaEpsilon * exp_mu[t] + (sigmaPsi + exp_var[t]) * stim[t] ) / 
      (sigmaEpsilon + sigmaPsi + exp_var[t]);
    
    perceptvar[t] = ( sigmaEpsilon * (sigmaPsi + exp_var[t]) ) / (sigmaEpsilon + sigmaPsi + exp_var[t])
    
    perceptprec[t]= 1/perceptvar[t]
    
    perceptmu[t] = ifelse(perceptmu[t] > 0.999, 0.999, ifelse(perceptmu[t] < 0.001, 0.001, perceptmu[t]))
    
    percept[t] = extraDistr::rprop(1,perceptprec[t],perceptmu[t])
    
    
    percept_bin[t] = rbinom(1,1,(perceptmu[t]^perceptprec[t])/((perceptmu[t]^perceptprec[t])+(1-perceptmu[t])^(perceptprec[t])))
    
    association[t+1] <- ((sigmaEpsilon + sigmaPsi) * association[t] + (exp_var[t] * u[t])) / 
      (sigmaEpsilon + sigmaPsi + exp_var[t])
    
    exp_var[t+1] <- ((sigmaEpsilon + sigmaPsi) * exp_var[t] / (sigmaEpsilon + sigmaPsi + exp_var[t])) + sigmaEta
    
  }
  
  
  return(data.frame(exp_mu = exp_mu[1:ntrials], perceptmu = perceptmu, pred = pred, percept = percept,
                    percept_bin = percept_bin, trial = 1:ntrials, stim = stim, u = u, association = association[1:ntrials], cue = cue,
                    sigmaPsi = sigmaPsi, sigmaEta = sigmaEta, sigmaEpsilon = sigmaEpsilon))
  
  
  
  
}


our_kalman_agent_v2 = function(sigmaEta, sigmaEpsilon, sigmaPsi){
  
  dd = get_experiment()
  
  cue = dd$cue
  stim = dd$stim
  u = dd$u
  bias = dd$bias
  
  
  
  ntrials = nrow(dd)
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  perceptvar = array(NA, ntrials)
  perceptprec = array(NA,ntrials)
  
  
  
  exp_mu = array(NA, ntrials)
  exp_var = array(NA, ntrials)
  exp_prec = array(NA, ntrials)  
  
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  pe = array(NA,ntrials)
  percept_bin = array(NA,ntrials)
  
  association = array(NA,ntrials)
  
  association[1] = 0.5
  exp_var[1] = 2
  
  for(t in 1:ntrials){
    
    exp_mu[t] = ifelse(cue[t] == 1, association[t], 1-association[t])
    
    exp_prec[t] = 1/exp_var[t]
    
    exp_mu[t] = ifelse(exp_mu[t] > 0.999, 0.999, ifelse(exp_mu[t] < 0.001, 0.001, exp_mu[t]))
    
    pred[t] = rbinom(1,1,(exp_mu[t]^exp_prec[t])/((exp_mu[t]^exp_prec[t])+(1-exp_mu[t])^(exp_prec[t])))
    
    
    perceptmu[t] =  (sigmaEpsilon * exp_mu[t] + (sigmaPsi + exp_var[t]) * stim[t] ) / 
      (sigmaEpsilon + sigmaPsi + exp_var[t]);
    
    perceptvar[t] = ( sigmaEpsilon * (sigmaPsi + exp_var[t]) ) / (sigmaEpsilon + sigmaPsi + exp_var[t])
    
    perceptprec[t]= 1/perceptvar[t]
    
    perceptmu[t] = ifelse(perceptmu[t] > 0.999, 0.999, ifelse(perceptmu[t] < 0.001, 0.001, perceptmu[t]))
    
    
    percept[t] = extraDistr::rprop(1,perceptprec[t],perceptmu[t])
    
    
    percept_bin[t] = rbinom(1,1,(perceptmu[t]^perceptprec[t])/((perceptmu[t]^perceptprec[t])+(1-perceptmu[t])^(perceptprec[t])))
    
    association[t+1] <- ((sigmaEpsilon + sigmaPsi) * association[t] + (exp_var[t] * u[t])) / 
      (sigmaEpsilon + sigmaPsi + exp_var[t])
    
    exp_var[t+1] <- ((sigmaEpsilon + sigmaPsi) * exp_var[t] / (sigmaEpsilon + sigmaPsi + exp_var[t])) + sigmaEta
    
    
    
  }
  
  
  return(data.frame(exp_mu = exp_mu[1:ntrials], perceptmu = perceptmu, pred = pred, percept = percept,
                    percept_bin = percept_bin, trial = 1:ntrials, stim = stim, u = u, association = association[1:ntrials], cue = cue, id = rnorm(1,0,1),
                    sigmaPsi = sigmaPsi, sigmaEta = sigmaEta, sigmaEpsilon = sigmaEpsilon))
  
  
  
  
}


our_hier_kalman_agent = function(parameters){
  
  nsubs = parameters$nsubs
  
  df = data.frame()
  for (s in 1:nsubs){
    
    df1 = our_kalman_agent_v2(sigmaEta = rlnorm(1,parameters$mu_sigmaEta,parameters$sd_sigmaEta),
                          sigmaEpsilon = rlnorm(1,parameters$mu_sigmaEpsilon,parameters$sd_sigmaEpsilon),
                          sigmaPsi = rlnorm(1,parameters$mu_sigmaPsi,parameters$sd_sigmaPsi)
    )
    
    df = rbind(df,df1)
    
  }
  
  
  return(list(df, parameters))
  
  
  
}



weighted_Bayes_f = function(parameters) {
  
  trials = c(40,40,40,40)
  bias = rep(c(0.1,0.9,0.1,0.9),trials)
  ntrials = sum(trials)
  
  cue = rbinom(sum(trials),1,0.5)
  
  stim = array(sum(trials, NA))
  for(i in 1:ntrials){
    stim[i] = ifelse(cue[i] == 1, rbinom(1,1,bias[i]), rbinom(1,1,(1-bias[i])))
    
  }
  
  u = ifelse(cue == stim, 1,0)
  dd = data.frame(stim = stim, cue = cue, u = u)
  
  dd %>% ggplot(aes(x = 1:nrow(.), y = u, col = as.factor(stim)))+geom_point()+theme_classic()
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  percept_bin = array(NA, ntrials)
  expect = array(NA, ntrials)
  pe = array(NA, ntrials)
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  
  association[1] = 0.5
  
  alpha = parameters$alpha
  w1 = parameters$w1
  w2 = parameters$w2
  percept_precision = parameters$percept_precision
  beta = parameters$beta
  
  
  for (i in seq(ntrials)){
    
    expect[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    
    pred[i] = rbinom(1,1,(expect[i]^beta)/((expect[i]^beta)+(1-expect[i])^(beta)))
    
    # correlated weight / Resterla Wagner paradigm
    #perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*expect[i] == 0, 0.01, ifelse(w1*stim[i]+(1-w1)*expect[i] == 1, 0.99,w1*stim[i]+(1-w1)*expect[i]))
    
    # Bayesian inference paradigm: weighted Bayes
    stim[i] = ifelse(stim[i] == 0, 0.01, ifelse(stim[i] == 1, 0.99, stim[i]))
    expect[i] = ifelse(expect[i] == 0, 0.01, ifelse(expect[i] == 1, 0.99, expect[i]))
    perceptmu[i] = inv_logit_scaled(w1 * logit_scaled(stim[i]) + w2 * logit_scaled(expect[i]))
    
    # agent's rating on the intensity of the stimulus
    percept[i] = extraDistr::rprop(1,percept_precision,perceptmu[i])
    # agent's response on burning or not burning
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -expect[i]),(-(perceptmu[i] -expect[i])))
    association[i+1] = association[i]+alpha*pe[i]
    
  }
  
  df = data.frame(u = u[1:ntrials], cue = cue[1:ntrials], stim = stim[1:ntrials], expect = expect[1:ntrials], pred = pred[1:ntrials], perceptmu = perceptmu[1:ntrials], percept_bin = percept_bin[1:ntrials], percept = percept[1:ntrials], association = association[1:ntrials],
                  pe = pe[1:ntrials], alpha = alpha, w1 = w1, w2 = w2, percept_precision = percept_precision, beta = beta, x = 1:ntrials, desired = rep(bias,1))
  
  return(df)
  
}



weighted_Bayes_f_v2 = function(alpha, w1, w2, percept_precision, beta) {
  
  dd = get_experiment()
  
  cue = dd$cue
  stim = dd$stim
  u = dd$u
  bias = dd$bias
  
  ntrials = nrow(dd)
  
  percept = array(NA, ntrials)
  perceptmu = array(NA, ntrials)
  percept_bin = array(NA, ntrials)
  expect = array(NA, ntrials)
  pe = array(NA, ntrials)
  association = array(NA, ntrials)
  pred = array(NA,ntrials)
  
  association[1] = 0.5
  
  for (i in seq(ntrials)){
    
    expect[i] = ifelse(cue[i] == 1, association[i], 1-association[i])
    
    pred[i] = rbinom(1,1,(expect[i]^beta)/((expect[i]^beta)+(1-expect[i])^(beta)))
    
    # correlated weight / Resterla Wagner paradigm
    #perceptmu[i] = ifelse(w1*stim[i]+(1-w1)*expect[i] == 0, 0.01, ifelse(w1*stim[i]+(1-w1)*expect[i] == 1, 0.99,w1*stim[i]+(1-w1)*expect[i]))
    
    # Bayesian inference paradigm: weighted Bayes 
    stim[i] = ifelse(stim[i] == 0, 0.01, ifelse(stim[i] == 1, 0.99, stim[i]))
    expect[i] = ifelse(expect[i] == 0, 0.01, ifelse(expect[i] == 1, 0.99, expect[i]))
    perceptmu[i] = inv_logit_scaled(w1 * logit_scaled(stim[i]) + w2 * logit_scaled(expect[i]))
    
    # agent's rating on the intensity of the stimulus
    percept[i] = extraDistr::rprop(1,percept_precision,perceptmu[i])
    # agent's response on burning or not burning
    percept_bin[i] = rbinom(1,1,(perceptmu[i]^beta)/((perceptmu[i]^beta)+(1-perceptmu[i])^(beta)))
    
    pe[i] = ifelse(cue[i] == 1, (perceptmu[i] -expect[i]),(-(perceptmu[i] -expect[i])))
    association[i+1] = association[i]+alpha*pe[i]
    
  }
  
  df = data.frame(u = u[1:ntrials], cue = cue[1:ntrials], stim = stim[1:ntrials], expect = expect[1:ntrials], pred = pred[1:ntrials], perceptmu = perceptmu[1:ntrials], percept_bin = percept_bin[1:ntrials], percept = percept[1:ntrials], association = association[1:ntrials],
                  pe = pe[1:ntrials], alpha = alpha, w1 = w1, w2 = w2, percept_precision = percept_precision, beta = beta, x = 1:ntrials, desired = rep(bias,1), id = rnorm(1,0,1))
  
  return(df)
  
}

hier_Bayes_f = function(parameters){
  
  nsubs = parameters$nsubs
  
  df = data.frame()
  for (s in 1:nsubs){
    
    df1 = weighted_Bayes_f_v2(w1 = extraDistr::rprop(1,parameters$kappa_w1,parameters$mu_w1),
                              w2 = extraDistr::rprop(1,parameters$kappa_w2,parameters$mu_w2),    
                              alpha = extraDistr::rprop(1,parameters$kappa_alpha,parameters$mu_alpha),
                              percept_precision = rexp(1,parameters$sd_percept_precision),
                              beta = rexp(1,parameters$sd_beta)
    )
    
    df = rbind(df,df1)
    
  }
  
  
  return(list(df, parameters))
  
  
  
}
