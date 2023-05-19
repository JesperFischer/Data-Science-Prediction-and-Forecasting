data {
  
  int<lower=1> nsubs; // number of subjects
  int<lower=1> ntrials; // number of subjects
  matrix[ntrials, nsubs] percept; // observations
  int expectPain[ntrials, nsubs]; // prediciton
  int percept_bin[ntrials, nsubs]; // prediciton
  
  matrix[ntrials, nsubs] stim; // observations
  matrix[ntrials, nsubs] cues; // observations
  matrix[ntrials, nsubs] u; // observations
  
  
}

parameters {
  vector<lower=0>[nsubs] sigmaEta;
  vector<lower=0>[nsubs] sigmaPsi;
  vector<lower=0>[nsubs] sigmaEpsilon;
  
  // Group-level parameters
  real  mu_sigmaEpsilon;
  real  mu_sigmaEta;
  real  mu_sigmaPsi;
  
  real <lower=0> sd_sigmaEpsilon;
  real <lower=0> sd_sigmaEta;
  real <lower=0> sd_sigmaPsi;
  
}

transformed parameters {
  matrix <lower=0, upper  = 1> [ntrials, nsubs] perceptmu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] perceptmuu;
  matrix <lower=0> [ntrials, nsubs] perceptvar;
  
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] association; 
  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] exp_mu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] exp_muu;
  matrix <lower=0> [ntrials+1, nsubs] exp_var;
  

  for (s in 1:nsubs) {
      association[1,s] = 0.5;
      exp_var[1,s] = sigmaEta[s];
    for (t in 1:(ntrials)){
      
      if(cues[t,s] == 1){
        exp_muu[t,s] = association[t,s];
      }else{
        exp_muu[t,s] = 1-association[t,s];
      }
      
      
      perceptmuu[t,s] =  (sigmaEpsilon[s] * exp_muu[t,s] + (sigmaPsi[s] + exp_var[t,s]) * stim[t,s] ) / 
                        (sigmaEpsilon[s] + sigmaPsi[s] + exp_var[t,s]);
                       
                       
                       
      perceptvar[t,s] = (sigmaEpsilon[s] * (sigmaPsi[s] + exp_var[t,s]) ) / 
                        (sigmaEpsilon[s] + sigmaPsi[s] + exp_var[t,s]) ;    
                        
                        
      
      association[t+1,s] = ((sigmaEpsilon[s] + sigmaPsi[s]) * exp_muu[t,s] + (exp_var[t,s] * u[t,s])) / 
                                       (sigmaEpsilon[s] + sigmaPsi[s] + exp_var[t,s]) ;
                                       
      exp_var[t+1,s] = ((sigmaEpsilon[s] + sigmaPsi[s]) * exp_var[t,s] / (sigmaEpsilon[s] + sigmaPsi[s] + exp_var[t,s])) + sigmaEta[s];
      
      if (exp_muu[t,s] > 0.9999){
                exp_mu[t,s] = 0.9999;
            } else if (exp_muu[t,s] < 0.0001) {
                exp_mu[t,s] = 0.0001;
            } else if (exp_muu[t,s] > 0.0001 && exp_muu[t,s] < 0.9999) {
                exp_mu[t,s] = exp_muu[t,s];
            } else {
                exp_mu[t,s] = 0.5;
            }
      
      if (perceptmuu[t,s] > 0.9999){
                perceptmu[t,s] = 0.9999;
            } else if (perceptmuu[t,s] < 0.0001) {
                perceptmu[t,s] = 0.0001;
            } else if (perceptmuu[t,s] > 0.0001 && perceptmuu[t,s] < 0.9999) {
                perceptmu[t,s] = perceptmuu[t,s];
            } else {
                perceptmu[t,s] = 0.5;
            }
      
      
    }
  }
}

model {
   for (s in 1:nsubs){

    target += lognormal_lpdf(sigmaEta[s] | mu_sigmaEta, sd_sigmaEta);
    target += lognormal_lpdf(sigmaPsi[s] | mu_sigmaPsi, sd_sigmaPsi);
    target += lognormal_lpdf(sigmaEpsilon[s] | mu_sigmaEpsilon, sd_sigmaEpsilon);
    
    
    
    for (t in 1:ntrials){
      target += beta_proportion_lpdf(percept[t,s] | perceptmu[t,s], 1/perceptvar[t,s]);
      
      target += bernoulli_lpmf(percept_bin[t,s] | (perceptmu[t,s]^(1/perceptvar[t,s]))/((perceptmu[t,s]^(1/perceptvar[t,s]))+(1-perceptmu[t,s])^(1/perceptvar[t,s])));

      target += bernoulli_lpmf(expectPain[t,s] |  (exp_mu[t,s]^(1/exp_var[t,s]))/((exp_mu[t,s]^(1/exp_var[t,s]))+(1-exp_mu[t,s])^(1/exp_var[t,s])));
  
    }
    
  }
  
  // Hierarchical Priors
  
  target += normal_lpdf(mu_sigmaEta | 0 , 2);
  target += normal_lpdf(mu_sigmaPsi | 0 , 2);
  target += normal_lpdf(mu_sigmaEpsilon | 0 , 2);
  
  target += lognormal_lpdf(sd_sigmaEta | 0 , 1);
  target += lognormal_lpdf(sd_sigmaPsi | 0 , 1);
  target += lognormal_lpdf(sd_sigmaEpsilon | 0 , 1);
  
}



generated quantities{
  vector<lower=0>[nsubs] prior_sigmaEta;
  vector<lower=0>[nsubs] prior_sigmaPsi;
  vector<lower=0>[nsubs] prior_sigmaEpsilon;
  
  // Group-level parameters
  real  prior_mu_sigmaEpsilon;
  real  prior_mu_sigmaEta;
  real  prior_mu_sigmaPsi;
  
  real <lower=0> prior_sd_sigmaEpsilon;
  real <lower=0> prior_sd_sigmaEta;
  real <lower=0> prior_sd_sigmaPsi;
  
  //trial level
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_perceptmu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_perceptmuu;
  matrix <lower=0> [ntrials, nsubs] prior_perceptvar;
  
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] prior_association; 
  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_exp_mu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_exp_muu;
  matrix <lower=0> [ntrials+1, nsubs] prior_exp_var;
  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_percept;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_percept_bin;  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_expectPain;
  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] post_percept;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] post_percept_bin;  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] post_expectPain;
  
  
  
  matrix[ntrials, nsubs] log_lik;

    
  
  
  prior_mu_sigmaEta = normal_rng(0 , 2);
  prior_mu_sigmaPsi = normal_rng(0 , 2);
  prior_mu_sigmaEpsilon = normal_rng(0 , 2);
  
  prior_sd_sigmaEta = lognormal_rng(0 , 1);
  prior_sd_sigmaPsi = lognormal_rng(0 , 1);
  prior_sd_sigmaEpsilon = lognormal_rng(0 , 1);
  
  
  
  
  for (s in 1:nsubs) {
      prior_sigmaEta[s] = lognormal_rng(prior_mu_sigmaEta, prior_sd_sigmaEta);
      prior_sigmaPsi[s] = lognormal_rng(prior_mu_sigmaPsi, prior_sd_sigmaPsi);
      prior_sigmaEpsilon[s] = lognormal_rng(prior_mu_sigmaEpsilon, prior_sd_sigmaEpsilon);
    
    
    
      prior_association[1,s] = 0.5;
      prior_exp_var[1,s] = prior_sigmaEta[s];
    for (t in 1:(ntrials)){
      
      if(cues[t,s] == 1){
        prior_exp_muu[t,s] = prior_association[t,s];
      }else{
        prior_exp_muu[t,s] = 1-prior_association[t,s];
      }
      
      
      prior_perceptmuu[t,s] =  (prior_sigmaEpsilon[s] * prior_exp_muu[t,s] + (prior_sigmaPsi[s] + prior_exp_var[t,s]) * stim[t,s] ) / 
                        (prior_sigmaEpsilon[s] + prior_sigmaPsi[s] + prior_exp_var[t,s]);
                       
                       
                       
      prior_perceptvar[t,s] = (prior_sigmaEpsilon[s] * (prior_sigmaPsi[s] + prior_exp_var[t,s]) ) / 
                        (prior_sigmaEpsilon[s] + prior_sigmaPsi[s] + prior_exp_var[t,s]) ;    
                        
                        
      
      prior_association[t+1,s] = ((prior_sigmaEpsilon[s] + prior_sigmaPsi[s]) * prior_exp_muu[t,s] + (prior_exp_var[t,s] * u[t,s])) / 
                                       (prior_sigmaEpsilon[s] + prior_sigmaPsi[s] + prior_exp_var[t,s]) ;
                                       
      prior_exp_var[t+1,s] = ((prior_sigmaEpsilon[s] + prior_sigmaPsi[s]) * prior_exp_var[t,s] /
                              (prior_sigmaEpsilon[s] + prior_sigmaPsi[s] + prior_exp_var[t,s])) + prior_sigmaEta[s];
      
      if (prior_exp_muu[t,s] > 0.9999){
                prior_exp_mu[t,s] = 0.9999;
            } else if (prior_exp_muu[t,s] < 0.0001) {
                prior_exp_mu[t,s] = 0.0001;
            } else if (prior_exp_muu[t,s] > 0.0001 && prior_exp_muu[t,s] < 0.9999) {
                prior_exp_mu[t,s] = prior_exp_muu[t,s];
            } else {
                prior_exp_mu[t,s] = 0.5;
            }
      
      if (prior_perceptmuu[t,s] > 0.9999){
                prior_perceptmu[t,s] = 0.9999;
            } else if (prior_perceptmuu[t,s] < 0.0001) {
                prior_perceptmu[t,s] = 0.0001;
            } else if (prior_perceptmuu[t,s] > 0.0001 && prior_perceptmuu[t,s] < 0.9999) {
                prior_perceptmu[t,s] = prior_perceptmuu[t,s];
            } else {
                prior_perceptmu[t,s] = 0.5;
            }
      
      prior_percept[t,s] = beta_proportion_rng(prior_perceptmu[t,s], 1/prior_perceptvar[t,s]);
      
      prior_percept_bin[t,s] = bernoulli_rng((prior_perceptmu[t,s]^(1/prior_perceptvar[t,s]))/((prior_perceptmu[t,s]^(1/prior_perceptvar[t,s]))+(1-prior_perceptmu[t,s])^(1/prior_perceptvar[t,s])));

      prior_expectPain[t,s] = bernoulli_rng((prior_exp_mu[t,s]^(1/prior_exp_var[t,s]))/((prior_exp_mu[t,s]^(1/prior_exp_var[t,s]))+(1-prior_exp_mu[t,s])^(1/prior_exp_var[t,s])));
      
      post_percept[t,s] = beta_proportion_rng(perceptmu[t,s], 1/perceptvar[t,s]);
      
      post_percept_bin[t,s] = bernoulli_rng((perceptmu[t,s]^(1/perceptvar[t,s]))/((perceptmu[t,s]^(1/perceptvar[t,s]))+(1-perceptmu[t,s])^(1/perceptvar[t,s])));

      post_expectPain[t,s] = bernoulli_rng((exp_mu[t,s]^(1/exp_var[t,s]))/((exp_mu[t,s]^(1/exp_var[t,s]))+(1-exp_mu[t,s])^(1/exp_var[t,s])));
      
      
      
      log_lik[t,s] = bernoulli_lpmf(percept_bin[t,s] | (perceptmu[t,s]^(1/perceptvar[t,s]))/((perceptmu[t,s]^(1/perceptvar[t,s]))+(1-perceptmu[t,s])^(1/perceptvar[t,s])))+
                     beta_proportion_lpdf(percept[t,s] | perceptmu[t,s], 1/perceptvar[t,s])+
                     bernoulli_lpmf(expectPain[t,s] |  (exp_mu[t,s]^(1/exp_var[t,s]))/((exp_mu[t,s]^(1/exp_var[t,s]))+(1-exp_mu[t,s])^(1/exp_var[t,s])));
      
    }
  }
  
  
  
  
  
}
