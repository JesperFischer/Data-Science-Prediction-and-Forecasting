

data {
  int<lower=1> nsubs;
  int<lower=0> ntrials;
  matrix[ntrials, nsubs] percept;
  int pred[ntrials, nsubs];
  int percept_bin[ntrials, nsubs];
  
  matrix[ntrials, nsubs] stim;
  matrix[ntrials, nsubs] cue;
  

}

parameters {
  
  vector<lower=0,upper=1>[nsubs] alpha;
  vector<lower=0>[nsubs] percept_precision;
  vector<lower=0>[nsubs] beta;
  vector<lower=0,upper=1>[nsubs] w1;
  vector<lower=0,upper=1>[nsubs] w2;

  // Group-level parameters
  real <lower=0> kappa_alpha;
  real <lower=0, upper  = 1> mu_alpha;
  real <lower=0,upper=1> mu_w1;
  real <lower=0> kappa_w1;
  real <lower=0,upper=1> mu_w2;
  real <lower=0> kappa_w2;
 // Group-level parameters
  real <lower=0> sd_beta;
  real <lower=0> sd_percept_precision;
}


transformed parameters{
  
  matrix<lower=0,upper=1>[ntrials, nsubs] perceptmu; 
  matrix<lower=0,upper=1>[ntrials+1, nsubs] expect;
  matrix<lower=0,upper=1>[ntrials+1, nsubs] association;
  matrix<lower=0,upper=1>[ntrials+1, nsubs] ppred;
  matrix<lower=0,upper=1>[ntrials, nsubs] p_percept_bin;
  

  matrix[ntrials, nsubs] pe;
  
  for (s in 1:nsubs){
    
    association[1,s] = 0.5;
    expect[161,s] = 0.5;
   
    ppred[161,s] = (expect[161,s]^beta[s])/((expect[161,s]^beta[s])+(1-expect[161,s])^(beta[s]));
      
  
    for (t in 1:ntrials){
      
      if(cue[t,s] == 1){
          expect[t,s] = association[t,s];
        }else{
          expect[t,s] = 1-association[t,s];
        }
      
     
      perceptmu[t,s] = inv_logit(w1[s] * logit(stim[t,s]) + w2[s] * logit(expect[t,s]));
      
      
      if(cue[t,s] == 1){
          pe[t,s] = (perceptmu[t,s] - expect[t,s]);
        }else{
          pe[t,s] = -(perceptmu[t,s] - expect[t,s]);
        }
        
      association[t+1,s] = association[t,s] + alpha[s] * pe[t,s];
      
      ppred[t,s] = (expect[t,s]^beta[s])/((expect[t,s]^beta[s])+(1-expect[t,s])^(beta[s]));

      p_percept_bin[t,s] = (perceptmu[t,s]^beta[s])/((perceptmu[t,s]^beta[s])+(1-perceptmu[t,s])^(beta[s]));



    }
  }
}


model {
  for (s in 1: nsubs){
    // generating data
    for (t in 1:ntrials){
      target += beta_proportion_lpdf(percept[t,s] | perceptmu[t,s], percept_precision[s]);
      target += bernoulli_lpmf(percept_bin[t,s] | p_percept_bin[t,s]);
      target += bernoulli_lpmf(pred[t,s] |  ppred[t,s]);
    }
    // generating subject-level parameters
    target += beta_proportion_lpdf(alpha[s] | mu_alpha , kappa_alpha);
    target += beta_proportion_lpdf(w1[s] | mu_w1 , kappa_w1);
    target += beta_proportion_lpdf(w2[s] | mu_w2 , kappa_w2);

    target += lognormal_lpdf(percept_precision[s] | log(10), sd_percept_precision);
    target += lognormal_lpdf(beta[s] | log(10), sd_beta);
    
  }
  
  
  
  // Hierarchical Priors
  target += beta_proportion_lpdf(mu_alpha | 0.1 , 10) ; 
  target += lognormal_lpdf(kappa_alpha | log(30),0.5); 
  
  target += beta_proportion_lpdf(mu_w1 | 0.1 , 10) ; 
  target += lognormal_lpdf(kappa_w1 | log(30),0.5);
  
  target += beta_proportion_lpdf(mu_w2 | 0.1 , 10) ; 
  target += lognormal_lpdf(kappa_w2 | log(30),0.5);
  
  target += exponential_lpdf(sd_percept_precision | 0.5);
  target += exponential_lpdf(sd_beta | 0.5);
  
  
  
  
  
  
  
}

generated quantities{
  real prior_sd_percept_precision;
  real prior_sd_beta;
  real prior_kappa_alpha;
  real prior_mu_alpha;
  
  real prior_mu_w1;
  real prior_kappa_w1;
  
  real prior_mu_w2;
  real prior_kappa_w2;
  
  //subject level
  
  real prior_alpha[nsubs];
  real prior_percept_precision[nsubs];
  real prior_beta[nsubs];
  real prior_w1[nsubs];
  real prior_w2[nsubs];
  
  //trial level:
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_perceptmu; 
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] prior_association; 
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] prior_expect;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_p_percept_bin;
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] prior_ppred;
  
  
  matrix[ntrials, nsubs] prior_pe;
  
  matrix[ntrials, nsubs] prior_percept;
  matrix[ntrials, nsubs] prior_percept_bin;  
  matrix[ntrials, nsubs] prior_expectPain;
  
  matrix[ntrials, nsubs] post_percept;
  matrix[ntrials, nsubs] post_percept_bin;  
  matrix[ntrials, nsubs] post_expectPain;
  
  
  matrix[ntrials, nsubs] log_lik;

  
  prior_mu_w1 = beta_proportion_rng(0.1,10);
  prior_kappa_w1 = lognormal_rng(log(30),0.5);
  
  prior_mu_w2 = beta_proportion_rng(0.1,10);
  prior_kappa_w2 = lognormal_rng(log(30),0.5);
  
  prior_mu_alpha = beta_proportion_rng(0.1,10);
  prior_kappa_alpha = lognormal_rng(log(30),0.5);
  
  prior_sd_percept_precision = exponential_rng(0.5);
  prior_sd_beta = exponential_rng(0.5);
  
  
  
  for (s in 1:nsubs){
    prior_alpha[s] = beta_proportion_rng(prior_mu_alpha , prior_kappa_alpha);
    prior_w1[s] = beta_proportion_rng(prior_mu_w1 , prior_kappa_w1);
    prior_w2[s] = beta_proportion_rng(prior_mu_w2 , prior_kappa_w2);
    
    prior_percept_precision[s] = lognormal_rng(0, prior_sd_percept_precision);
    prior_beta[s] = lognormal_rng(0, prior_sd_beta);
    

    prior_association[1, s] = 0.5;
    prior_expect[161, s] = 0.5;
    

    prior_ppred[161,s] = (prior_expect[161,s]^prior_beta[s])/((prior_expect[161,s]^prior_beta[s])+(1-prior_expect[161,s])^(prior_beta[s]));
    
    for (t in 1:ntrials){
      
      
      if(cue[t,s] == 1){
        prior_expect[t,s] = prior_association[t,s];
      }else{
        prior_expect[t,s] = 1-prior_association[t,s];
      }
      
      
      prior_perceptmu[t,s] = inv_logit(prior_w1[s] * logit(stim[t,s]) + prior_w2[s] * logit(prior_expect[t,s]));
      
      
      if(cue[t,s] == 1){
        prior_pe[t,s] = (prior_perceptmu[t,s] - prior_expect[t,s]);
      }else{
        prior_pe[t,s] = -(prior_perceptmu[t,s] - prior_expect[t,s]);
      }
  
      prior_association[t+1,s] = prior_association[t,s] + prior_alpha[s] * prior_pe[t,s];
    
      prior_ppred[t,s] = (prior_expect[t,s]^prior_beta[s])/((prior_expect[t,s]^prior_beta[s])+(1-prior_expect[t,s])^(prior_beta[s]));
      

      prior_p_percept_bin[t,s] = (prior_perceptmu[t,s]^prior_beta[s])/((prior_perceptmu[t,s]^prior_beta[s])+(1-prior_perceptmu[t,s])^(prior_beta[s]));

    
    
      prior_percept[t,s] = beta_proportion_rng(prior_perceptmu[t,s], prior_percept_precision[s]);
      
      prior_percept_bin[t,s] = bernoulli_rng(prior_p_percept_bin[t,s]);

      prior_expectPain[t,s] = bernoulli_rng(prior_ppred[t,s]);
      
      
      post_percept[t,s] = beta_proportion_rng(perceptmu[t,s], percept_precision[s]);
      
      post_percept_bin[t,s] = bernoulli_rng(p_percept_bin[t,s]);

      post_expectPain[t,s] = bernoulli_rng(ppred[t,s]);
      
      
      log_lik[t,s] = bernoulli_lpmf(percept_bin[t,s] | (perceptmu[t,s]^beta[s])/((perceptmu[t,s]^beta[s])+(1-perceptmu[t,s])^(beta[s])))+
                     beta_proportion_lpdf(percept[t,s] | perceptmu[t,s], percept_precision[s])+
                     bernoulli_lpmf(pred[t,s] |  (expect[t,s]^beta[s])/((expect[t,s]^beta[s])+(1-expect[t,s])^(beta[s])));
    
    }
  } 
}
