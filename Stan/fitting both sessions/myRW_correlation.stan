data {

  int<lower=1> nsubs; // number of subjects
  int<lower=1> ntrials; // number of subjects
  int<lower=1> sesions;
  array[ntrials, nsubs,sesions] real <lower=0, upper=1> percept; // observations
  array[ntrials, nsubs,sesions] int expectPain; // prediciton
  array[ntrials, nsubs,sesions] int percept_bin; // prediciton
  
  array[ntrials, nsubs,sesions] int stim; // observations
  array[ntrials, nsubs,sesions] int cues; // observations
}

parameters {
  array[nsubs,sesions] real <lower=0, upper = 1> alpha;
  array[nsubs,sesions] real <lower=0> precision_percept;
  array[nsubs,sesions] real <lower=0> beta;
  
  array[nsubs,sesions] real <lower=0, upper = 1> w1;

  // Group-level parameters
  array[sesions] real <lower=0> kappa_alpha;
  array[sesions] real <lower=0, upper = 1> mu_alpha;
  array[sesions] real <lower=0, upper = 1> mu_w1;
  array[sesions] real <lower=0> kappa_w1;
 // Group-level parameters
  array[sesions] real <lower=0> sd_beta;
  array[sesions] real <lower=0> sd_precision_percept;
}

transformed parameters {
  array[ntrials, nsubs,sesions] real <lower=0, upper = 1> painMu; 
  array[ntrials+1, nsubs,sesions] real <lower=0, upper = 1> association; 
  array[ntrials+1, nsubs,sesions] real <lower=0, upper = 1> expectMu;
  array[ntrials, nsubs,sesions] real predErr;


  for (sess in 1:sesions){
    for (s in 1:nsubs) {
      association[1, s,sess] = 0.5;
      expectMu[161, s,sess] = 0.5;
        
      for (t in 1:ntrials){
        
        
        if(cues[t,s,sess] == 1){
          expectMu[t,s,sess] = association[t,s,sess];
        }else{
          expectMu[t,s,sess] = 1-association[t,s,sess];
        }
  
       if(w1[s,sess]*stim[t,s,sess]+(1-w1[s,sess])*expectMu[t,s,sess] == 0){
          painMu[t,s,sess] = 0.01;
        }else if (w1[s,sess]*stim[t,s,sess]+(1-w1[s,sess])*expectMu[t,s,sess] == 1){
          painMu[t,s,sess] = 0.99;
        }else{
          painMu[t,s,sess] = w1[s,sess]*stim[t,s,sess]+(1-w1[s,sess])*expectMu[t,s,sess];
          }
      
        
        if(cues[t,s,sess] == 1){
          predErr[t,s,sess] = (painMu[t,s,sess] - expectMu[t,s,sess]);
        }else{
          predErr[t,s,sess] = -(painMu[t,s,sess] - expectMu[t,s,sess]);
        }
        
        
        association[t+1,s,sess] = association[t,s,sess] + alpha[s,sess] * predErr[t,s,sess];
      
      }
    }
  }
}


model {
  
  for (sess in 1:sesions){
    for (s in 1:nsubs){
  
      for (t in 1:ntrials){
        target += beta_proportion_lpdf(percept[t,s,sess] | painMu[t,s,sess], precision_percept[s,sess]);
        
        target += bernoulli_lpmf(percept_bin[t,s,sess] | (painMu[t,s,sess]^beta[s,sess])/((painMu[t,s,sess]^beta[s,sess])+(1-painMu[t,s,sess])^(beta[s,sess])));
  
        target += bernoulli_lpmf(expectPain[t,s,sess] |  (expectMu[t,s,sess]^beta[s,sess])/((expectMu[t,s,sess]^beta[s,sess])+(1-expectMu[t,s,sess])^(beta[s,sess])));
        
      }
      
      target += beta_proportion_lpdf(alpha[s,sess] | mu_alpha[sess] , kappa_alpha[sess]);
      target += beta_proportion_lpdf(w1[s,sess] | mu_w1[sess] , kappa_w1[sess]);
  
      target += lognormal_lpdf(precision_percept[s,sess] | log(10), sd_precision_percept[sess]);
      target += lognormal_lpdf(beta[s,sess] | log(10), sd_beta[sess]);
    
  }
    // Hierarchical Priors
    target += beta_proportion_lpdf(mu_alpha[sess] | 0.3 , 3) ; 
    target += lognormal_lpdf(kappa_alpha[sess] | log(30) , 0.5); 
    
    target += beta_proportion_lpdf(mu_w1[sess] | 0.3 , 3) ; 
    target += lognormal_lpdf(kappa_w1[sess] | log(30) , 0.5);
    
    target += exponential_lpdf(sd_precision_percept[sess] | 1);
    target += exponential_lpdf(sd_beta[sess] | 1);
  }
}




generated quantities{
  
      array[sesions] real <lower=0> prior_sd_precision_percept;
      array[sesions] real <lower=0> prior_sd_beta;
      
      array[sesions] real <lower=0> prior_kappa_alpha;
      array[sesions] real <lower=0, upper=1> prior_mu_alpha;
      
      array[sesions] real <lower=0, upper=1> prior_mu_w1;
      array[sesions] real <lower=0> prior_kappa_w1;
      
      //subject level
      
      array[nsubs,sesions] real <lower=0, upper = 1> prior_alpha;
      array[nsubs,sesions] real <lower=0> prior_precision_percept;
      array[nsubs,sesions] real <lower=0> prior_beta;
      array[nsubs,sesions] real <lower=0, upper = 1> prior_w1;
      
      //trial level:
      array[ntrials, nsubs,sesions] real <lower=0, upper  = 1> prior_painMu; 
      array[ntrials+1, nsubs,sesions] real <lower=0, upper  = 1> prior_association; 
      array[ntrials+1, nsubs,sesions] real <lower=0, upper  = 1> prior_expectMu;
      array[ntrials, nsubs,sesions] real prior_predErr;
      
      array[ntrials, nsubs,sesions] real <lower=0, upper  = 1> prior_percept;
      array[ntrials, nsubs,sesions] real <lower=0, upper  = 1> prior_percept_bin;  
      array[ntrials, nsubs,sesions] real <lower=0, upper  = 1> prior_expectPain;

      array[ntrials, nsubs,sesions] real <lower=0, upper  = 1> post_percept;
      array[ntrials, nsubs,sesions] real <lower=0, upper  = 1> post_percept_bin;  
      array[ntrials, nsubs,sesions] real <lower=0, upper  = 1> post_expectPain;
      
      array[ntrials, nsubs,sesions] real log_lik;
      
      
      

      for(sess in 1:sesions){
        prior_mu_w1[sess] = beta_proportion_rng(0.3,3);
        prior_kappa_w1[sess] = lognormal_rng(log(30) , 0.5);
        
        prior_mu_alpha[sess] = beta_proportion_rng(0.3,3);
        prior_kappa_alpha[sess] = lognormal_rng(log(30) , 0.5);
        
        prior_sd_precision_percept[sess] = exponential_rng(1);
        prior_sd_beta[sess] = exponential_rng(1);
        
        
        
        for (s in 1:nsubs){
          prior_alpha[s,sess] = beta_proportion_rng(prior_mu_alpha[sess] , prior_kappa_alpha[sess]);
          prior_w1[s,sess] = beta_proportion_rng(prior_mu_w1[sess] , prior_kappa_w1[sess]);
          
          prior_precision_percept[s,sess] = lognormal_rng(log(20), prior_sd_precision_percept[sess]);
          prior_beta[s,sess] = lognormal_rng(log(20), prior_sd_beta[sess]);
          
      
          prior_association[1, s,sess] = 0.5;
          prior_expectMu[161, s,sess] = 0.5;
            
          for (t in 1:ntrials){
    
            if(cues[t,s,sess] == 1){
              prior_expectMu[t,s,sess] = prior_association[t,s,sess];
            }else{
              prior_expectMu[t,s,sess] = 1-prior_association[t,s,sess];
            }
      
            
            
             if(prior_w1[s,sess]*stim[t,s,sess]+(1-prior_w1[s,sess])*prior_expectMu[t,s,sess] == 0){
                prior_painMu[t,s,sess] = 0.01;
              }else if (prior_w1[s,sess]*stim[t,s,sess]+(1-prior_w1[s,sess])*prior_expectMu[t,s,sess] == 1){
                prior_painMu[t,s,sess] = 0.99;
              }else{
                prior_painMu[t,s,sess] = prior_w1[s,sess]*stim[t,s,sess]+(1-prior_w1[s,sess])*prior_expectMu[t,s,sess];
                }
            
            
            if(cues[t,s,sess] == 1){
              prior_predErr[t,s,sess] = (prior_painMu[t,s,sess] - prior_expectMu[t,s,sess]);
            }else{
              prior_predErr[t,s,sess] = -(prior_painMu[t,s,sess] - prior_expectMu[t,s,sess]);
            }
        
            prior_association[t+1,s,sess] = prior_association[t,s,sess] + prior_alpha[s,sess] * prior_predErr[t,s,sess];
        
            prior_percept[t,s,sess] = beta_proportion_rng(prior_painMu[t,s,sess], prior_precision_percept[s,sess]);
            
            prior_percept_bin[t,s,sess] = bernoulli_rng((prior_painMu[t,s,sess]^prior_beta[s,sess])/((prior_painMu[t,s,sess]^prior_beta[s,sess])+(1-prior_painMu[t,s,sess])^(prior_beta[s,sess])));
      
            prior_expectPain[t,s,sess] = bernoulli_rng((prior_expectMu[t,s,sess]^prior_beta[s,sess])/((prior_expectMu[t,s,sess]^prior_beta[s,sess])+(1-prior_expectMu[t,s,sess])^(prior_beta[s,sess])));
          
            
            post_percept[t,s,sess] = beta_proportion_rng(painMu[t,s,sess], precision_percept[s,sess]);
            
            post_percept_bin[t,s,sess] = bernoulli_rng((painMu[t,s,sess]^beta[s,sess])/((painMu[t,s,sess]^beta[s,sess])+(1-painMu[t,s,sess])^(beta[s,sess])));
      
            post_expectPain[t,s,sess] = bernoulli_rng((expectMu[t,s,sess]^beta[s,sess])/((expectMu[t,s,sess]^beta[s,sess])+(1-expectMu[t,s,sess])^(beta[s,sess])));
            
            
            log_lik[t,s,sess] = bernoulli_lpmf(percept_bin[t,s,sess] | (painMu[t,s,sess]^beta[s,sess])/((painMu[t,s,sess]^beta[s,sess])+(1-painMu[t,s,sess])^(beta[s,sess])))+
                           beta_proportion_lpdf(percept[t,s,sess] | painMu[t,s,sess], precision_percept[s,sess])+
                           bernoulli_lpmf(expectPain[t,s,sess] |  (expectMu[t,s,sess]^beta[s,sess])/((expectMu[t,s,sess]^beta[s,sess])+(1-expectMu[t,s,sess])^(beta[s,sess])));
          
          }
        }   
    }
  }

