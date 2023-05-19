data {

  int<lower=1> nsubs; // number of subjects
  int<lower=1> ntrials; // number of subjects
  matrix[ntrials, nsubs] percept; // observations
  int expectPain[ntrials, nsubs]; // prediciton
  int percept_bin[ntrials, nsubs]; // prediciton
  
  matrix[ntrials, nsubs] stim; // observations
  matrix[ntrials, nsubs] cues; // observations
}

parameters {
  vector<lower=0,upper=1>[nsubs] alpha;
  vector<lower=0>[nsubs] precision_percept;
  vector<lower=0>[nsubs] beta;
  
  vector<lower=0, upper=1>[nsubs] w1;

  // Group-level parameters
  real <lower=0> kappa_alpha;
  real <lower=0, upper  = 1> mu_alpha;
  real <lower=0, upper = 1> mu_w1;
  real <lower=0> kappa_w1;
 // Group-level parameters
  real <lower=0> sd_beta;
  real <lower=0> sd_precision_percept;
}

transformed parameters {
  matrix <lower=0, upper  = 1> [ntrials, nsubs] painMu; 
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] association; 
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] expectMu;
  matrix[ntrials, nsubs] predErr;


  for (s in 1:nsubs) {
    association[1, s] = 0.5;
    expectMu[161, s] = 0.5;
      
    for (t in 1:ntrials){
      
      
      if(cues[t,s] == 1){
        expectMu[t,s] = association[t,s];
      }else{
        expectMu[t,s] = 1-association[t,s];
      }

     if(w1[s]*stim[t,s]+(1-w1[s])*expectMu[t,s] == 0){
        painMu[t,s] = 0.01;
      }else if (w1[s]*stim[t,s]+(1-w1[s])*expectMu[t,s] == 1){
        painMu[t,s] = 0.99;
      }else{
        painMu[t,s] = w1[s]*stim[t,s]+(1-w1[s])*expectMu[t,s];
        }
    
      
      if(cues[t,s] == 1){
        predErr[t,s] = (painMu[t,s] - expectMu[t,s]);
      }else{
        predErr[t,s] = -(painMu[t,s] - expectMu[t,s]);
      }
      
      
      association[t+1,s] = association[t,s] + alpha[s] * predErr[t,s];
    
    }
  }
}


model {
  for (s in 1:nsubs){

    for (t in 1:ntrials){
      target += beta_proportion_lpdf(percept[t,s] | painMu[t,s], precision_percept[s]);
      
      target += bernoulli_lpmf(percept_bin[t,s] | (painMu[t,s]^beta[s])/((painMu[t,s]^beta[s])+(1-painMu[t,s])^(beta[s])));

      target += bernoulli_lpmf(expectPain[t,s] |  (expectMu[t,s]^beta[s])/((expectMu[t,s]^beta[s])+(1-expectMu[t,s])^(beta[s])));
      
    }
    
    target += beta_proportion_lpdf(alpha[s] | mu_alpha , kappa_alpha);
    target += beta_proportion_lpdf(w1[s] | mu_w1 , kappa_w1);

    target += lognormal_lpdf(precision_percept[s] | log(10), sd_precision_percept);
    target += lognormal_lpdf(beta[s] | log(10), sd_beta);
    
  }
  
  // Hierarchical Priors
  target += beta_proportion_lpdf(mu_alpha | 0.3 , 3) ; 
  target += lognormal_lpdf(kappa_alpha | log(30) , 0.5); 
  
  target += beta_proportion_lpdf(mu_w1 | 0.3 , 3) ; 
  target += lognormal_lpdf(kappa_w1 | log(30) , 0.5);
  
  target += exponential_lpdf(sd_precision_percept | 1);
  target += exponential_lpdf(sd_beta | 1);
}

generated quantities{
  
  
      real <lower=0> prior_sd_precision_percept;
      real <lower=0> prior_sd_beta;
      
      real <lower=0> prior_kappa_alpha;
      real <lower=0, upper = 1> prior_mu_alpha;
      
      real <lower=0, upper = 1> prior_mu_w1;
      real <lower=0> prior_kappa_w1;
      
      //subject level
      
      real <lower=0, upper = 1> prior_alpha[nsubs];
      real <lower=0> prior_precision_percept[nsubs];
      real <lower=0> prior_beta[nsubs];
      real <lower=0, upper = 1> prior_w1[nsubs];
      
      //trial level:
      matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_painMu; 
      matrix <lower=0, upper  = 1> [ntrials+1, nsubs] prior_association; 
      matrix <lower=0, upper  = 1> [ntrials+1, nsubs] prior_expectMu;
      matrix[ntrials, nsubs] prior_predErr;
      
      matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_percept;
      matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_percept_bin;  
      matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_expectPain;

      matrix <lower=0, upper  = 1> [ntrials, nsubs] post_percept;
      matrix <lower=0, upper  = 1> [ntrials, nsubs] post_percept_bin;  
      matrix <lower=0, upper  = 1> [ntrials, nsubs] post_expectPain;
      
      
      matrix[ntrials, nsubs] log_lik;
      

  
      prior_mu_w1 = beta_proportion_rng(0.3,3);
      prior_kappa_w1 = lognormal_rng(log(30) , 0.5);
      
      prior_mu_alpha = beta_proportion_rng(0.3,3);
      prior_kappa_alpha = lognormal_rng(log(30) , 0.5);
      
      prior_sd_precision_percept = exponential_rng(1);
      prior_sd_beta = exponential_rng(1);
      
      
      
      for (s in 1:nsubs){
        prior_alpha[s] = beta_proportion_rng(prior_mu_alpha , prior_kappa_alpha);
        prior_w1[s] = beta_proportion_rng(prior_mu_w1 , prior_kappa_w1);
        
        prior_precision_percept[s] = lognormal_rng(log(20), prior_sd_precision_percept);
        prior_beta[s] = lognormal_rng(log(20), prior_sd_beta);
        
    
        prior_association[1, s] = 0.5;
        prior_expectMu[161, s] = 0.5;
          
        for (t in 1:ntrials){
          
          
          if(cues[t,s] == 1){
            prior_expectMu[t,s] = prior_association[t,s];
          }else{
            prior_expectMu[t,s] = 1-prior_association[t,s];
          }
    
          
          
           if(prior_w1[s]*stim[t,s]+(1-prior_w1[s])*prior_expectMu[t,s] == 0){
              prior_painMu[t,s] = 0.01;
            }else if (prior_w1[s]*stim[t,s]+(1-prior_w1[s])*prior_expectMu[t,s] == 1){
              prior_painMu[t,s] = 0.99;
            }else{
              prior_painMu[t,s] = prior_w1[s]*stim[t,s]+(1-prior_w1[s])*prior_expectMu[t,s];
              }
          
          
          if(cues[t,s] == 1){
            prior_predErr[t,s] = (prior_painMu[t,s] - prior_expectMu[t,s]);
          }else{
            prior_predErr[t,s] = -(prior_painMu[t,s] - prior_expectMu[t,s]);
          }
      
          prior_association[t+1,s] = prior_association[t,s] + prior_alpha[s] * prior_predErr[t,s];
      
          prior_percept[t,s] = beta_proportion_rng(prior_painMu[t,s], prior_precision_percept[s]);
          
          prior_percept_bin[t,s] = bernoulli_rng((prior_painMu[t,s]^prior_beta[s])/((prior_painMu[t,s]^prior_beta[s])+(1-prior_painMu[t,s])^(prior_beta[s])));
    
          prior_expectPain[t,s] = bernoulli_rng((prior_expectMu[t,s]^prior_beta[s])/((prior_expectMu[t,s]^prior_beta[s])+(1-prior_expectMu[t,s])^(prior_beta[s])));
        
          
          post_percept[t,s] = beta_proportion_rng(painMu[t,s], precision_percept[s]);
          
          post_percept_bin[t,s] = bernoulli_rng((painMu[t,s]^beta[s])/((painMu[t,s]^beta[s])+(1-painMu[t,s])^(beta[s])));
    
          post_expectPain[t,s] = bernoulli_rng((expectMu[t,s]^beta[s])/((expectMu[t,s]^beta[s])+(1-expectMu[t,s])^(beta[s])));
          
          
          log_lik[t,s] = bernoulli_lpmf(percept_bin[t,s] | (painMu[t,s]^beta[s])/((painMu[t,s]^beta[s])+(1-painMu[t,s])^(beta[s])))+
                         beta_proportion_lpdf(percept[t,s] | painMu[t,s], precision_percept[s])+
                         bernoulli_lpmf(expectPain[t,s] |  (expectMu[t,s]^beta[s])/((expectMu[t,s]^beta[s])+(1-expectMu[t,s])^(beta[s])));
        
        }
      }   
    
  }


