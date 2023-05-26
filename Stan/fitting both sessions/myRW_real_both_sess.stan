data {

  int<lower=1> nsubs; // number of subjects
  int<lower=1> ntrials; // number of subjects
  int<lower=1> sesions; // number of subjects
  
  array [ntrials, nsubs,sesions] real percept; // observations
  array [ntrials, nsubs,sesions] int expectPain; // prediciton
  array [ntrials, nsubs,sesions] int percept_bin; // prediciton
  
  array [ntrials, nsubs,sesions] int stim; // observations
  array [ntrials, nsubs,sesions] int cues; // observations
}

parameters {
  array[nsubs,sesions] real <lower=0, upper = 1> alpha;
  array[nsubs,sesions] real <lower=0> precision_percept;
  
  array[nsubs,sesions] real <lower=0, upper = 1> w1;
  // Group-level parameters
  array[sesions] real <lower=0> kappa_alpha;
  array[sesions] real <lower=0, upper = 1> mu_alpha;
  array[sesions] real <lower=0, upper = 1> mu_w1;
  array[sesions] real <lower=0> kappa_w1;
 // Group-level parameters
 
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
        
        target += bernoulli_lpmf(percept_bin[t,s,sess] | painMu[t,s,sess]);
  
        target += bernoulli_lpmf(expectPain[t,s,sess] |  expectMu[t,s,sess]);
        
      }
      
      target += beta_proportion_lpdf(alpha[s,sess] | mu_alpha[sess] , kappa_alpha[sess]);
      target += beta_proportion_lpdf(w1[s,sess] | mu_w1[sess] , kappa_w1[sess]);
  
      target += lognormal_lpdf(precision_percept[s,sess] | log(10), sd_precision_percept[sess]);
    
  }
    // Hierarchical Priors
    target += beta_proportion_lpdf(mu_alpha[sess] | 0.3 , 3) ; 
    target += lognormal_lpdf(kappa_alpha[sess] | log(30) , 0.5); 
    
    target += beta_proportion_lpdf(mu_w1[sess] | 0.3 , 3) ; 
    target += lognormal_lpdf(kappa_w1[sess] | log(30) , 0.5);
    
    target += exponential_lpdf(sd_precision_percept[sess] | 0.3);
  }
}


generated quantities{
    
      array[sesions] real <lower=0> prior_sd_precision_percept;
      
      array[sesions] real <lower=0> prior_kappa_alpha;
      array[sesions] real <lower=0, upper=1> prior_mu_alpha;
      
      array[sesions] real <lower=0, upper=1> prior_mu_w1;
      array[sesions] real <lower=0> prior_kappa_w1;
      
      //subject level
      
      array[nsubs,sesions] real <lower=0, upper = 1> prior_alpha;
      array[nsubs,sesions] real <lower=0> prior_precision_percept;
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
          
      real correlation_alpha;
      real sum_squared_diff_alpha;
      real sum_product_diff_alpha;

      array[nsubs] real diff_s1_alpha;
      array[nsubs] real diff_s2_alpha;
      array[nsubs] real sq_diff_s1_alpha;
      array[nsubs] real sq_diff_s2_alpha;

      array[nsubs] real product_diff_alpha;
      
      
      real correlation_w1;
      real sum_squared_diff_w1;
      real sum_product_diff_w1;

      array[nsubs] real diff_s1_w1;
      array[nsubs] real diff_s2_w1;
      array[nsubs] real sq_diff_s1_w1;
      array[nsubs] real sq_diff_s2_w1;

      array[nsubs] real product_diff_w1;
      
      for (s in 1:nsubs){
        diff_s1_w1[s] = mean(w1[,1])-w1[s,1];
        diff_s2_w1[s] = mean(w1[,2])-w1[s,2];
        
        product_diff_w1[s] = diff_s1_w1[s] * diff_s2_w1[s];
        
        sq_diff_s1_w1[s] = (mean(w1[,1])-w1[s,1])^2;
        sq_diff_s2_w1[s] = (mean(w1[,2])-w1[s,2])^2;
      }
      
      sum_squared_diff_w1 = sum(sq_diff_s1_w1) * sum(sq_diff_s2_w1);
      
      sum_product_diff_w1 = sum(product_diff_w1);

      correlation_w1 = sum_product_diff_w1 / sqrt(sum_squared_diff_w1);
      
      
      
      
      for (s in 1:nsubs){
        diff_s1_alpha[s] = mean(alpha[,1])-alpha[s,1];
        diff_s2_alpha[s] = mean(alpha[,2])-alpha[s,2];
        
        product_diff_alpha[s] = diff_s1_alpha[s] * diff_s2_alpha[s];
        
        sq_diff_s1_alpha[s] = (mean(alpha[,1])-alpha[s,1])^2;
        sq_diff_s2_alpha[s] = (mean(alpha[,2])-alpha[s,2])^2;
      }
      
      sum_squared_diff_alpha = sum(sq_diff_s1_alpha) * sum(sq_diff_s2_alpha);
      
      sum_product_diff_alpha = sum(product_diff_alpha);

      correlation_alpha = sum_product_diff_alpha / sqrt(sum_squared_diff_alpha);
      

      for(sess in 1:sesions){
        prior_mu_w1[sess] = beta_proportion_rng(0.3,3);
        prior_kappa_w1[sess] = lognormal_rng(log(30) , 0.5);
        
        prior_mu_alpha[sess] = beta_proportion_rng(0.3,3);
        prior_kappa_alpha[sess] = lognormal_rng(log(30) , 0.5);
        
        prior_sd_precision_percept[sess] = exponential_rng(0.3);
        
        
        
        for (s in 1:nsubs){
          prior_alpha[s,sess] = beta_proportion_rng(prior_mu_alpha[sess] , prior_kappa_alpha[sess]);
          prior_w1[s,sess] = beta_proportion_rng(prior_mu_w1[sess] , prior_kappa_w1[sess]);
          
          prior_precision_percept[s,sess] = lognormal_rng(log(20), prior_sd_precision_percept[sess]);

      
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
            
            prior_percept_bin[t,s,sess] = bernoulli_rng(prior_painMu[t,s,sess]);
      
            prior_expectPain[t,s,sess] = bernoulli_rng(prior_expectMu[t,s,sess]);
          
            
            post_percept[t,s,sess] = beta_proportion_rng(painMu[t,s,sess], precision_percept[s,sess]);
            
            post_percept_bin[t,s,sess] = bernoulli_rng(painMu[t,s,sess]);
      
            post_expectPain[t,s,sess] = bernoulli_rng(expectMu[t,s,sess]);
            
            
            log_lik[t,s,sess] = bernoulli_lpmf(percept_bin[t,s,sess] | painMu[t,s,sess])+
                           beta_proportion_lpdf(percept[t,s,sess] | painMu[t,s,sess], precision_percept[s,sess])+
                           bernoulli_lpmf(expectPain[t,s,sess] |  expectMu[t,s,sess]);
          
          }
        }   
    }
  }

