data {
  int<lower=1> ntrials; //number of ntrials
  int<lower=1> nsubs; //number of nsubs
  int<lower=1> sesions; //number of nsubs
  
  
  array[ntrials,nsubs,sesions] int u;
  array [ntrials, nsubs,sesions] real percept; // observations
  array [ntrials, nsubs,sesions] int pred; // prediciton
  array [ntrials, nsubs,sesions] int percept_bin; // prediciton
  
  
  array [ntrials, nsubs,sesions] int stim; // observations
  array [ntrials, nsubs,sesions] int cues; // observations
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  array[nsubs,sesions] real  b0;
  array[nsubs,sesions] real  b1;
  array[nsubs,sesions] real  b2;
  
  
  array[nsubs,sesions] real  b0_b;
  array[nsubs,sesions] real  b1_b;
  array[nsubs,sesions] real  b2_b;
  
  
  array[nsubs,sesions] real <lower = 0, upper = 1> alpha;
  
  array[nsubs,sesions] real <lower = 0> percept_precision;
  
  real <lower = 0, upper = 1> mu_alpha[sesions];
  real <lower = 0 > kappa_alpha[sesions];
  
  
  real <lower=0> sd_percept_precision[sesions];
  
  real mu_b0[sesions];
  real <lower = 0> sd_b0[sesions];
  
  real  mu_b1[sesions];
  real <lower = 0> sd_b1[sesions];
  
  real  mu_b2[sesions];
  real <lower = 0> sd_b2[sesions];
  
  
  real mu_b0_b[sesions];
  real <lower = 0> sd_b0_b[sesions];
  
  real  mu_b1_b[sesions];
  real <lower = 0> sd_b1_b[sesions];
  
  real  mu_b2_b[sesions];
  real <lower = 0> sd_b2_b[sesions];
  
}


transformed parameters{
  array[ntrials,nsubs,sesions] real <lower = 0, upper = 1> belief;
  array[ntrials,nsubs,sesions] real <lower = 0, upper = 1> expect;
  
  
  for (sess in 1:sesions){
    for (s in 1:nsubs){
      belief[1,s,sess] = 0.5;
      expect[1,s,sess] = 0.5;
      
  
      for (t in 2:ntrials){
        
        
        belief[t,s,sess] = belief[t-1,s,sess]+alpha[s,sess]*(u[t-1,s,sess]-belief[t-1,s,sess]);
        
        if(cues[t,s,sess] == 1){
          expect[t,s,sess] = belief[t,s,sess];
        }else{
          expect[t,s,sess] = 1-belief[t,s,sess];
        }
        
      
      }
    }
  }
  
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  for (sess in 1:sesions){
    target += beta_lpdf(mu_alpha[sess] | 1,1);
    target += lognormal_lpdf(kappa_alpha[sess] | log(20),0.6);
    
    target +=normal_lpdf(mu_b0[sess] | 0,2);
    target +=lognormal_lpdf(sd_b0[sess] | 0,0.5);
    
    target +=normal_lpdf(mu_b1[sess] | 0,2);
    target +=lognormal_lpdf(sd_b1[sess] | 0,0.5);
    
    target +=normal_lpdf(mu_b2[sess] | 0,2);
    target +=lognormal_lpdf(sd_b2[sess] | 0,0.5);
    
    
    
    target +=normal_lpdf(mu_b0_b[sess] | 0,2);
    target +=lognormal_lpdf(sd_b0_b[sess] | 0,0.5);
    
    target +=normal_lpdf(mu_b1_b[sess] | 0,2);
    target +=lognormal_lpdf(sd_b1_b[sess] | 0,0.5);
    
    
    target +=normal_lpdf(mu_b2_b[sess] | 0,2);
    target +=lognormal_lpdf(sd_b2_b[sess] | 0,0.5);
    
    
    
    target += exponential_lpdf(sd_percept_precision[sess] | 0.5);
  
      
    for (s in 1:nsubs){
      
      target +=beta_proportion_lpdf(alpha[s,sess] | mu_alpha[sess],kappa_alpha[sess]);
      
      target += lognormal_lpdf(percept_precision[s,sess] | log(20), sd_percept_precision[sess]);
  
      
      target +=normal_lpdf(b0[s] | mu_b0[sess],sd_b0[sess]);
      target +=normal_lpdf(b1[s] | mu_b1[sess],sd_b1[sess]);
      target +=normal_lpdf(b2[s] | mu_b2[sess],sd_b2[sess]);
      
      target +=normal_lpdf(b0_b[s,sess] | mu_b0_b[sess],sd_b0_b[sess]);
      target +=normal_lpdf(b1_b[s,sess] | mu_b1_b[sess],sd_b1_b[sess]);
      target +=normal_lpdf(b2_b[s,sess] | mu_b2_b[sess],sd_b2_b[sess]);    
      
      for (t in 1:ntrials){
      
        target +=bernoulli_lpmf(pred[t,s,sess] | belief[t,s,sess]);
        target +=beta_proportion_lpdf(percept[t,s,sess] | inv_logit(b0[s,sess]+b1[s,sess]*stim[t,s,sess]+b2[s,sess]*expect[t,s,sess]), percept_precision[s,sess]);
        target += bernoulli_lpmf(percept_bin[t,s,sess] |  inv_logit(b0_b[s,sess]+b1_b[s,sess]*stim[t,s,sess]+b2_b[s,sess]*expect[t,s,sess]));
  
      }
    }
  }

}



generated quantities{
  array[ntrials, nsubs,sesions] real log_lik;
      
  
  real correlation_alpha;
  real sum_squared_diff_alpha;
  real sum_product_diff_alpha;

  array[nsubs] real diff_s1_alpha;
  array[nsubs] real diff_s2_alpha;
  array[nsubs] real sq_diff_s1_alpha;
  array[nsubs] real sq_diff_s2_alpha;

  array[nsubs] real product_diff_alpha;
  
  
  real correlation_b2;
  real sum_squared_diff_b2;
  real sum_product_diff_b2;

  array[nsubs] real diff_s1_b2;
  array[nsubs] real diff_s2_b2;
  array[nsubs] real sq_diff_s1_b2;
  array[nsubs] real sq_diff_s2_b2;

  array[nsubs] real product_diff_b2;
  
  
  real correlation_b1;
  real sum_squared_diff_b1;
  real sum_product_diff_b1;

  array[nsubs] real diff_s1_b1;
  array[nsubs] real diff_s2_b1;
  array[nsubs] real sq_diff_s1_b1;
  array[nsubs] real sq_diff_s2_b1;

  array[nsubs] real product_diff_b1;
  
  
  real correlation_b1_b;
  real sum_squared_diff_b1_b;
  real sum_product_diff_b1_b;

  array[nsubs] real diff_s1_b1_b;
  array[nsubs] real diff_s2_b1_b;
  array[nsubs] real sq_diff_s1_b1_b;
  array[nsubs] real sq_diff_s2_b1_b;

  array[nsubs] real product_diff_b1_b;
  
  
    real correlation_b2_b;
  real sum_squared_diff_b2_b;
  real sum_product_diff_b2_b;

  array[nsubs] real diff_s1_b2_b;
  array[nsubs] real diff_s2_b2_b;
  array[nsubs] real sq_diff_s1_b2_b;
  array[nsubs] real sq_diff_s2_b2_b;

  array[nsubs] real product_diff_b2_b;
  
  for (s in 1:nsubs){
    diff_s1_b2_b[s] = mean(b2_b[,1])-b2_b[s,1];
    diff_s2_b2_b[s] = mean(b2_b[,2])-b2_b[s,2];
    
    product_diff_b2_b[s] = diff_s1_b2_b[s] * diff_s2_b2_b[s];
    
    sq_diff_s1_b2_b[s] = (mean(b2_b[,1])-b2_b[s,1])^2;
    sq_diff_s2_b2_b[s] = (mean(b2_b[,2])-b2_b[s,2])^2;
  }
  
  sum_squared_diff_b2_b = sum(sq_diff_s1_b2_b) * sum(sq_diff_s2_b2_b);
  
  sum_product_diff_b2_b = sum(product_diff_b2_b);

  correlation_b2_b = sum_product_diff_b2_b / sqrt(sum_squared_diff_b2_b);
  
  for (s in 1:nsubs){
    diff_s1_b1_b[s] = mean(b1_b[,1])-b1_b[s,1];
    diff_s2_b1_b[s] = mean(b1_b[,2])-b1_b[s,2];
    
    product_diff_b1_b[s] = diff_s1_b1_b[s] * diff_s2_b1_b[s];
    
    sq_diff_s1_b1_b[s] = (mean(b1_b[,1])-b1_b[s,1])^2;
    sq_diff_s2_b1_b[s] = (mean(b1_b[,2])-b1_b[s,2])^2;
  }
  
  sum_squared_diff_b1_b = sum(sq_diff_s1_b1_b) * sum(sq_diff_s2_b1_b);
  
  sum_product_diff_b1_b = sum(product_diff_b1_b);

  correlation_b1_b = sum_product_diff_b1_b / sqrt(sum_squared_diff_b1_b);
  
  
  for (s in 1:nsubs){
    diff_s1_b1[s] = mean(b1[,1])-b1[s,1];
    diff_s2_b1[s] = mean(b1[,2])-b1[s,2];
    
    product_diff_b1[s] = diff_s1_b1[s] * diff_s2_b1[s];
    
    sq_diff_s1_b1[s] = (mean(b1[,1])-b1[s,1])^2;
    sq_diff_s2_b1[s] = (mean(b1[,2])-b1[s,2])^2;
  }
  
  sum_squared_diff_b1 = sum(sq_diff_s1_b1) * sum(sq_diff_s2_b1);
  
  sum_product_diff_b1 = sum(product_diff_b1);

  correlation_b1 = sum_product_diff_b1 / sqrt(sum_squared_diff_b1);
  
  
  
  
  for (s in 1:nsubs){
    diff_s1_b2[s] = mean(b2[,1])-b2[s,1];
    diff_s2_b2[s] = mean(b2[,2])-b2[s,2];
    
    product_diff_b2[s] = diff_s1_b2[s] * diff_s2_b2[s];
    
    sq_diff_s1_b2[s] = (mean(b2[,1])-b2[s,1])^2;
    sq_diff_s2_b2[s] = (mean(b2[,2])-b2[s,2])^2;
  }
  
  sum_squared_diff_b2 = sum(sq_diff_s1_b2) * sum(sq_diff_s2_b2);
  
  sum_product_diff_b2 = sum(product_diff_b2);

  correlation_b2 = sum_product_diff_b2 / sqrt(sum_squared_diff_b2);
  
      
      
        
        
        
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
  

  
  
  // array[nsubs] real  prior_b0;
  // array[nsubs] real  prior_b1;
  // array[nsubs] real  prior_b2;
  // 
  // array[nsubs] real  prior_b0_b;
  // array[nsubs] real  prior_b1_b;
  // array[nsubs] real  prior_b2_b;
  // 
  // array[nsubs] real <lower = 0, upper = 1> prior_alpha;
  // 
  // vector<lower=0>[nsubs] prior_precision_percept;
  // 
  // real <lower = 0, upper = 1> prior_mu_alpha;
  // real <lower = 0, upper = 1> prior_sd_alpha;
  // 
  // 
  // 
  // real <lower=0> prior_sd_precision_percept;
  // 
  // real prior_mu_b0;
  // real <lower = 0> prior_sd_b0;
  // 
  // real  prior_mu_b1;
  // real <lower = 0> prior_sd_b1;
  // 
  // real  prior_mu_b2;
  // real <lower = 0> prior_sd_b2;
  // 
  // real prior_mu_b0_b;
  // real <lower = 0> prior_sd_b0_b;
  // 
  // real  prior_mu_b1_b;
  // real <lower = 0> prior_sd_b1_b;
  // 
  // real  prior_mu_b2_b;
  // real <lower = 0> prior_sd_b2_b;
  // 
  // array[ntrials,nsubs] real <lower = 0, upper = 1> prior_belief;
  // array[ntrials,nsubs] real <lower = 0, upper = 1> prior_expect;
  // 
  // array[ntrials,nsubs] real <lower = 0, upper = 1> prior_percept;
  // array[ntrials,nsubs] int <lower = 0, upper = 1> prior_percept_bin;
  // array[ntrials,nsubs] int <lower = 0, upper = 1> prior_pred;
  // 
  // array[ntrials,nsubs] real <lower = 0, upper = 1> post_percept;
  // array[ntrials,nsubs] int <lower = 0, upper = 1> post_percept_bin;
  // array[ntrials,nsubs] int <lower = 0, upper = 1> post_pred;
  // 
  // 
  // 
  for (sess in 1:sesions){
  // 
  // prior_mu_alpha = beta_rng(1,1);
  // prior_sd_alpha = beta_rng(1,1);
  // 
  // prior_mu_b0 =normal_rng(0,1);
  // prior_sd_b0 =lognormal_rng(0,0.5);
  // 
  // prior_mu_b1 =normal_rng(0,1);
  // prior_sd_b1 =lognormal_rng(0,0.5);
  // 
  // prior_mu_b2 =normal_rng(0,1);
  // prior_sd_b2 =lognormal_rng(0,0.5);
  // 
  // 
  // prior_mu_b0_b =normal_rng(0,1);
  // prior_sd_b0_b =lognormal_rng(0,0.5);
  // 
  // prior_mu_b1_b =normal_rng(0,1);
  // prior_sd_b1_b =lognormal_rng(0,0.5);
  // 
  // 
  // prior_mu_b2_b =normal_rng(0,1);
  // prior_sd_b2_b =lognormal_rng(0,0.5);
  // 
  // prior_sd_precision_percept = exponential_rng(0.1);
  
  
    for (s in 1:nsubs){
    
    // prior_precision_percept[s] = lognormal_rng(log(20), prior_sd_precision_percept);
    // 


    // prior_b0[s] = normal_rng(prior_mu_b0,prior_sd_b0);
    // prior_b1[s] = normal_rng(prior_mu_b1,prior_sd_b1);
    // prior_b2[s] = normal_rng(prior_mu_b2,prior_sd_b2);
    // 
    // prior_b0_b[s] = normal_rng(prior_mu_b0_b,prior_sd_b0_b);
    // prior_b1_b[s] = normal_rng(prior_mu_b1_b,prior_sd_b1_b);
    // prior_b2_b[s] = normal_rng(prior_mu_b2_b,prior_sd_b2_b);    
    // 
    // 
    // prior_belief[1,s] = 0.5;
    // prior_expect[1,s] = 0.5;

      for (t in 1:ntrials){
      
   //    
   //    prior_belief[t,s] = prior_belief[t-1,s]+prior_alpha[s]*(u[t-1,s]-prior_belief[t-1,s]);
   //    
   //    if(cues[t,s] == 1){
   //      prior_expect[t,s] = prior_belief[t,s];
   //    }else{
   //      prior_expect[t,s] = 1-prior_belief[t,s];
   //    }
   //  
   // 
   //  
   //  
   // prior_percept[t,s] = beta_proportion_rng(inv_logit(prior_b0[s]+prior_b1[s]*stim[t,s]+prior_b2[s]*prior_expect[t,s]), prior_precision_percept[s]);
   //        
   // prior_percept_bin[t,s] = bernoulli_rng(inv_logit(prior_b0_b[s]+prior_b1_b[s]*stim[t,s]+prior_b1_b[s]*prior_expect[t,s]));
   //  
   // prior_pred[t,s] = bernoulli_rng(prior_belief[t,s]);
   //  
   //  
   // post_percept[t,s] = beta_proportion_rng(inv_logit(b0[s]+b1[s]*stim[t,s]+b2[s]*expect[t,s]), precision_percept[s]);
   //        
   // post_percept_bin[t,s] = bernoulli_rng(inv_logit(b0_b[s]+b1_b[s]*stim[t,s]+b2_b[s]*expect[t,s]));
   //  
   // post_pred[t,s] = bernoulli_rng(belief[t,s]);
   //  
    
     log_lik[t,s,sess] = bernoulli_lpmf(percept_bin[t,s,sess] | inv_logit(b0_b[s,sess]+b1_b[s,sess]*stim[t,s,sess]+b2_b[s,sess]*expect[t,s,sess]))+
                     beta_proportion_lpdf(percept[t,s,sess] | inv_logit(b0[s,sess]+b1[s,sess]*stim[t,s,sess]+b2[s,sess]*expect[t,s,sess]), percept_precision[s,sess])+
                     bernoulli_lpmf(pred[t,s,sess] |  belief[t,s,sess]);
      
      
      
      
      
      }
    }
  
  }
  
  
}

