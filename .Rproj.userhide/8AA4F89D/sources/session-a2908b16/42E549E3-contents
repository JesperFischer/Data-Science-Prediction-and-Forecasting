data {
  int<lower=1> ntrials; //number of ntrials
  int<lower=1> nsubs; //number of nsubs
  
  array[ntrials,nsubs] int pred;
  array[ntrials,nsubs] int u;
  matrix[ntrials, nsubs] percept; // observations
  
  matrix[ntrials, nsubs] stim; // observations
  matrix[ntrials, nsubs] cues; // observations
  int percept_bin[ntrials, nsubs];

  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  array[nsubs] real  b0;
  array[nsubs] real  b1;
  array[nsubs] real  b2;
  array[nsubs] real  b0_b;
  array[nsubs] real  b1_b;
  array[nsubs] real  b2_b;
  array[nsubs] real  b3_b;
  
  array[nsubs] real <lower = 0, upper = 1> alpha;
  
  vector<lower=0>[nsubs] precision_percept;
  vector<lower=0>[nsubs] beta;
  
  
  real <lower = 0, upper = 1> mu_alpha;
  real <lower = 0, upper = 1> sd_alpha;
  
  
  real <lower=0> sd_beta;
  real <lower=0> sd_precision_percept;
  
  real mu_b0;
  real <lower = 0> sd_b0;
  
  real  mu_b1;
  real <lower = 0> sd_b1;
  
  real mu_b2;
  real <lower = 0> sd_b2;
  
  real mu_b0_b;
  real <lower = 0> sd_b0_b;
  
  real  mu_b1_b;
  real <lower = 0> sd_b1_b;
  
  real mu_b2_b;
  real <lower = 0> sd_b2_b;
  
  real mu_b3_b;
  real <lower = 0> sd_b3_b;
  
}


transformed parameters{
  array[ntrials,nsubs] real <lower = 0, upper = 1> belief;
  array[ntrials,nsubs] real <lower = 0, upper = 1> expect;
  

  for (s in 1:nsubs){
    belief[1,s] = 0.5;
    expect[1,s] = 0.5;

    for (t in 2:ntrials){
      
      
      belief[t,s] = belief[t-1,s]+alpha[s]*(u[t-1,s]-belief[t-1,s]);
      
      
      if(cues[t,s] == 1){
        expect[t,s] = belief[t,s];
      }else{
        expect[t,s] = 1-belief[t,s];
      }
    
    }
  }
  
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target += beta_lpdf(mu_alpha | 1,1);
  target += beta_lpdf(sd_alpha | 1,1);
  
  target +=normal_lpdf(mu_b0 | 0,1);
  target +=lognormal_lpdf(sd_b0 | 0,0.5);
  
  target +=normal_lpdf(mu_b1 | 0,1);
  target +=lognormal_lpdf(sd_b1 | 0,0.5);
  
  target +=normal_lpdf(mu_b2 | 0,1);
  target +=lognormal_lpdf(sd_b2 | 0,0.5);
  
  
  target +=normal_lpdf(mu_b0_b | 0,1);
  target +=lognormal_lpdf(sd_b0_b | 0,0.5);
  
  target +=normal_lpdf(mu_b1_b | 0,1);
  target +=lognormal_lpdf(sd_b1_b | 0,0.5);
  
  target +=normal_lpdf(mu_b2_b | 0,1);
  target +=lognormal_lpdf(sd_b2_b | 0,0.5);
  
  target +=normal_lpdf(mu_b3_b | 0,1);
  target +=lognormal_lpdf(sd_b3_b | 0,0.5);
  
  
  target += exponential_lpdf(sd_precision_percept | 1);
  target += exponential_lpdf(sd_beta | 1);
  
    
  for (s in 1:nsubs){
    
    
    target +=beta_lpdf(alpha[s] | mu_alpha,sd_alpha);
    
    target += lognormal_lpdf(precision_percept[s] | log(20), sd_precision_percept);
    target += lognormal_lpdf(beta[s] | log(20), sd_beta);

    
    target +=normal_lpdf(b0[s] | mu_b0,sd_b0);
    target +=normal_lpdf(b1[s] | mu_b1,sd_b1);
    target +=normal_lpdf(b2[s] | mu_b2,sd_b2);
    
    target +=normal_lpdf(b0_b[s] | mu_b0_b,sd_b0_b);
    target +=normal_lpdf(b1_b[s] | mu_b1_b,sd_b1_b);
    target +=normal_lpdf(b2_b[s] | mu_b2_b,sd_b2_b);
    
    
    for (t in 1:ntrials){
    
      target +=bernoulli_lpmf(pred[t,s] | (belief[t,s]^beta[s])/((belief[t,s]^beta[s])+(1-belief[t,s])^(beta[s])));
    ,  target +=beta_proportion_lpdf(percept[t,s] | inv_logit(b0[s]+b1[s]*stim[t,s]+b2[s]*expect[t,s]), precision_percept[s]);
      target += bernoulli_lpmf(percept_bin[t,s] |  inv_logit(b0_b[s]+b1_b[s]*percept[t,s]+b2_b[s]*stim[t,s]+b3_b[s]*expect[t,s]));

    }
    
}

}



generated quantities{
  
  array[nsubs] real  prior_b0;
  array[nsubs] real  prior_b1;
  array[nsubs] real  prior_b2;
  array[nsubs] real  prior_b0_b;
  array[nsubs] real  prior_b1_b;
  array[nsubs] real  prior_b2_b;
  array[nsubs] real  prior_b3_b;
  
  array[nsubs] real <lower = 0, upper = 1> prior_alpha;
  
  vector<lower=0>[nsubs] prior_precision_percept;
  vector<lower=0>[nsubs] prior_beta;
  
  
  real <lower = 0, upper = 1> prior_mu_alpha;
  real <lower = 0, upper = 1> prior_sd_alpha;
  
  
  real <lower=0> prior_sd_beta;
  real <lower=0> prior_sd_precision_percept;
  
  real prior_mu_b0;
  real <lower = 0> prior_sd_b0;
  
  real  prior_mu_b1;
  real <lower = 0> prior_sd_b1;
  
  real prior_mu_b2;
  real <lower = 0> prior_sd_b2;
  
  real prior_mu_b0_b;
  real <lower = 0> prior_sd_b0_b;
  
  real  prior_mu_b1_b;
  real <lower = 0> prior_sd_b1_b;
  
  real prior_mu_b2_b;
  real <lower = 0> prior_sd_b2_b;
  
  real prior_mu_b3_b;
  real <lower = 0> prior_sd_b3_b;
  
  array[ntrials,nsubs] real <lower = 0, upper = 1> prior_belief;
  array[ntrials,nsubs] real <lower = 0, upper = 1> prior_expect;
  
  matrix[ntrials, nsubs] log_lik;
  
  
  mu_alpha = beta_rng(1,1);
  sd_alpha = beta_rng(1,1);
  
  mu_b0 =normal_rng(0,1);
  sd_b0 =lognormal_rng(0,0.5);
  
  mu_b1 =normal_rng(0,1);
  sd_b1 =lognormal_rng(0,0.5);
  
  mu_b2 =normal_rng(0,1);
  sd_b2 =lognormal_rng(0,0.5);
  
  
  mu_b0_b =normal_rng(0,1);
  sd_b0_b =lognormal_rng(0,0.5);
  
  mu_b1_b =normal_rng(0,1);
  sd_b1_b =lognormal_rng(0,0.5);
  
  mu_b2_b =normal_rng(0,1);
  sd_b2_b =lognormal_rng(0,0.5);
  
  mu_b3_b =normal_rng(0,1);
  sd_b3_b =lognormal_rng(0,0.5);
  
  
  
  for (s in 1:nsubs){
    belief_prior[1,s] = 0.5;
    expect_prior[1,s] = 0.5;

    for (t in 2:ntrials){
      
      
      belief_prior[t,s] = belief_prior[t-1,s]+alpha_prior[s]*(u[t-1,s]-belief_prior[t-1,s]);
      
      
      if(cues[t,s] == 1){
        expect_prior[t,s] = belief_prior[t,s];
      }else{
        expect_prior[t,s] = 1-belief_prior[t,s];
      }
    
    

    
    
   prior_percept[t,s] = beta_proportion_rng(inv_logit(prior_b0[s]+prior_b1[s]*stim[t,s]+prior_b2[s]*expect[t,s]), prior_precision_percept[s]));
          
   prior_percept_bin[t,s] = bernoulli_rng(inv_logit(prior_b0_b[s]+prior_b1_b[s]*percept[t,s]+prior_b2_b[s]*stim[t,s]+prior_b3_b[s]*prior_expect[t,s]));
    
   prior_pred[t,s] = bernoulli_rng((prior_belief[t,s]^prior_beta[s])/((prior_belief[t,s]^prior_beta[s])+(1-prior_belief[t,s])^(prior_beta[s])));
    
    
   post_percept[t,s] = beta_proportion_rng(inv_logit(b0[s]+b1[s]*stim[t,s]+b2[s]*expect[t,s]), precision_percept[s]));
          
   post_percept_bin[t,s] = bernoulli_rng(inv_logit(b0_b[s]+b1_b[s]*percept[t,s]+b2_b[s]*stim[t,s]+b3_b[s]*expect[t,s]));
    
   post_pred[t,s] = bernoulli_rng((belief[t,s]^beta[s])/((belief[t,s]^beta[s])+(1-belief[t,s])^(beta[s])));
    
    
   log_lik[t,s] = bernoulli_lpmf(percept_bin[t,s] | inv_logit(b0_b[s]+b1_b[s]*percept[t,s]+b2_b[s]*stim[t,s]+b3_b[s]*expect[t,s]))+
                   beta_proportion_lpdf(percept[t,s] | inv_logit(b0[s]+b1[s]*stim[t,s]+b2[s]*expect[t,s]), precision_percept[s])+
                   bernoulli_lpmf(pred[t,s] |  (belief[t,s]^beta[s])/((belief[t,s]^beta[s])+(1-belief[t,s])^(beta[s])));
    
    
    
    
    
    }
  }
  
  
  
  
}

