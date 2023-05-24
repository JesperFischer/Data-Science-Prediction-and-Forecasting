data {
  int<lower=1> sessions;
  int<lower=1> ntrials; //number of ntrials
  int<lower=1> nsubs; //number of nsubs
  
  array[ntrials,nsubs,sessions] int pred;
  array[ntrials,nsubs,sessions] int u;
  

  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  
  array[nsubs,sessions] real <lower = 0, upper = 1> alpha;
  
  array[sessions] real <lower = 0, upper = 1> mu_alpha;
  array[sessions] real <lower = 0, upper = 1> sd_alpha;
  
  real <lower = 0, upper = 1> mu_mu_alpha;
  real <lower = 0, upper = 1> sd_mu_alpha;
  real <lower = 0, upper = 1> sd_sd_alpha;
  real <lower = 0, upper = 1> mu_sd_alpha;
  
}


transformed parameters{
  array[ntrials,nsubs,sessions] real <lower = 0, upper = 1> belief;
  
  for (ss in 1:sessions){
    for (s in 1:nsubs){
      belief[1,s,ss] = 0.5;
      
  
      for (t in 2:ntrials){
        
        
        belief[t,s,ss] = belief[t-1,s,ss]+alpha[s,ss]*(u[t-1,s,ss]-belief[t-1,s,ss]);
        

      }
    }
  }
  
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  target += beta_lpdf(mu_mu_alpha | 1,1);
  target += beta_lpdf(sd_mu_alpha | 1,1);
    
      
  target += beta_lpdf(sd_sd_alpha | 1,1);
  target += beta_lpdf(mu_sd_alpha | 1,1);
    
  
  for (ss in 1:sessions){
    target += beta_lpdf(mu_alpha[ss] | mu_mu_alpha,sd_mu_alpha);
    target += beta_lpdf(sd_alpha[ss] | mu_sd_alpha,sd_sd_alpha);
    
    
    for (s in 1:nsubs){
      
      target +=beta_lpdf(alpha[s,ss] | mu_alpha[ss],sd_alpha[ss]);
      
      
      for (t in 1:ntrials){
      
        target +=bernoulli_lpmf(pred[t,s,ss] | belief[t,s,ss]);

      }
    
    }
  }
}



generated quantities{
  
  
  
  
}
