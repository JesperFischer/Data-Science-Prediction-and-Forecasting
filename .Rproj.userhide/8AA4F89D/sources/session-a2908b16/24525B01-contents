data {
  int<lower=1> ntrials; //number of ntrials
  int<lower=1> nsubs; //number of nsubs
  
  array[ntrials,nsubs] int pred;
  array[ntrials,nsubs] int u;
  

  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  
  array[nsubs] real <lower = 0, upper = 1> alpha;
  
  real <lower = 0, upper = 1> mu_alpha;
  real <lower = 0, upper = 1> sd_alpha;
  
}


transformed parameters{
  array[ntrials,nsubs] real <lower = 0, upper = 1> belief;
  
    for (s in 1:nsubs){
      belief[1,s] = 0.5;
      
  
      for (t in 2:ntrials){
        
        
        belief[t,s] = belief[t-1,s]+alpha[s]*(u[t-1,s]-belief[t-1,s]);
        

      }
    }
  
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  
    target += beta_lpdf(mu_alpha | 1,1);
    target += beta_lpdf(sd_alpha | 1,1);
    
    
    for (s in 1:nsubs){
      
      target +=beta_lpdf(alpha[s] | mu_alpha,sd_alpha);
      
      
      for (t in 1:ntrials){
      
        target +=bernoulli_lpmf(pred[t,s] | belief[t,s]);

      }
    
  }
}
