data {
  int<lower=1> ntrials; //number of ntrials
  int<lower=1> nsubs; //number of nsubs
  array[ntrials,nsubs] int pred;
  
  array[ntrials,nsubs] int u;
  
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  
  array[nsubs] real <lower = 0, upper = 1> bias_1;
  array[nsubs] real <lower = 0, upper = 1> alpha_1;
  real <lower = 0, upper = 1> bias_1_mu;
  real <lower = 0, upper = 1> bias_1_sd;
  
  real <lower = 0, upper = 1> alpha_1_mu;
  real <lower = 0, upper = 1> alpha_1_sd;
  
  
}


transformed parameters{
  array[ntrials,nsubs] real <lower = 0, upper = 1> belief_1;

  for (s in 1:nsubs){
    belief_1[1,s] = bias_1[s];

    for (t in 2:ntrials){
      belief_1[t,s] = belief_1[t-1,s]+alpha_1[s]*(u[t-1,s]-belief_1[t-1,s]);
  }
}
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target += beta_lpdf(alpha_1_mu | 1,1);
  target += beta_lpdf(alpha_1_sd | 1,1);
  
  
  target +=beta_lpdf(bias_1_mu | 1,1);
  target +=beta_lpdf(bias_1_sd | 1,1);
  
    
  for (s in 1:nsubs){
    target +=beta_lpdf(bias_1[s] | bias_1_mu, bias_1_sd);
    
    target +=beta_lpdf(alpha_1[s] | alpha_1_mu,alpha_1_sd);
    
    for (t in 1:ntrials){
    
      target +=bernoulli_lpmf(pred[t,s] | belief_1[t,s]);
    }
    
}

}



