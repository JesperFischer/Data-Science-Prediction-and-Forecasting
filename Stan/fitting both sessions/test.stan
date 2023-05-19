
data {
  int<lower=0> N;
  vector[N] y;
  
  int<lower=1> nsubs; // number of subjects
  int<lower=1> ntrials; // number of subjects
  int<lower=1> sesions;
  array[ntrials, nsubs,sesions] real  percept; // observations
  array[ntrials, nsubs,sesions] int expectPain; // prediciton
  array[ntrials, nsubs,sesions] int percept_bin; // prediciton
  
  array[ntrials, nsubs,sesions] int stim; // observations
  array[ntrials, nsubs,sesions] int cues; // observations

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
  
  array[nsubs,sesions] real <lower=0, upper = 1> alpha;
  array[nsubs,sesions] real precision_percept;
  array[nsubs,sesions] real beta;
  
  array[nsubs,sesions] real w1;

  // Group-level parameters
  array[sesions] real kappa_alpha;
  array[sesions] real mu_alpha;
  array[sesions] real mu_w1;
  array[sesions] real kappa_w1;
 // Group-level parameters
  array[sesions] real sd_beta;
  array[sesions] real sd_precision_percept;
}

transformed parameters{
  array[ntrials, nsubs,sesions] real  painMu;
  array[ntrials+1, nsubs,sesions] real association;

  
  array[9, 10, 11] real w;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y ~ normal(mu, sigma);
}

