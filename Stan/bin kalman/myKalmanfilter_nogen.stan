data {
  
  int<lower=1> nsubs; // number of subjects
  int<lower=1> ntrials; // number of subjects
  matrix[ntrials, nsubs] percept; // observations
  int expectPain[ntrials, nsubs]; // prediciton
  int percept_bin[ntrials, nsubs]; // prediciton
  
  matrix[ntrials, nsubs] stim; // observations
  matrix[ntrials, nsubs] cues; // observations
  matrix[ntrials, nsubs] u; // observations
  
  
}

parameters {
  vector<lower=0>[nsubs] sigmaEta;
  vector<lower=0>[nsubs] sigmaPsi;
  vector<lower=0>[nsubs] sigmaEpsilon;
  
  // Group-level parameters
  real  mu_sigmaEpsilon;
  real  mu_sigmaEta;
  real  mu_sigmaPsi;
  
  real <lower=0> sd_sigmaEpsilon;
  real <lower=0> sd_sigmaEta;
  real <lower=0> sd_sigmaPsi;
  
}

transformed parameters {
  matrix <lower=0, upper  = 1> [ntrials, nsubs] perceptmu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] perceptmuu;
  matrix <lower=0> [ntrials, nsubs] perceptvar;
  
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] association; 
  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] exp_mu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] exp_muu;
  matrix <lower=0> [ntrials+1, nsubs] exp_var;
  

  for (s in 1:nsubs) {
      association[1,s] = 0.5;
      exp_var[1,s] = sigmaEta[s];
    for (t in 1:(ntrials)){
      
      if(cues[t,s] == 1){
        exp_muu[t,s] = association[t,s];
      }else{
        exp_muu[t,s] = 1-association[t,s];
      }
      
      
      perceptmuu[t,s] =  (sigmaEpsilon[s] * exp_muu[t,s] + (sigmaPsi[s] + exp_var[t,s]) * stim[t,s] ) / 
                        (sigmaEpsilon[s] + sigmaPsi[s] + exp_var[t,s]);
                       
                       
                       
      perceptvar[t,s] = (sigmaEpsilon[s] * (sigmaPsi[s] + exp_var[t,s]) ) / 
                        (sigmaEpsilon[s] + sigmaPsi[s] + exp_var[t,s]) ;    
                        
                        
      
      association[t+1,s] = ((sigmaEpsilon[s] + sigmaPsi[s]) * exp_muu[t,s] + (exp_var[t,s] * u[t,s])) / 
                                       (sigmaEpsilon[s] + sigmaPsi[s] + exp_var[t,s]) ;
                                       
      exp_var[t+1,s] = ((sigmaEpsilon[s] + sigmaPsi[s]) * exp_var[t,s] / (sigmaEpsilon[s] + sigmaPsi[s] + exp_var[t,s])) + sigmaEta[s];
      
      if (exp_muu[t,s] > 0.9999){
                exp_mu[t,s] = 0.9999;
            } else if (exp_muu[t,s] < 0.0001) {
                exp_mu[t,s] = 0.0001;
            } else if (exp_muu[t,s] > 0.0001 && exp_muu[t,s] < 0.9999) {
                exp_mu[t,s] = exp_muu[t,s];
            } else {
                exp_mu[t,s] = 0.5;
            }
      
      if (perceptmuu[t,s] > 0.9999){
                perceptmu[t,s] = 0.9999;
            } else if (perceptmuu[t,s] < 0.0001) {
                perceptmu[t,s] = 0.0001;
            } else if (perceptmuu[t,s] > 0.0001 && perceptmuu[t,s] < 0.9999) {
                perceptmu[t,s] = perceptmuu[t,s];
            } else {
                perceptmu[t,s] = 0.5;
            }
      
      
    }
  }
}

model {
   for (s in 1:nsubs){

    target += lognormal_lpdf(sigmaEta[s] | mu_sigmaEta, sd_sigmaEta);
    target += lognormal_lpdf(sigmaPsi[s] | mu_sigmaPsi, sd_sigmaPsi);
    target += lognormal_lpdf(sigmaEpsilon[s] | mu_sigmaEpsilon, sd_sigmaEpsilon);
    
    
    
    for (t in 1:ntrials){
      target += beta_proportion_lpdf(percept[t,s] | perceptmu[t,s], 1/perceptvar[t,s]);
      
      target += bernoulli_lpmf(percept_bin[t,s] | (perceptmu[t,s]^(1/perceptvar[t,s]))/((perceptmu[t,s]^(1/perceptvar[t,s]))+(1-perceptmu[t,s])^(1/perceptvar[t,s])));

      target += bernoulli_lpmf(expectPain[t,s] |  (exp_mu[t,s]^(1/exp_var[t,s]))/((exp_mu[t,s]^(1/exp_var[t,s]))+(1-exp_mu[t,s])^(1/exp_var[t,s])));
  
    }
    
  }
  
  // Hierarchical Priors
  
  target += normal_lpdf(mu_sigmaEta | 0 , 2);
  target += normal_lpdf(mu_sigmaPsi | 0 , 2);
  target += normal_lpdf(mu_sigmaEpsilon | 0 , 2);
  
  target += lognormal_lpdf(sd_sigmaEta | 0 , 1);
  target += lognormal_lpdf(sd_sigmaPsi | 0 , 1);
  target += lognormal_lpdf(sd_sigmaEpsilon | 0 , 1);
  
}

