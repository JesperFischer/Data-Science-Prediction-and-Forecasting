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
  vector[nsubs] sigmaEtaID_Z;
  vector[nsubs] sigmaPsiID_Z;
  vector[nsubs] sigmaEpsilonID_Z;
  
  vector<lower=0>[nsubs] precision_percept;
  
  
  
  // Group-level parameters
  real<lower=0>mu_sigmaEpsilon;
  real<lower=0>mu_sigmaEta;
  real<lower=0>mu_sigmaPsi;
  
  real <lower=0> sd_sigmaEpsilon;
  real <lower=0> sd_sigmaEta;
  real <lower=0> sd_sigmaPsi;
  
  
  real <lower=0> sd_precision_percept;
  
}

transformed parameters {
  matrix <lower=0, upper  = 1> [ntrials, nsubs] perceptmu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] perceptmuu;
  matrix <lower=0> [ntrials, nsubs] perceptvar;
  
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] association; 
  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] exp_mu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] exp_muu;
  matrix <lower=0> [ntrials+1, nsubs] exp_var;
  
  
  vector[nsubs] sigmaEtaID;
  vector[nsubs] sigmaPsiID;
  vector[nsubs] sigmaEpsilonID;
  
  
  sigmaEtaID = exp(sigmaEtaID_Z)*sd_sigmaEta;
  sigmaPsiID = exp(sigmaPsiID_Z)*sd_sigmaPsi;
  sigmaEpsilonID = exp(sigmaEpsilonID_Z)*sd_sigmaEpsilon;
  


  for (s in 1:nsubs) {
      association[1,s] = 0.5;
      exp_var[1,s] = mu_sigmaEta+sigmaEtaID[s];
      
    for (t in 1:(ntrials)){
      
      if(cues[t,s] == 1){
        exp_muu[t,s] = association[t,s];
      }else{
        exp_muu[t,s] = 1-association[t,s];
      }
      
      
      perceptmuu[t,s] =  ((mu_sigmaEpsilon+sigmaEpsilonID[s]) * exp_muu[t,s] + ((mu_sigmaPsi+sigmaPsiID[s]) + exp_var[t,s]) * stim[t,s] ) / 
                        ((mu_sigmaEpsilon+sigmaEpsilonID[s]) + (mu_sigmaPsi+sigmaPsiID[s]) + exp_var[t,s]);
                       
                       
                       
      perceptvar[t,s] = ((mu_sigmaEpsilon+sigmaEpsilonID[s]) * ((mu_sigmaPsi+sigmaPsiID[s]) + exp_var[t,s]) ) / 
                        ((mu_sigmaEpsilon+sigmaEpsilonID[s]) + (mu_sigmaPsi+sigmaPsiID[s]) + exp_var[t,s]) ;    
                        
                        
      
      association[t+1,s] = (((mu_sigmaEpsilon+sigmaEpsilonID[s]) + (mu_sigmaPsi+sigmaPsiID[s])) * exp_muu[t,s] + (exp_var[t,s] * u[t,s])) / 
                                       ((mu_sigmaEpsilon+sigmaEpsilonID[s]) + (mu_sigmaPsi+sigmaPsiID[s]) + exp_var[t,s]) ;
                                       
      exp_var[t+1,s] = (((mu_sigmaEpsilon+sigmaEpsilonID[s]) + (mu_sigmaPsi+sigmaPsiID[s])) * exp_var[t,s] / 
                        ((mu_sigmaEpsilon+sigmaEpsilonID[s]) + (mu_sigmaPsi+sigmaPsiID[s]) + exp_var[t,s])) +  (mu_sigmaEta+sigmaEtaID[s]);
      
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

    target += std_normal_lpdf(to_vector(sigmaEtaID_Z));
    target += std_normal_lpdf(to_vector(sigmaPsiID_Z)); 
    target += std_normal_lpdf(to_vector(sigmaEpsilonID_Z)); 
    

    
    
    target += lognormal_lpdf(precision_percept[s] | log(20), sd_precision_percept);

    
    
    for (t in 1:ntrials){
      target += beta_proportion_lpdf(percept[t,s] | perceptmu[t,s], precision_percept[s]);
      
      target += bernoulli_lpmf(percept_bin[t,s] | perceptmu[t,s]);

      target += bernoulli_lpmf(expectPain[t,s] |  exp_mu[t,s]);
  
    }
    
  }
  
  // Hierarchical Priors
  
  target += lognormal_lpdf(mu_sigmaEta | log(20),0.6);
  target += lognormal_lpdf(mu_sigmaPsi | log(20),0.6);
  target += lognormal_lpdf(mu_sigmaEpsilon | log(20),0.6);
  
  target += exponential_lpdf(sd_sigmaEta | 1);
  target += exponential_lpdf(sd_sigmaPsi | 1);
  target += exponential_lpdf(sd_sigmaEpsilon | 1);
  
  target += exponential_lpdf(sd_precision_percept | 1);

}




generated quantities{
  
  
  vector[nsubs] prior_sigmaEtaID;
  vector[nsubs] prior_sigmaPsiID;
  vector[nsubs] prior_sigmaEpsilonID;
  
  
  vector[nsubs] prior_sigmaEtaID_Z;
  vector[nsubs] prior_sigmaPsiID_Z;
  vector[nsubs] prior_sigmaEpsilonID_Z;
  
  vector<lower=0>[nsubs] prior_precision_percept;

  real  prior_mu_sigmaEpsilon;
  real  prior_mu_sigmaEta;
  real  prior_mu_sigmaPsi;
  
  real <lower=0> prior_sd_sigmaEpsilon;
  real <lower=0> prior_sd_sigmaEta;
  real <lower=0> prior_sd_sigmaPsi;

  real <lower=0> prior_sd_precision_percept;

  
    //trial level
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_perceptmu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_perceptmuu;
  matrix <lower=0> [ntrials, nsubs] prior_perceptvar;
  
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] prior_association; 
  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_exp_mu;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_exp_muu;
  matrix <lower=0> [ntrials+1, nsubs] prior_exp_var;
  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_percept;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_percept_bin;  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] prior_expectPain;
  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] post_percept;
  matrix <lower=0, upper  = 1> [ntrials, nsubs] post_percept_bin;  
  matrix <lower=0, upper  = 1> [ntrials, nsubs] post_expectPain;
  
  
  
  matrix[ntrials, nsubs] log_lik;
  
  
  prior_mu_sigmaEta = lognormal_rng(log(20),0.6);
  prior_mu_sigmaPsi = lognormal_rng(log(20),0.6);
  prior_mu_sigmaEpsilon = lognormal_rng(log(20),0.6);
  
  prior_sd_sigmaEta = exponential_rng(1);
  prior_sd_sigmaPsi = exponential_rng(1);
  prior_sd_sigmaEpsilon = exponential_rng(1);
  
  
  prior_sd_precision_percept = exponential_rng(1);

  
  
  
  
    for (s in 1:nsubs) {
      prior_sigmaEtaID_Z[s] = std_normal_rng();
      prior_sigmaPsiID_Z[s] = std_normal_rng();
      prior_sigmaEpsilonID_Z[s] = std_normal_rng();
      
      
      
      prior_sigmaEtaID[s] = exp(prior_sigmaEtaID_Z[s])*prior_sd_sigmaEta;
      prior_sigmaPsiID[s] = exp(prior_sigmaPsiID_Z[s])*prior_sd_sigmaPsi;
      prior_sigmaEpsilonID[s] = exp(prior_sigmaEpsilonID_Z[s])*prior_sd_sigmaEpsilon;

      prior_precision_percept[s] = lognormal_rng(log(20), prior_sd_precision_percept);

    
      prior_association[1,s] = 0.5;
      prior_exp_var[1,s] = prior_mu_sigmaEta+prior_sigmaEtaID[s];
      
    for (t in 1:(ntrials)){
      
      if(cues[t,s] == 1){
        prior_exp_muu[t,s] = prior_association[t,s];
      }else{
        prior_exp_muu[t,s] = 1-prior_association[t,s];
      }
      
      
      prior_perceptmuu[t,s] =  ((prior_mu_sigmaEpsilon+prior_sigmaEpsilonID[s]) * prior_exp_muu[t,s] + ((prior_mu_sigmaPsi+prior_sigmaPsiID[s]) + prior_exp_var[t,s]) * stim[t,s] ) / 
                        ((prior_mu_sigmaEpsilon+prior_sigmaEpsilonID[s]) + (prior_mu_sigmaPsi+prior_sigmaPsiID[s]) + prior_exp_var[t,s]);
                       
                       
                       
      prior_perceptvar[t,s] = ((prior_mu_sigmaEpsilon+prior_sigmaEpsilonID[s]) * ((prior_mu_sigmaPsi+prior_sigmaPsiID[s]) + prior_exp_var[t,s]) ) / 
                        ((prior_mu_sigmaEpsilon+prior_sigmaEpsilonID[s]) + prior_mu_sigmaPsi + prior_exp_var[t,s]) ;    
                        
                        
      
      prior_association[t+1,s] = (((prior_mu_sigmaEpsilon+prior_sigmaEpsilonID[s]) + (prior_mu_sigmaPsi+prior_sigmaPsiID[s])) * prior_exp_muu[t,s] + (prior_exp_var[t,s] * u[t,s])) / 
                                       ((prior_mu_sigmaEpsilon+prior_sigmaEpsilonID[s]) +(prior_mu_sigmaPsi+prior_sigmaPsiID[s]) + prior_exp_var[t,s]);
                                       
      prior_exp_var[t+1,s] = ((((prior_mu_sigmaEpsilon+prior_sigmaEpsilonID[s]) + (prior_mu_sigmaPsi+prior_sigmaPsiID[s])) * prior_exp_var[t,s]) /
                              ((prior_mu_sigmaEpsilon+prior_sigmaEpsilonID[s]) + (prior_mu_sigmaPsi+prior_sigmaPsiID[s]) + prior_exp_var[t,s])) + (prior_mu_sigmaEta+prior_sigmaEtaID[s]);
      
      if (prior_exp_muu[t,s] > 0.9999){
                prior_exp_mu[t,s] = 0.9999;
            } else if (prior_exp_muu[t,s] < 0.0001) {
                prior_exp_mu[t,s] = 0.0001;
            } else if (prior_exp_muu[t,s] > 0.0001 && prior_exp_muu[t,s] < 0.9999) {
                prior_exp_mu[t,s] = prior_exp_muu[t,s];
            } else {
                prior_exp_mu[t,s] = 0.5;
            }
      
      if (prior_perceptmuu[t,s] > 0.9999){
                prior_perceptmu[t,s] = 0.9999;
            } else if (prior_perceptmuu[t,s] < 0.0001) {
                prior_perceptmu[t,s] = 0.0001;
            } else if (prior_perceptmuu[t,s] > 0.0001 && prior_perceptmuu[t,s] < 0.9999) {
                prior_perceptmu[t,s] = prior_perceptmuu[t,s];
            } else {
                prior_perceptmu[t,s] = 0.5;
            }
      
      prior_percept[t,s] = beta_proportion_rng(prior_perceptmu[t,s], prior_precision_percept[s]);
      
      prior_percept_bin[t,s] = bernoulli_rng(prior_perceptmu[t,s]);

      prior_expectPain[t,s] = bernoulli_rng(prior_exp_mu[t,s]);
      
      post_percept[t,s] = beta_proportion_rng(perceptmu[t,s], precision_percept[s]);
      
      post_percept_bin[t,s] = bernoulli_rng(perceptmu[t,s]);

      post_expectPain[t,s] = bernoulli_rng(exp_mu[t,s]);
      
      
      
      log_lik[t,s] = bernoulli_lpmf(percept_bin[t,s] | perceptmu[t,s])+
                     beta_proportion_lpdf(percept[t,s] | perceptmu[t,s], precision_percept[s])+
                     bernoulli_lpmf(expectPain[t,s] |  exp_mu[t,s]);
      
    }
  }
  
  
  
}
