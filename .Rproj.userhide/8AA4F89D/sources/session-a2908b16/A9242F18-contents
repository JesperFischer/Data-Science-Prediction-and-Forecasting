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
  vector<lower=0>[nsubs] beta;
  
  
  // Group-level parameters
  real<lower=0>mu_sigmaEpsilon;
  real<lower=0>mu_sigmaEta;
  real<lower=0>mu_sigmaPsi;
  
  real <lower=0> sd_sigmaEpsilon;
  real <lower=0> sd_sigmaEta;
  real <lower=0> sd_sigmaPsi;
  
  real <lower=0> sd_beta;
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
      
      // print("eta");
      // print(mu_sigmaEta+sigmaEtaID[s]);
      // print("ep");
      // print(mu_sigmaEpsilon+sigmaEpsilonID[s]);
      // print("psi");
      // print(mu_sigmaPsi+sigmaPsiID[s]);
    for (t in 1:(ntrials)){
      
      if(cues[t,s] == 1){
        exp_muu[t,s] = association[t,s];
      }else{
        exp_muu[t,s] = 1-association[t,s];
      }
      
      
      perceptmuu[t,s] =  (mu_sigmaEpsilon+sigmaEpsilonID[s] * exp_muu[t,s] + (mu_sigmaPsi+sigmaPsiID[s] + exp_var[t,s]) * stim[t,s] ) / 
                        (mu_sigmaEpsilon+sigmaEpsilonID[s] + mu_sigmaPsi+sigmaPsiID[s] + exp_var[t,s]);
                       
                       
                       
      perceptvar[t,s] = (mu_sigmaEpsilon+sigmaEpsilonID[s] * (mu_sigmaPsi+sigmaPsiID[s] + exp_var[t,s]) ) / 
                        (mu_sigmaEpsilon+sigmaEpsilonID[s] + mu_sigmaPsi+sigmaPsiID[s] + exp_var[t,s]) ;    
                        
                        
      
      association[t+1,s] = ((mu_sigmaEpsilon+sigmaEpsilonID[s] + mu_sigmaPsi+sigmaPsiID[s]) * exp_muu[t,s] + (exp_var[t,s] * u[t,s])) / 
                                       (mu_sigmaEpsilon+sigmaEpsilonID[s] + mu_sigmaPsi+sigmaPsiID[s] + exp_var[t,s]) ;
                                       
      exp_var[t+1,s] = ((mu_sigmaEpsilon+sigmaEpsilonID[s] + mu_sigmaPsi+sigmaPsiID[s]) * exp_var[t,s] / 
                        (mu_sigmaEpsilon+sigmaEpsilonID[s] + mu_sigmaPsi+sigmaPsiID[s] + exp_var[t,s])) +  mu_sigmaEta+sigmaEtaID[s];
      
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
    target += lognormal_lpdf(beta[s] | log(20), sd_beta);
    
    
    
    for (t in 1:ntrials){
      target += beta_proportion_lpdf(percept[t,s] | perceptmu[t,s], precision_percept[s]);
      
      target += bernoulli_lpmf(percept_bin[t,s] | (perceptmu[t,s]^beta[s])/((perceptmu[t,s]^beta[s])+(1-perceptmu[t,s])^beta[s]));

      target += bernoulli_lpmf(expectPain[t,s] |  (exp_mu[t,s]^beta[s])/((exp_mu[t,s]^beta[s])+(1-exp_mu[t,s])^beta[s]));
  
    }
    
  }
  
  // Hierarchical Priors
  
  target += lognormal_lpdf(mu_sigmaEta | 0 , 1);
  target += lognormal_lpdf(mu_sigmaPsi | 0 , 1);
  target += lognormal_lpdf(mu_sigmaEpsilon | 0 , 1);
  
  target += lognormal_lpdf(sd_sigmaEta | 0 , 1);
  target += lognormal_lpdf(sd_sigmaPsi | 0 , 1);
  target += lognormal_lpdf(sd_sigmaEpsilon | 0 , 1);
  
  target += exponential_lpdf(sd_precision_percept | 0.1);
  target += exponential_lpdf(sd_beta | 0.1);
  
}


