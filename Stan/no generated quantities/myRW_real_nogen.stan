data {

  int<lower=1> nsubs; // number of subjects
  int<lower=1> ntrials; // number of subjects
  matrix[ntrials, nsubs] percept; // observations
  int expectPain[ntrials, nsubs]; // prediciton
  int percept_bin[ntrials, nsubs]; // prediciton
  
  matrix[ntrials, nsubs] stim; // observations
  matrix[ntrials, nsubs] cues; // observations
}

parameters {
  vector<lower=0,upper=1>[nsubs] alpha;
  vector<lower=0>[nsubs] precision_percept;
  vector<lower=0>[nsubs] beta;
  
  vector<lower=0, upper=1>[nsubs] w1;

  // Group-level parameters
  real <lower=0> kappa_alpha;
  real <lower=0, upper  = 1> mu_alpha;
  real <lower=0, upper = 1> mu_w1;
  real <lower=0> kappa_w1;
 // Group-level parameters
  real <lower=0> sd_beta;
  real <lower=0> sd_precision_percept;
}

transformed parameters {
  matrix <lower=0, upper  = 1> [ntrials, nsubs] painMu; 
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] association; 
  matrix <lower=0, upper  = 1> [ntrials+1, nsubs] expectMu;
  matrix[ntrials, nsubs] predErr;


  for (s in 1:nsubs) {
    association[1, s] = 0.5;
    expectMu[161, s] = 0.5;
      
    for (t in 1:ntrials){
      
      
      if(cues[t,s] == 1){
        expectMu[t,s] = association[t,s];
      }else{
        expectMu[t,s] = 1-association[t,s];
      }

     if(w1[s]*stim[t,s]+(1-w1[s])*expectMu[t,s] == 0){
        painMu[t,s] = 0.01;
      }else if (w1[s]*stim[t,s]+(1-w1[s])*expectMu[t,s] == 1){
        painMu[t,s] = 0.99;
      }else{
        painMu[t,s] = w1[s]*stim[t,s]+(1-w1[s])*expectMu[t,s];
        }
    
      
      if(cues[t,s] == 1){
        predErr[t,s] = (painMu[t,s] - expectMu[t,s]);
      }else{
        predErr[t,s] = -(painMu[t,s] - expectMu[t,s]);
      }
      
      
      association[t+1,s] = association[t,s] + alpha[s] * predErr[t,s];
    
    }
  }
}


model {
  for (s in 1:nsubs){

    for (t in 1:ntrials){
      target += beta_proportion_lpdf(percept[t,s] | painMu[t,s], precision_percept[s]);
      
      target += bernoulli_lpmf(percept_bin[t,s] | (painMu[t,s]^beta[s])/((painMu[t,s]^beta[s])+(1-painMu[t,s])^(beta[s])));

      target += bernoulli_lpmf(expectPain[t,s] |  (expectMu[t,s]^beta[s])/((expectMu[t,s]^beta[s])+(1-expectMu[t,s])^(beta[s])));
      
    }
    
    target += beta_proportion_lpdf(alpha[s] | mu_alpha , kappa_alpha);
    target += beta_proportion_lpdf(w1[s] | mu_w1 , kappa_w1);

    target += lognormal_lpdf(precision_percept[s] | log(10), sd_precision_percept);
    target += lognormal_lpdf(beta[s] | log(10), sd_beta);
    
  }
  
  // Hierarchical Priors
  target += beta_proportion_lpdf(mu_alpha | 0.3 , 3) ; 
  target += lognormal_lpdf(kappa_alpha | log(30) , 0.5); 
  
  target += beta_proportion_lpdf(mu_w1 | 0.3 , 3) ; 
  target += lognormal_lpdf(kappa_w1 | log(30) , 0.5);
  
  target += exponential_lpdf(sd_precision_percept | 1);
  target += exponential_lpdf(sd_beta | 1);
}

