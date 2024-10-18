
data {
  int<lower=1> N;  // total number of observations (all subjects)
  
  int<lower=0> S; //Total number of subjects
  array[N] int S_id;  // total number of observations (per subjects)

  array[N] int Y;  // response variable
  matrix[N, 2] X;  // design matrix (first column being intercept i.e. 1)

}

parameters {
  // hierarchical group level means 
  vector [3] gm;
  // hierarchical group level deviations
  vector<lower = 0>[3]  tau_u;
  // Subject-level estimate matrix 
  matrix[1, S] z_expo;
  
  vector[S] alpha;
  vector<lower=0>[S] beta;

  }



transformed parameters{
  //getting subject level estimates for ease.
  // vector[S] alpha = rep_vector(0, S);
  // vector<lower = 0>[S] beta = rep_vector(0, S);
  vector<lower = 0, upper = 1>[S] lapse = rep_vector(0, S);
  //trial level parameters
  vector[N] X_normal = rep_vector(0, N);
  vector[N] mu = rep_vector(0, N);
  

  // alpha = (gm[1]+(tau_u[1] * z_expo[1,]))';
  // beta = exp(-2 * (gm[2]+(tau_u[2] * z_expo[1,])))';
  lapse = (inv_logit(gm[3]+(tau_u[3] * z_expo[1,])) / 2)';
  

  

  for (n in 1:N){
    
      X_normal[n] = 0.5+0.5*erf(sqrt(beta[S_id[n]]/2)*(X[n,2]-alpha[S_id[n]]));
      
      mu[n] = lapse[S_id[n]] + (1 - 2 * lapse[S_id[n]]) * X_normal[n];
    }
  }
  
  



model{
  
  //priors
  target += normal_lpdf(gm[1] | 0,50);
  target += normal_lpdf(gm[2] |  -1.3,1.3);
  target += normal_lpdf(gm[3] | -4, 2);
  // 
  // // target += normal_lpdf(gm[4] | -3, 2);
  target += std_normal_lpdf(to_vector(z_expo));
  target += normal_lpdf(tau_u[1] | 0, 50)-normal_lccdf(0 | 0, 50);
  target += normal_lpdf(tau_u[2] | 0, 1.3)-normal_lccdf(0 | 0, 1.3);
  target += normal_lpdf(tau_u[3] | 0, 2)-normal_lccdf(0 | 0, 2);
  

  target += normal_lpdf(alpha | gm[1], tau_u[1]);
  target += normal_lpdf(beta | gm[2], tau_u[2]);


  //target += bernoulli_lpmf(Y | mu);
  target += bernoulli_lpmf(Y |mu);
}


generated quantities{


  vector[N] log_lik;
  
  for (n in 1:N){
    log_lik[n] = bernoulli_lpmf(Y[n] | mu[n]);
  
  }
  
}
