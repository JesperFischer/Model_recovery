
data {
  int<lower=1> N;  // total number of observations (all subjects)
  array[N] int npx;  // total number of observations per X
  
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
  vector[N] X_logistic = rep_vector(0, N);
  vector[N] mu = rep_vector(0, N);
  

  // alpha = (gm[1]+(tau_u[1] * z_expo[1,]))';
  // beta = exp(gm[2]+(tau_u[2] * z_expo[2,]))';
  lapse = (inv_logit(gm[3]+(tau_u[3] * z_expo[1,])) / 2)';
  

  
  //model
  for (n in 1:N){
      X_logistic[n] = 1/(1+exp(-beta[S_id[n]]*(X[n,2]-alpha[S_id[n]])));
      
      mu[n] = lapse[S_id[n]] + (1 - 2 * lapse[S_id[n]]) * X_logistic[n];
    }
}
  
  



model{
  
  //priors
  target += normal_lpdf(gm[1] |  0 , 50);
  target += normal_lpdf(gm[2] |  0.5 , 3);
  target += normal_lpdf(gm[3] | -4 , 2);
  
  // target += normal_lpdf(gm[4] | -3, 2);
  target += std_normal_lpdf(to_vector(z_expo));
  target += normal_lpdf(tau_u[1] | 0, 20)-normal_lccdf(0 | 0, 20);
  target += normal_lpdf(tau_u[2] | 0, 5)-normal_lccdf(0 | 0, 5);
  target += normal_lpdf(tau_u[3] | 0, 5)-normal_lccdf(0 | 0, 5);

  target += normal_lpdf(alpha | gm[1], tau_u[1]);
  target += normal_lpdf(beta | gm[2], tau_u[2]);


  //target += bernoulli_lpmf(Y | mu);
  target += binomial_lpmf(Y |npx ,mu);
}


generated quantities{

  vector[N] log_lik;
  
  for (n in 1:N){
    log_lik[n] = binomial_lpmf(Y[n] |npx[n] ,mu[n]);
  
  }
  
}
