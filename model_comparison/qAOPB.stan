functions {
  array[] real ode(real t, array[] real y, array[] real theta,
                   array[] real x_r, array[] int x_i) {
    // Define the state variables
    real MIE = y[1];
    real KE1 = y[2];
    real KE2 = y[3];
    
    // Define the parameter names (optional)
    
    real tau = theta[1];
    real eps = theta[2];
    real d1 = theta[3];
    real k_12 = theta[4];
    real d2 = theta[5];
    
    /** Below you can fill in the differental equations of the model*/
    
    real dMIE = -theta[1] * y[1];
    real dKE1 = theta[2] + y[1] - theta[3] * y[2];
    real dKE2 = theta[4] * y[2] - theta[5] * y[3];
    return {dMIE, dKE1, dKE2};
  }
}
//Here we define the input data
data {
  int<lower=2> N; // Number of time points
  real t0; // Initial time
  array[N] real ts; // List of time points
  array[N, 3] real y_dose1; // Observed states for first dose
  //array[N, 3] real y_dose2; // Observed states for second dose
  array[N, 3] real<lower=0> sigma_dose1; // Known standard deviation of the observations for first dose
  //array[N, 3] real<lower=0> sigma_dose2; // Known standard deviation of the observations for second dose
  real KE10;
  real KE20;
  //real<lower=0> S1;
  //real<lower=0> S2;
}
transformed data {
  array[0] real x_r;
  array[0] int x_i;
}
// Here we define the parameters and/or initial conditions that we want to estimate
parameters {
  real<lower=0> tau; //in case we also want to give an upper bound: "real<lower=0, upper = *> tau;", with * being the upper bound
  real<lower=0> eps;
  real<lower=0> d1;
  real<lower=0> k_12;
  real<lower=0> d2;
  real<lower=0> S0; //if commented, needs to be defined in the  data section
  //real<lower=0> S2; //if commented, needs to be defined in the  data section
}
transformed parameters {
  //here we solve the ODE for the two different dose levels
  array[N, 3] real y_hat1; //solution dose 1
  //array[N, 3] real y_hat2; //solution dose 2
  {
    array[5] real theta = {tau, eps, d1, k_12, d2};
    array[3] real y1 = {S0, KE10, KE20}; //initial condition dose 1
    //array[3] real y2 = {S2, KE10, KE20}; //initial condition dose 1
    y_hat1 = integrate_ode_rk45(ode, y1, t0, ts, theta, x_r, x_i); //solve for dose 1
    //y_hat2 = integrate_ode_rk45(ode, y2, t0, ts, theta, x_r, x_i); //solve for dose 2
  }
}
model {
  // We first define prior distributions for parameters & initial conditions (to be estimated, in this example we estimate all of them)
  
  tau ~ normal(0.05, 0.5);
  eps ~ normal(0.2, 0.5);
  d1 ~ normal(0.5, 0.25);
  k_12 ~ normal(0.3, 0.25);
  d2 ~ normal(0.4, 0.25);
  
  S0 ~ normal(0.5, 0.025); //IC for KE1, dose 1
  //S2 ~ normal(1, 0.01); //IC for KE1, dose 2
  
  // Likelihood for dose 1
  for (n in 1 : N) {
    y_dose1[n, 1] ~ normal(y_hat1[n, 1], sigma_dose1[n, 1]);
    y_dose1[n, 2] ~ normal(y_hat1[n, 2], sigma_dose1[n, 2]);
    y_dose1[n, 3] ~ normal(y_hat1[n, 3], sigma_dose1[n, 3]);
  }
  
 /* // Likelihood for dose 2
  for (n in 1 : N) {
    y_dose2[n, 1] ~ normal(y_hat2[n, 1], sigma_dose2[n, 1]);
    y_dose2[n, 2] ~ normal(y_hat2[n, 2], sigma_dose2[n, 2]);
    y_dose2[n, 3] ~ normal(y_hat2[n, 3], sigma_dose2[n, 3]);
  }*/
}
generated quantities {
  vector[N] logLikelihood;
  for (n in 1 : N) {
    logLikelihood[n] = normal_lpdf(y_dose1[n, 1] | y_hat1[n, 1], sigma_dose1[n, 1]) + normal_lpdf(y_dose1[n, 2] | y_hat1[n, 2], sigma_dose1[n, 2]) + normal_lpdf(y_dose1[n, 3] | y_hat1[n, 3], sigma_dose1[n, 3]);
  }
}
