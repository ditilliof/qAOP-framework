functions {
  array[] real ode(real t, array[] real y, array[] real theta,
                   array[] real x_r, array[] int x_i) {
    // Define the state variables
    real MIE2 = y[1];
    real KE3 = y[2];
    real KE2 = y[3];
    
    // Define the parameter names (optional)
    
    real tau2 = theta[1];
    real d3 = theta[2];
    real k_32 = theta[3];
    real d2 = theta[4];
  
    /** Below you can fill in the differental equations of the model*/
    
    real dMIE2 = -theta[1] * y[1];
    real dKE3 = y[1] - theta[2] * y[2];
    real dKE2 = theta[3] * y[2] - theta[4] * y[3];
    return {dMIE2, dKE3, dKE2};
  }
}
//Here we define the input data
data {
  int<lower=2> N; // Number of time points
  real t0; // Initial time
  array[N] real ts; // List of time points
  array[N, 3] real y; // Observed states for first dose
  array[N, 3] real<lower=0> sigma; // Known standard deviation of the observations for first dose
  real KE30;
  real KE20;
  real tau2;
  
}
transformed data {
  array[0] real x_r;
  array[0] int x_i;
}
// Here we define the parameters and/or initial conditions that we want to estimate
parameters {
  //real<lower=0.0001, upper = 0.1> tau2;
  real<lower=0.05> d3;
  real<lower=0> k_32;
  real<lower=0> d2;
  real<lower=0, upper= 0.6> S0; //if commented, needs to be defined in the  data section
}
transformed parameters {
  //here we solve the ODE for the two different dose levels
  array[N, 3] real y_hat; //solution dose 1
  {
    array[4] real theta = {tau2, d3, k_32, d2};
    array[3] real y0 = {S0, KE30, KE20}; //initial condition
    y_hat = integrate_ode_rk45(ode, y0, t0, ts, theta, x_r, x_i); //solve ODE
  }
}
model {
  // We first define prior distributions for parameters & initial conditions (to be estimated, in this example we estimate all of them)
  //tau2 ~ normal(0.01, 0.01);
  d3 ~ normal(0.06, 0.03);
  k_32 ~ normal(0.08, 0.03);
  d2 ~ normal(0.5, 0.05);
  
  
  S0 ~ normal(0.5, 0.1); //IC for KE1, dose 1

  // Likelihood 
  for (n in 1 : N) {
    y[n, 1] ~ normal(y_hat[n, 1], sigma[n, 1]);
    y[n, 2] ~ normal(y_hat[n, 2], sigma[n, 2]);
    y[n, 3] ~ normal(y_hat[n, 3], sigma[n, 3]);
  }
}
generated quantities {
  vector[N] logLikelihood;
  for (n in 1 : N) {
    logLikelihood[n] = normal_lpdf(y[n, 1] | y_hat[n, 1], sigma[n, 1]) + normal_lpdf(y[n, 2] | y_hat[n, 2], sigma[n, 2]) + normal_lpdf(y[n, 3] | y_hat[n, 3], sigma[n, 3]);
  }
}

