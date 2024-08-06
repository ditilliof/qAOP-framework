functions {
  array[] real ode(real t, array[] real y, array[] real theta,
                   array[] real x_r, array[] int x_i) {
    // Define the state variables
    real MIE = y[1];
    real KE1 = y[2];
    real KE2 = y[3];
    
    // Define the parameter names (optional)
    
    real tau1 = theta[1];
    real d1 = theta[2];
    real k_12 = theta[3];
    real h_12 = theta[4];
    real d2 = theta[5];
    
    /** Below you can fill in the differental equations of the model*/
    
    real dMIE1 = -theta[1] * y[1];
    real dKE1 = y[1] - theta[2] * y[2];
    real dKE2 = theta[3] * y[2] / (theta[4] + y[2]) - theta[5] * y[3];
    return {dMIE1, dKE1, dKE2};
  }
}
//Here we define the input data
data {
  int<lower=2> N; // Number of time points
  real t0; // Initial time
  array[N] real ts; // List of time points
  array[N, 3] real y; // Observed states for first dose
  array[N, 3] real<lower=0> sigma; // Known standard deviation of the observations for first dose
  real KE10;
  real KE20;
  
}
transformed data {
  array[0] real x_r;
  array[0] int x_i;
}
// Here we define the parameters and/or initial conditions that we want to estimate
parameters {
  real<lower=0> tau1; //in case we also want to give an upper bound: "real<lower=0, upper = *> tau;", with * being the upper bound
  real<lower=0> d1;
  real<lower=0> k_12;
  real<lower=0> h_12;
  real<lower=0> d2;
  
  real<lower=0> S0; //if commented, needs to be defined in the  data section
}
transformed parameters {
  //here we solve the ODE for the two different dose levels
  array[N, 3] real y_hat; //solution dose 1
  {
    array[5] real theta = {tau1, d1, k_12, h_12, d2};
    array[3] real y0 = {S0, KE10, KE20}; //initial condition
    y_hat = integrate_ode_rk45(ode, y0, t0, ts, theta, x_r, x_i); //solve ODE
  }
}
model {
  // We first define prior distributions for parameters & initial conditions (to be estimated, in this example we estimate all of them)
  
  tau1 ~ normal(0.05, 0.01);
  d1 ~ normal(0.05, 0.01);
  k_12 ~ normal(0.5, 0.025);
  d2 ~ normal(0.04, 0.01);
  h_12 ~ normal(5, 1);
  //pars = c(tau_1=0.05, tau_2 = 0, d1 = 0.05, d3 = 0.06, alpha = 0.8, k_32 = 0.08, k_12 = 0.07, h_12 = 5, d2 = 0.1)
  //c(tau_1=0.05, tau_2 = 0, d1 = 0.05, d3 = 0.06, alpha = 0.2, k_32 = 0.08, k_12 = 0.7, h_12 = 1, d2 = 0.5)
  S0 ~ normal(0.1, 0.05); //IC for KE1, dose 1

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


