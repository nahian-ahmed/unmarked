// Negative Log-Likelihood function for the occuN model
// This is the statistical core of the model, written in C++ using TMB.

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // --- 1. Data and Parameters ---

  // Data passed from R
  DATA_MATRIX(y); // Detection/non-detection data (sites x surveys)
  DATA_MATRIX(X); // State design matrix (cells x covariates)
  DATA_MATRIX(V); // Detection design matrix (observations x covariates)
  DATA_MATRIX(w); // Weights matrix (sites x cells)

  // Parameters to be estimated by TMB
  PARAMETER_VECTOR(alpha); // Detection coefficients (logit-scale)
  PARAMETER_VECTOR(beta);  // State (abundance) coefficients (log-scale)

  // Initialize the negative log-likelihood (nll)
  Type nll = 0.0;

  // --- 2. State Process (Abundance Intensity) ---

  // Calculate cell-level intensity (lambda) on the log scale
  // This is the linear predictor: log(lambda_j) = X_j * beta
  vector<Type> log_lambda_j = X * beta;

  // Exponentiate to get the intensity for each cell
  // lambda_j = exp(log(lambda_j))
  vector<Type> lambda_j = exp(log_lambda_j);

  // --- 3. Link State (Abundance) to Occupancy ---

  // Calculate site-level expected abundance (lambda_tilde) using the weights matrix
  // This is the core of the occuN model: lambda_tilde_i = sum(w_ij * lambda_j)
  vector<Type> lambda_tilde_i = w * lambda_j;

  // Calculate site-level occupancy probability (psi) using the cloglog link
  // psi_i = 1 - exp(-lambda_tilde_i)
  vector<Type> psi_i = Type(1.0) - exp(-lambda_tilde_i);

  // --- 4. Observation Process and Likelihood Calculation ---

  // Get model dimensions
  int M = y.rows(); // Number of sites
  int J = y.cols(); // Number of surveys per site

  // Calculate the linear predictor for detection probability for all observations
  // logit(p_it) = V_it * alpha
  vector<Type> logit_p = V * alpha;

  // Loop through each site 'i' to calculate its contribution to the likelihood
  for (int i = 0; i < M; i++) {
    
    // Calculate the probability of the observed detection history for site i,
    // given that the site is occupied.
    Type log_prob_y_given_occupied = 0.0;
    for (int t = 0; t < J; t++) {
      // Back-transform from logit scale to probability scale for each observation
      Type p_it = Type(1.0) / (Type(1.0) + exp(-logit_p(i * J + t)));
      
      // Add the log-probability of this specific observation (Bernoulli trial)
      // This uses R's dbinom function in C++: dbinom(observation, size, prob, give_log)
      log_prob_y_given_occupied += dbinom(y(i, t), Type(1.0), p_it, true);
    }

    // Check if the species was ever detected at this site
    if (y.row(i).sum() > 0) {
      // If detected, the site MUST be occupied.
      // The log-likelihood is log(psi) + log(P(detection history | occupied))
      nll -= log(psi_i(i)) + log_prob_y_given_occupied;
    } else {
      // If not detected, there are two possibilities:
      // 1. Site is occupied, but species was not detected (prob = psi * P(no detects))
      // 2. Site is unoccupied (prob = 1 - psi)
      Type prob_occupied_missed = psi_i(i) * exp(log_prob_y_given_occupied);
      Type prob_unoccupied = Type(1.0) - psi_i(i);
      
      // The total probability is the sum of these two mutually exclusive events.
      // The log-likelihood is the log of this total probability.
      nll -= log(prob_occupied_missed + prob_unoccupied);
    }
  }

  return nll;
}
