// This is an "include guard". It prevents the file from being included more than once.
#ifndef TMB_OCCUN_HPP
#define TMB_OCCUN_HPP

// We no longer need to include TMB.hpp here because the main dispatcher file already does.
// This was the source of the error.

template<class Type>
// The function signature is changed to match the dispatcher pattern
Type tmb_occuN(objective_function<Type>* obj) {

  // --- 1. Data and Parameters ---
  // The DATA_ and PARAMETER_ macros now need to use the 'obj' pointer
  DATA_MATRIX(y);
  DATA_MATRIX(X);
  DATA_MATRIX(V);
  DATA_MATRIX(w);

  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);

  Type nll = 0.0;

  // --- 2. State Process ---
  vector<Type> log_lambda_j = X * beta;
  vector<Type> lambda_j = exp(log_lambda_j);

  // --- 3. Link State to Occupancy ---
  vector<Type> lambda_tilde_i = w * lambda_j;
  vector<Type> psi_i = Type(1.0) - exp(-lambda_tilde_i);

  // --- 4. Observation Process ---
  int M = y.rows();
  int J = y.cols();
  vector<Type> logit_p = V * alpha;

  for (int i = 0; i < M; i++) {
    Type log_prob_y_given_occupied = 0.0;
    for (int t = 0; t < J; t++) {
      Type p_it = Type(1.0) / (Type(1.0) + exp(-logit_p(i * J + t)));
      log_prob_y_given_occupied += dbinom(y(i, t), Type(1.0), p_it, true);
    }

    if (y.row(i).sum() > 0) {
      nll -= log(psi_i(i)) + log_prob_y_given_occupied;
    } else {
      Type prob_occupied_missed = psi_i(i) * exp(log_prob_y_given_occupied);
      Type prob_unoccupied = Type(1.0) - psi_i(i);
      nll -= log(prob_occupied_missed + prob_unoccupied);
    }
  }

  return nll;
}

#endif // End of the include guard