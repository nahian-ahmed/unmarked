#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type tmb_occuN(objective_function<Type>* obj) {

  // Data and Parameters
  DATA_MATRIX(y);
  DATA_MATRIX(X);
  DATA_MATRIX(V);
  DATA_MATRIX(w);

  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);

  Type nll = 0.0;

  // Model Logic
  vector<Type> log_lambda_j = X * beta;
  vector<Type> lambda_j = exp(log_lambda_j);
  vector<Type> lambda_tilde_i = w * lambda_j;
  vector<Type> psi_i = Type(1.0) - exp(-lambda_tilde_i);
  int M = y.rows();
  int J = y.cols();
  vector<Type> logit_p = V * alpha;


  for (int i = 0; i < M; i++) {
    Type log_prob_y_given_occupied = 0.0;
    for (int t = 0; t < J; t++) {
      if(y(i,t) == y(i,t)) { // Check for NA
        Type logit_p_it = logit_p(i * J + t); 

        Type log_p_it = -logspace_add(Type(0.0), -logit_p_it);
        Type log_one_minus_p_it = -logit_p_it + log_p_it;

        if (y(i, t) == 1.0) {
            log_prob_y_given_occupied += log_p_it;
        } else {
            log_prob_y_given_occupied += log_one_minus_p_it;
        }
      }
    }

    
    Type lambda_tilde_i_current = lambda_tilde_i(i);
    
    // Calculate log(psi_i) = log(1 - exp(-lambda_tilde_i))
    // logspace_sub(0, -x) = log(exp(0) - exp(-x)) = log(1 - exp(-x))
    // The epsilon is required to prevent log(0) when lambda_tilde_i is exactly 0.
    // Type log_psi_i = logspace_sub(Type(0.0), -(lambda_tilde_i_current + 1e-15));

    // Calculate log(psi_i) = log(1 - exp(-lambda_tilde_i))
    // Use TMB's numerically stable log1m_exp(x) which calculates log(1-exp(x)) for x < 0.
    Type log_psi_i = log1m_exp(-lambda_tilde_i_current);

    // Calculate log(1 - psi_i) = log(exp(-lambda_tilde_i))
    Type log_one_minus_psi_i = -lambda_tilde_i_current;

    
    if (y.row(i).sum() > 0) {
      // Site was occupied and detected
      nll -= log_psi_i + log_prob_y_given_occupied;

    } else {
      // Site was not detected
      Type log_prob_occupied_missed = log_psi_i + log_prob_y_given_occupied;
      
      // Use logspace_add() for robustly adding probabilities in log-space
      nll -= logspace_add(log_prob_occupied_missed, log_one_minus_psi_i);
    }
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this