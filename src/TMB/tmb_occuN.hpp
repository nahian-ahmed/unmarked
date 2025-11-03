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

  // for (int i = 0; i < M; i++) {
  //   Type log_prob_y_given_occupied = 0.0;
  //   for (int t = 0; t < J; t++) {
  //     if(y(i,t) == y(i,t)) {
  //       Type logit_p_it = logit_p(i * J + t); 
  //       // Type p_it = Type(1.0) / (Type(1.0) + exp(-logit_p_it));
        
  //       // log_prob_y_given_occupied += dbinom(y(i, t), Type(1.0), p_it, true);


  //       // NEW: Robust log-likelihood calculation
  //       // This calculates dbinom(y, 1, plogis(logit_p_it), log=TRUE)
  //       // without ever producing Inf.

  //       // log(p) = -log(1 + exp(-logit))
  //       Type log_p_it = -logspace_add(Type(0.0), -logit_p_it);
  //       // log(1-p) = -logit - log(1 + exp(-logit))
  //       Type log_one_minus_p_it = -logit_p_it + log_p_it;

  //       if (y(i, t) == 1.0) {
  //           log_prob_y_given_occupied += log_p_it;
  //       } else {
  //           log_prob_y_given_occupied += log_one_minus_p_it;
  //       }


  //     }
  //   }

  //   if (y.row(i).sum() > 0) {
  //     nll -= log(psi_i(i)) + log_prob_y_given_occupied;
  //   } else {
  //     Type prob_occupied_missed = psi_i(i) * exp(log_prob_y_given_occupied);
  //     Type prob_unoccupied = Type(1.0) - psi_i(i);
  //     nll -= log(prob_occupied_missed + prob_unoccupied);
  //   }
  // }


  for (int i = 0; i < M; i++) {
    Type log_prob_y_given_occupied = 0.0;
    for (int t = 0; t < J; t++) {
      if(y(i,t) == y(i,t)) { // Check for NA
        Type logit_p_it = logit_p(i * J + t); 

        // Robust detection log-likelihood (This part is correct)
        Type log_p_it = -logspace_add(Type(0.0), -logit_p_it);
        Type log_one_minus_p_it = -logit_p_it + log_p_it;

        if (y(i, t) == 1.0) {
            log_prob_y_given_occupied += log_p_it;
        } else {
            log_prob_y_given_occupied += log_one_minus_p_it;
        }
      }
    }

    // --- ROBUST STATE LIKELIHOOD ---
    
    Type lambda_tilde_i_current = lambda_tilde_i(i);
    
    // Calculate log(psi_i) = log(1 - exp(-lambda_tilde_i))
    // Use log1mexp() which is robust.
    //
    // --- THE FIX IS HERE: ---
    // Add a tiny epsilon (1e-15) to prevent log(0) if lambda_tilde_i is exactly 0.
    Type log_psi_i = log1mexp(lambda_tilde_i_current + 1e-15);

    // Calculate log(1 - psi_i) = log(exp(-lambda_tilde_i))
    // This part was already stable.
    Type log_one_minus_psi_i = -lambda_tilde_i_current;

    // --- End robust part ---

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