// // Include guard to prevent multiple inclusions
// #ifndef TMB_OCCUN_HPP
// #define TMB_OCCUN_HPP

// template<class Type>
// Type tmb_occuN(objective_function<Type>* obj) {
  
//   // --- This is the crucial fix ---
//   // Safely save the original definition of the macro
//   #pragma push_macro("TMB_OBJECTIVE_PTR")
//   #undef TMB_OBJECTIVE_PTR
//   // Redefine the macro for THIS FUNCTION ONLY.
//   #define TMB_OBJECTIVE_PTR obj

//   // --- Data and Parameters ---
//   DATA_MATRIX(y); 
//   DATA_MATRIX(X); 
//   DATA_MATRIX(V); 
//   DATA_MATRIX(w); 

//   PARAMETER_VECTOR(alpha); 
//   PARAMETER_VECTOR(beta);  

//   Type nll = 0.0;

//   // --- Model Logic ---
//   vector<Type> log_lambda_j = X * beta;
//   vector<Type> lambda_j = exp(log_lambda_j);
//   vector<Type> lambda_tilde_i = w * lambda_j;
//   vector<Type> psi_i = Type(1.0) - exp(-lambda_tilde_i);
//   int M = y.rows();
//   int J = y.cols();
//   vector<Type> logit_p = V * alpha;

//   for (int i = 0; i < M; i++) {
//     Type log_prob_y_given_occupied = 0.0;
//     for (int t = 0; t < J; t++) {
//       Type p_it = Type(1.0) / (Type(1.0) + exp(-logit_p(i * J + t)));
//       log_prob_y_given_occupied += dbinom(y(i, t), Type(1.0), p_it, true);
//     }

//     if (y.row(i).sum() > 0) {
//       nll -= log(psi_i(i)) + log_prob_y_given_occupied;
//     } else {
//       Type prob_occupied_missed = psi_i(i) * exp(log_prob_y_given_occupied);
//       Type prob_unoccupied = Type(1.0) - psi_i(i);
//       nll -= log(prob_occupied_missed + prob_unoccupied);
//     }
//   }

//   // --- Restore the original macro definition ---
//   #pragma pop_macro("TMB_OBJECTIVE_PTR")

//   return -nll;
// }

// #endif // End of the include guard


#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type tmb_occuN(objective_function<Type>* obj) {
  // ++ DATA ++ //
  DATA_MATRIX(y);            // Detections
  DATA_MATRIX(X_state);      // Occupancy covs
  DATA_MATRIX(X_det);        // Detection covs
  DATA_MATRIX(W);            // Area weights

  // ++ PARAMETERS ++ //
  PARAMETER_VECTOR(beta_state);
  PARAMETER_VECTOR(beta_det);
  PARAMETER_VECTOR(log_lambda); // NNGP parameters

  // ++ NEGATIVE LOG-LIKELIHOOD ++ //
  Type nll = 0.0;

  // ++ PROCESS ++ //
  int M = y.rows(); // Number of sites
  int T = y.cols(); // Number of samples
  int J = W.cols(); // Number of raster cells

  // -- Detection model -- //
  matrix<Type> p(M, T);
  vector<Type> logit_p = X_det * beta_det;
  int p_counter = 0;
  for (int i = 0; i < M; i++){
    for (int j = 0; j < T; j++){
      p(i,j) = invlogit(logit_p(p_counter));
      p_counter++;
    }
  }

  // -- Occupancy model -- //
  vector<Type> lambda = exp(log_lambda);
  vector<Type> lambda_tilde = W * lambda;

  vector<Type> psi(M);
  for (int i = 0; i < M; i++){
    psi(i) = 1.0 - exp(-lambda_tilde(i));
  }

  // -- Likelihood -- //
  vector<Type> site_ndets(M);
  for (int i = 0; i < M; i++){
    site_ndets(i) = y.row(i).sum();
  }

  for (int i = 0; i < M; i++){

    // P(y_i|Z_i=1)
    Type log_lik_y_present = 0.0;
    for (int j = 0; j < T; j++) {
      if(y(i,j) == 0){
        log_lik_y_present += log(1.0 - p(i,j));
      } else {
        log_lik_y_present += log(p(i,j));
      }
    }

    Type psi_i = psi(i);

    // This is the corrected likelihood calculation
    Type log_psi_i = log(psi_i);

    if(site_ndets(i) == 0){
      // Site could be occupied with no detects, or unoccupied
      Type log_lik_present = log_psi_i + log_lik_y_present;
      Type log_lik_absent = log(1.0 - psi_i);
      // Use logspace_add for numerical stability: log(exp(a) + exp(b))
      nll -= logspace_add(log_lik_present, log_lik_absent);
    } else {
      // Site must be occupied. The nll is -(log(psi) + log(P(y|Z=1)))
      // The original code was effectively nll -= log(psi) + log(P(y|Z=1)), which was wrong.
      // Correct way is to add the negative log-likelihoods.
      nll -= log_psi_i;
      nll -= log_lik_y_present;
    }
  }


  // ++ REPORTING ++ //
  REPORT(psi);
  REPORT(p);
  REPORT(lambda);

  return nll;
}

#undef TMB_OBJECTIVE_PTR