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
  // Names now match the data list sent from R/occuN.R
  DATA_MATRIX(y);
  DATA_MATRIX(X_lambda);
  DATA_MATRIX(X_p);
  DATA_SPARSE_MATRIX(w_mat);

  // ++ PARAMETERS ++ //
  // Names now match the parameters list sent from R/occuN.R
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(alpha);

  // ++ PROCESS ++ //
  int M = y.rows(); // Number of sites
  int T = y.cols(); // Number of observations
  
  // -- Detection model -- //
  matrix<Type> p(M, T);
  vector<Type> logit_p_vec = X_p * alpha;
  for (int i = 0; i < M; i++){
    for (int j = 0; j < T; j++){
      p(i,j) = invlogit(logit_p_vec(i*T + j));
    }
  }

  // -- Occupancy model -- //
  // CORRECT: Calculate log_lambda from covariates (X_lambda) and coefficients (beta)
  vector<Type> log_lambda = X_lambda * beta;
  vector<Type> lambda = exp(log_lambda);
  
  // Re-define sparse matrix type for iterator
  typedef Eigen::SparseMatrix<Type> spmat_t;
  vector<Type> lambda_tilde(M);
  lambda_tilde.setZero();
  for (int i = 0; i < M; i++) {
    for (typename spmat_t::InnerIterator it(w_mat, i); it; ++it) {
      lambda_tilde(i) += it.value() * lambda(it.col());
    }
  }

  vector<Type> psi = 1.0 - exp(-lambda_tilde);

  // -- NEGATIVE LOG-LIKELIHOOD -- //
  Type nll = 0.0;
  
  vector<Type> site_ndets(M);
  for (int i = 0; i < M; i++){
    site_ndets(i) = y.row(i).sum();
  }

  for (int i = 0; i < M; i++){
    Type log_lik_y_present = 0.0;
    for (int j = 0; j < T; j++) {
      log_lik_y_present += dbinom(y(i,j), Type(1.0), p(i,j), true);
    }

    if(site_ndets(i) == 0){
      nll -= log(psi(i) * exp(log_lik_y_present) + (1.0 - psi(i)));
    } else {
      nll -= log(psi(i)) + log_lik_y_present;
    }
  }

  // ++ REPORTING ++ //
  REPORT(psi);
  REPORT(p);
  REPORT(lambda);

  return nll;
}

#undef TMB_OBJECTIVE_PTR