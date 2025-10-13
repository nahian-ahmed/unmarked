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
      if(y(i,t) == y(i,t)) {
        Type p_it = Type(1.0) / (Type(1.0) + exp(-logit_p(i * J + t)));
        log_prob_y_given_occupied += dbinom(y(i, t), Type(1.0), p_it, true);
      }
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

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this