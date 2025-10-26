#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type tmb_occuN(objective_function<Type>* obj) {

  // Data and Parameters
  DATA_MATRIX(y);
  DATA_MATRIX(X);
  DATA_MATRIX(V);
  DATA_MATRIX(w);

  DATA_SCALAR(lambda_reg_beta);
  DATA_SCALAR(lambda_reg_alpha);

  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta);

  PARAMETER(alpha_lambda);

  Type nll = 0.0;

  // Model Logic
  vector<Type> log_lambda_j = X * beta;
  vector<Type> lambda_j = exp(log_lambda_j);
  vector<Type> lambda_tilde_i = w * lambda_j;
  vector<Type> psi_i = Type(1.0) - exp(-lambda_tilde_i);
  int M = y.rows();
  int J = y.cols();
  vector<Type> logit_p_base = V * alpha;

  for (int i = 0; i < M; i++) {
    Type log_prob_y_given_occupied = 0.0;
    Type log_lambda_i = log(lambda_tilde_i(i));
    for (int t = 0; t < J; t++) {
      if(y(i,t) == y(i,t)) {
        Type logit_p_final_it = logit_p_base(i * J + t) + alpha_lambda * log_lambda_i;
        Type p_it = Type(1.0) / (Type(1.0) + exp(-logit_p_final_it));
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


  if(beta.size() > 1) {
    nll += lambda_reg_beta * (beta.segment(1, beta.size() - 1).square().sum());
  }


  Type det_penalty_sumsq = 0.0;
  
  if(alpha.size() > 1) {
    det_penalty_sumsq += alpha.segment(1, alpha.size() - 1).square().sum();
  }
  
  det_penalty_sumsq += pow(alpha_lambda, 2.0);
  
  nll += lambda_reg_alpha * det_penalty_sumsq;


  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this