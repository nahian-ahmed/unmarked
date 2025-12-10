# A Latent Abundance Framework for Modeling Occupancy at Sites of Arbitrary Geometry

#' @export
setClass("unmarkedFrameOccuN",
     slots = c(w = "dgCMatrix", 
           cellCovs = "data.frame"),
     contains = "unmarkedFrameOccu")



#' @export
setClass("unmarkedFitOccuN",
     contains = "unmarkedFitOccu")


#' @export
unmarkedFrameOccuN <- function(y, siteCovs = NULL, obsCovs = NULL,
                 cellCovs, w, mapInfo = NULL) {

  if(is.null(siteCovs)) {
    siteCovs <- data.frame(site = seq_len(nrow(y)))
  }
  
  if(!inherits(w, "dgCMatrix")) {
    w <- methods::as(w, "dgCMatrix")
  }

  
  parentFrame <- unmarkedFrameOccu(y = y, siteCovs = siteCovs,
                   obsCovs = obsCovs, mapInfo = mapInfo)
  umf <- new("unmarkedFrameOccuN", parentFrame,
         w = w,
         cellCovs = cellCovs)
  return(umf)
}


# Custom getDesign method
setMethod("getDesign", "unmarkedFrameOccuN",
  function(umf, formula, na.rm = TRUE) {
    
    M <- numSites(umf)
    J <- obsNum(umf)
    
    det_formula <- as.formula(formula[[2]])
    state_formula <- as.formula(paste("~", formula[3], sep=""))
    sc <- umf@siteCovs
    
    if(nrow(sc) > 0) 
      sc <- sc[rep(1:M, each = J), , drop = FALSE]
    
    oc <- as.data.frame(lapply(umf@obsCovs, as.vector))
    det_data <- cbind(sc, oc)
    
    mf_det <- model.frame(det_formula, det_data, na.action = na.pass)
    V_design <- model.matrix(det_formula, mf_det)
    mf_state <- model.frame(state_formula, umf@cellCovs, na.action = na.pass)
    X_design <- model.matrix(state_formula, mf_state)
    
    y <- getY(umf)
    
    return(list(y = y, X = X_design, V = V_design))
})



#' @export
occuN <- function(formula, data,
                  starts, method = "BFGS", control = list(), se = TRUE,
                  lower = -Inf, upper = Inf) {

  if(!is(data, "unmarkedFrameOccuN")) {
    stop("Data is not an object of class unmarkedFrameOccuN.")
  }

  print("NEW ONE WITH SE")
  designMats <- getDesign(data, formula)
  
  tmb_data <- list(model = "tmb_occuN",
                   y = designMats$y, X = designMats$X,
                   V = designMats$V, w = data@w)

  n_alpha <- ncol(designMats$V)
  n_beta <- ncol(designMats$X)
  n_pars <- n_alpha + n_beta

  if(missing(starts)) {
    starts <- rep(0, n_pars)
  }

  tmb_params <- list(alpha = starts[1:n_alpha],
                     beta = starts[(n_alpha + 1):n_pars])

  obj <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params,
                        DLL = "unmarked_TMBExports", silent = TRUE)

  # Optimization with bounds
  if (method == "nlminb") {
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = control, 
                  lower = lower, upper = upper)
  } else if (method %in% c("Nelder-Mead", "SANN")) {
    opt <- optim(obj$par, obj$fn, method = method, control = control, 
                 lower = lower, upper = upper)
  } else {
    opt <- optim(obj$par, obj$fn, obj$gr, method = method, control = control, 
                 lower = lower, upper = upper)
  }

  # extract nll and parameter estimates
  nll <- if (method == "nlminb") opt$objective else opt$value
  ests <- opt$par
  names(ests) <- names(obj$par)

  # Calculate SE only if requested
  if(se) {
    if (method != "nlminb") {
      obj$fn(opt$par) 
    }
    
    sd_rep <- TMB::sdreport(obj)
    est_mat <- summary(sd_rep, "fixed")
    
    ests_alpha   <- est_mat[1:n_alpha, 1]
    se_alpha     <- est_mat[1:n_alpha, 2]
    covMat_alpha <- sd_rep$cov.fixed[1:n_alpha, 1:n_alpha]
    
    ests_beta    <- est_mat[(n_alpha + 1):n_pars, 1]
    se_beta      <- est_mat[(n_alpha + 1):n_pars, 2]
    covMat_beta  <- sd_rep$cov.fixed[(n_alpha + 1):n_pars, (n_alpha + 1):n_pars]
    
  } else {
    # If SE is FALSE, fill with NAs to save time
    ests_alpha   <- ests[1:n_alpha]
    se_alpha     <- rep(NA, n_alpha)
    covMat_alpha <- matrix(NA, n_alpha, n_alpha)
    
    ests_beta    <- ests[(n_alpha + 1):n_pars]
    se_beta      <- rep(NA, n_beta)
    covMat_beta  <- matrix(NA, n_beta, n_beta)
  }

  # Create estimate objects
  # 'state' (beta) is second in the parameter list
  state_est <- unmarkedEstimate(name = "State", short.name = "lam",
                                estimates = ests_beta,
                                covMat = covMat_beta,
                                invlink = "exp", invlinkGrad = "exp")

  # 'det' (alpha) is first in the parameter list
  det_est <- unmarkedEstimate(name = "Detection", short.name = "p",
                              estimates = ests_alpha,
                              covMat = covMat_alpha,
                              invlink = "logistic", invlinkGrad = "logistic.grad")


  fit <- new("unmarkedFitOccuN",
             fitType = "occuN",
             call = match.call(),
             formula = formula,
             data = data,
             sitesRemoved = numeric(0),
             estimates = unmarkedEstimateList(list(state=state_est, det=det_est)),
             AIC = 2 * nll + 2 * n_pars,
             opt = opt,
             negLogLike = nll,
             nllFun = obj$fn)

  return(fit)
}