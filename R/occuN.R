
# # Load TMB, which is needed for the model fitting
# library(TMB)

# #' @export
# setClass("unmarkedFrameOccuN",
#          slots = c(w = "matrix",
#                    cellCovs = "data.frame"),
#          contains = "unmarkedFrameOccu")

# # Define the new output object for your model fit
# #' @export
# setClass("unmarkedFitOccuN",
#          contains = "unmarkedFitOccu")


# #' @export
# unmarkedFrameOccuN <- function(y, siteCovs = NULL, obsCovs = NULL,
#                                cellCovs, w, mapInfo = NULL) {

#     if(is.null(siteCovs)) {
#         siteCovs <- data.frame(site = 1:nrow(y))
#     }
#     parentFrame <- unmarkedFrameOccu(y = y, siteCovs = siteCovs,
#                                      obsCovs = obsCovs, mapInfo = mapInfo)
#     umf <- new("unmarkedFrameOccuN", parentFrame,
#                w = w,
#                cellCovs = cellCovs)
#     return(umf)
# }


# # Custom getDesign method
# setMethod("getDesign", "unmarkedFrameOccuN",
#     function(umf, formula, na.rm = TRUE) {
#     M <- numSites(umf)
#     J <- obsNum(umf)
#     det_formula <- as.formula(formula[[2]])
#     state_formula <- as.formula(paste("~", formula[3], sep=""))
#     sc <- umf@siteCovs
#     if(nrow(sc) > 0) sc <- sc[rep(1:M, each = J), , drop = FALSE]
#     oc <- as.data.frame(lapply(umf@obsCovs, as.vector))
#     det_data <- cbind(sc, oc)
#     mf_det <- model.frame(det_formula, det_data, na.action = na.pass)
#     V_design <- model.matrix(det_formula, mf_det)
#     mf_state <- model.frame(state_formula, umf@cellCovs, na.action = na.pass)
#     X_design <- model.matrix(state_formula, mf_state)
#     y <- getY(umf)
#     return(list(y = y, X = X_design, V = V_design))
# })


# #' @export
# occuN <- function(formula, data,
#                   starts, method = "BFGS", control = list(), se = TRUE) {

#     if(!is(data, "unmarkedFrameOccuN")) {
#         stop("Data is not an object of class unmarkedFrameOccuN.")
#     }

#     designMats <- getDesign(data, formula)
    
#     tmb_data <- list(model = "tmb_occuN",
#                      y = designMats$y, X = designMats$X,
#                      V = designMats$V, w = data@w)
    
#     n_alpha <- ncol(designMats$V)
#     n_beta <- ncol(designMats$X)
#     n_pars <- n_alpha + n_beta

#     if(missing(starts)) {
#         starts <- rep(0, n_pars)
#     }
#     tmb_params <- list(alpha = starts[1:n_alpha],
#                        beta = starts[(n_alpha + 1):n_pars])

#     obj <- TMB::MakeADFun(data = tmb_data, parameters = tmb_params,
#                           DLL = "unmarked_TMBExports", silent = TRUE)

#     opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
    
#     sd_rep <- TMB::sdreport(obj)
#     est_mat <- summary(sd_rep)
    
#     state_est <- unmarkedEstimate(name = "State", short.name = "lam",
#                                   estimates = est_mat[1:n_beta, 1],
#                                   covMat = sd_rep$cov.fixed[1:n_beta, 1:n_beta],
#                                   invlink = "exp", invlinkGrad = "exp")

#     det_est <- unmarkedEstimate(name = "Detection", short.name = "p",
#                                 estimates = est_mat[(n_beta + 1):n_pars, 1],
#                                 covMat = sd_rep$cov.fixed[(n_beta + 1):n_pars, (n_beta + 1):n_pars],
#                                 invlink = "logistic", invlinkGrad = "logistic.grad")

#     fit <- new("unmarkedFitOccuN",
#                fitType = "occuN",
#                call = match.call(),
#                formula = formula,
#                data = data,
#                sitesRemoved = numeric(0), # <-- FINAL FIX IS HERE
#                estimates = unmarkedEstimateList(list(state=state_est, det=det_est)),
#                AIC = 2 * opt$objective + 2 * n_pars,
#                opt = opt,
#                negLogLike = opt$objective,
#                nllFun = obj$fn)

#     return(fit)
# }

# setMethod("predict", "unmarkedFitOccuN",
#     function(object, type, newdata = NULL, backTransform = TRUE, ...) {

#     # Check for valid prediction type
#     if(!type %in% c("state", "det", "lambda")){
#       stop("Type must be 'state', 'det', or 'lambda'")
#     }

#     # If no newdata is provided, predict for the original data
#     if(is.null(newdata)){
#       newdata <- object@data
#     }

#     if(type == "state" || type == "lambda"){
#       # --- Predict State Process: Occupancy (psi) or Expected Abundance (lambda) ---

#       # 1. Get the state estimate object from the fitted model
#       state_est <- object@estimates['state']

#       # 2. Build the cell-level design matrix from the newdata's cellCovs
#       state_formula <- as.formula(paste("~", object@formula[3], sep=""))
#       X <- model.matrix(state_formula, newdata@cellCovs)

#       # 3. Calculate cell-level lambda on the log scale (linear predictor)
#       log_lambda_j_lc <- unmarked::linearComb(X, coef(state_est), ...)

#       # 4. Calculate site-level expected abundance (lambda_tilde)
#       # This requires matrix multiplication on the predictions and their variance
#       w <- newdata@w
#       lambda_tilde_lc <- unmarked::linearComb(w, log_lambda_j_lc, ...)

#       # 5. Handle back-transformation based on prediction type
#       if(type == "lambda"){
#         # Return site-level expected abundance
#         if(backTransform){
#           lambda_tilde_lc@estimate <- exp(lambda_tilde_lc@estimate)
#         }
#         return(lambda_tilde_lc)

#       } else { # type == "state"
#         # Return site-level occupancy probability (psi)
#         if(backTransform){
#           # Apply the cloglog-inverse transformation: psi = 1 - exp(-lambda_tilde)
#           lambda_tilde <- exp(lambda_tilde_lc@estimate)
#           psi <- 1 - exp(-lambda_tilde)
          
#           # Use the delta method to get SE on the probability scale
#           # The derivative of psi w.r.t lambda_tilde is exp(-lambda_tilde)
#           se <- sqrt( (exp(-lambda_tilde))^2 * diag(vcov(lambda_tilde_lc)) )
          
#           out <- data.frame(Predicted = psi, SE = se)
#           rownames(out) <- rownames(newdata@y)
#           return(out)
#         } else {
#           # On the linear scale, psi is equivalent to lambda_tilde
#           return(lambda_tilde_lc)
#         }
#       }
#     }

#     if(type == "det"){
#       # --- Predict Detection Process (p) ---
#       det_est <- object@estimates['det']
      
#       # Prepare detection data from newdata
#       M <- numSites(newdata)
#       J <- obsNum(newdata)
#       sc <- newdata@siteCovs
#       if(nrow(sc)>0) sc <- sc[rep(1:M, each=J),,drop=FALSE]
#       oc <- as.data.frame(lapply(newdata@obsCovs, as.vector))
#       det_data <- cbind(sc, oc)

#       # Build detection design matrix
#       det_formula <- as.formula(object@formula[[2]])
#       V <- model.matrix(det_formula, det_data)
      
#       # Get predictions
#       det_preds <- unmarked::linearComb(V, coef(det_est), ...)
      
#       if(backTransform){
#         det_preds@estimate <- plogis(det_preds@estimate)
#       }
#       return(det_preds)
#     }
# })
# Load TMB, which is needed for the model fitting
library(TMB)

#' @export
setClass("unmarkedFrameOccuN",
         slots = c(w = "matrix",
                   cellCovs = "data.frame"),
         contains = "unmarkedFrameOccu")

# Define the new output object for your model fit
#' @export
setClass("unmarkedFitOccuN",
         contains = "unmarkedFitOccu")


#' @export
unmarkedFrameOccuN <- function(y, siteCovs = NULL, obsCovs = NULL,
                               cellCovs, w, mapInfo = NULL) {

    if(is.null(siteCovs)) {
        siteCovs <- data.frame(site = 1:nrow(y))
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
    if(nrow(sc) > 0) sc <- sc[rep(1:M, each = J), , drop = FALSE]
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
                  starts, method = "BFGS", control = list(), se = TRUE) {

    if(!is(data, "unmarkedFrameOccuN")) {
        stop("Data is not an object of class unmarkedFrameOccuN.")
    }

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

    opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
    
    sd_rep <- TMB::sdreport(obj)
    est_mat <- summary(sd_rep, "fixed")
    
    state_est <- unmarkedEstimate(name = "State", short.name = "lam",
                                  estimates = est_mat[1:n_beta, 1],
                                  covMat = sd_rep$cov.fixed[1:n_beta, 1:n_beta],
                                  invlink = "exp", invlinkGrad = "exp")

    det_est <- unmarkedEstimate(name = "Detection", short.name = "p",
                                estimates = est_mat[(n_beta + 1):n_pars, 1],
                                covMat = sd_rep$cov.fixed[(n_beta + 1):n_pars, (n_beta + 1):n_pars],
                                invlink = "logistic", invlinkGrad = "logistic.grad")

    fit <- new("unmarkedFitOccuN",
               fitType = "occuN",
               call = match.call(),
               formula = formula,
               data = data,
               sitesRemoved = numeric(0),
               estimates = unmarkedEstimateList(list(state=state_est, det=det_est)),
               AIC = 2 * opt$objective + 2 * n_pars,
               opt = opt,
               negLogLike = opt$objective,
               nllFun = obj$fn)

    return(fit)
}


# --------------------------------------------------------------------------
# FINAL PREDICT METHOD WITH MULTI-LEVEL FUNCTIONALITY
# --------------------------------------------------------------------------

setMethod("predict", "unmarkedFitOccuN",
    function(object, type, newdata = NULL, backTransform = TRUE, ...) {

    valid_types <- c("state", "det", "lambda", "intensity")
    if(!type %in% valid_types){
      stop(paste("Type must be one of:", paste(valid_types, collapse=", ")))
    }

    if(is.null(newdata)){
      if(type == "intensity") stop("'newdata' is required for type='intensity'")
      newdata <- object@data
    }

    if(type %in% c("state", "lambda")){
      if(!is(newdata, "unmarkedFrameOccuN")) stop("'newdata' must be an unmarkedFrameOccuN for this prediction type")
      
      # Use '$' to correctly extract the unmarkedEstimate object
      state_est <- object@estimates$state
      log_lambda_j_lc <- linearComb(state_est, newdata = newdata@cellCovs)

      w <- newdata@w
      log_lambda_tilde_est <- w %*% log_lambda_j_lc@estimate
      vcov_log_lambda_tilde <- w %*% vcov(log_lambda_j_lc) %*% t(w)
      
      lambda_tilde_lc <- new("unmarkedLinComb", 
                             estimate = as.vector(log_lambda_tilde_est),
                             vcov = vcov_log_lambda_tilde)

      if(type == "lambda"){
        if(backTransform) lambda_tilde_lc@estimate <- exp(lambda_tilde_lc@estimate)
        return(lambda_tilde_lc)

      } else { # type == "state"
        if(backTransform){
          lambda_tilde <- exp(lambda_tilde_lc@estimate)
          psi <- 1 - exp(-lambda_tilde)
          se <- sqrt(diag( (exp(-lambda_tilde))^2 * vcov(lambda_tilde_lc) ))
          out <- data.frame(Predicted = psi, SE = se)
          rownames(out) <- 1:nrow(out)
          return(out)
        } else {
          return(lambda_tilde_lc)
        }
      }
    }

    if(type == "intensity"){
        if(!is(newdata, "data.frame")) stop("'newdata' must be a data.frame for type='intensity'")
        
        # Use '$' to correctly extract the unmarkedEstimate object
        state_est <- object@estimates$state
        preds <- linearComb(state_est, newdata = newdata)
        
        if(backTransform) preds@estimate <- exp(preds@estimate)
        return(preds)
    }

    if(type == "det"){
      if(!is(newdata, "unmarkedFrameOccuN")) stop("'newdata' must be an unmarkedFrameOccuN for this prediction type")

      # Use '$' to correctly extract the unmarkedEstimate object
      det_est <- object@estimates$det
      
      M <- numSites(newdata)
      J <- obsNum(newdata)
      sc <- newdata@siteCovs
      if(nrow(sc)>0) sc <- sc[rep(1:M, each=J),,drop=FALSE]
      oc <- as.data.frame(lapply(newdata@obsCovs, as.vector))
      det_data <- cbind(sc, oc)
      
      det_preds <- linearComb(det_est, newdata = det_data)
      
      if(backTransform) det_preds@estimate <- plogis(det_preds@estimate)
      return(det_preds)
    }
})