# A Latent Abundance Framework for Modeling Occupancy at Sites of Arbitrary Geometry

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
        siteCovs <- data.frame(site = seq_len(nrow(y)))
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



    if (method == "nlminb") {
        opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)
    } else if (method %in% c("Nelder-Mead", "SANN")) {
        # Methods that do not use gradients
        opt <- optim(obj$par, obj$fn, method = method, control = control)
    } else {
        # Methods that do use gradients (e.g., "BFGS", "L-BFGS-B", "CG")
        opt <- optim(obj$par, obj$fn, obj$gr, method = method, control = control)
    }

    # Handle slightly different output formats from nlminb and optim
    # The negative log-likelihood is in 'objective' for nlminb and 'value' for optim
    nll <- if (method == "nlminb") opt$objective else opt$value
    
    sd_rep <- TMB::sdreport(obj)
    est_mat <- summary(sd_rep, "fixed")


    # 'state' (beta) is second in the parameter list
    state_est <- unmarkedEstimate(name = "State", short.name = "lam",
                                  estimates = est_mat[(n_alpha + 1):n_pars, 1],
                                  covMat = sd_rep$cov.fixed[(n_alpha + 1):n_pars, (n_alpha + 1):n_pars],
                                  invlink = "exp", invlinkGrad = "exp")

    # 'det' (alpha) is first in the parameter list
    det_est <- unmarkedEstimate(name = "Detection", short.name = "p",
                                estimates = est_mat[1:n_alpha, 1],
                                covMat = sd_rep$cov.fixed[1:n_alpha, 1:n_alpha],
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