#' @export
setClass("unmarkedFrameOccuN",
         slots = c(w = "matrix",
                   cellCovs = "data.frame"),
         contains = "unmarkedFrameOccu")


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


# Custom getDesign method for the unmarkedFrameOccuN class
setMethod("getDesign", "unmarkedFrameOccuN",
    function(umf, formula, na.rm = TRUE) {

    M <- numSites(umf)
    J <- obsNum(umf)

    formula <- as.formula(formula)
    form_parts <- unmarked:::split_formula(formula)
    det_formula <- form_parts$det
    state_formula <- form_parts$state

    # --- Robustly Prepare Detection Data (Final Version) ---
    # 1. Create a base data frame with the correct number of rows for observations.
    det_data <- data.frame(matrix(NA, nrow = M * J, ncol = 0))

    # 2. Add site-level covariates, expanded to the observation level.
    if(ncol(umf@siteCovs) > 0) {
      site_covs_expanded <- umf@siteCovs[rep(1:M, each = J), , drop = FALSE]
      det_data <- cbind(det_data, site_covs_expanded)
    }

    # 3. Add observation-level covariates.
    if(length(umf@obsCovs) > 0) {
      obs_covs_df <- as.data.frame(lapply(umf@obsCovs, as.vector))
      det_data <- cbind(det_data, obs_covs_df)
    }
    
    rownames(det_data) <- NULL

    # Build the design matrices
    V_design <- model.matrix(det_formula, det_data)
    X_design <- model.matrix(state_formula, umf@cellCovs)

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
    X <- designMats$X
    V <- designMats$V
    y <- designMats$y
    w <- data@w

    cat("SUCCESS! occuN function is running with the correct custom getDesign method.\n")
    cat("State design matrix (X) has", nrow(X), "rows (cells).\n")
    cat("Detection design matrix (V) has", nrow(V), "rows (sites x obs).\n")
}