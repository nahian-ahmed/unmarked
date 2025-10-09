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

    # --- Definitive Data Preparation for Detection Model ---
    # 1. Expand site-level covariates to the observation level
    sc <- umf@siteCovs
    if(nrow(sc) > 0) {
      sc <- sc[rep(1:M, each = J), , drop = FALSE]
    }

    # 2. Convert the observation-level covariates list to a data frame
    oc <- as.data.frame(lapply(umf@obsCovs, as.vector))

    # 3. Combine them into a single, valid data frame
    det_data <- cbind(sc, oc)

    # 4. Use this clean data frame to build the model matrix
    V_design <- model.matrix(det_formula, det_data)

    # --- Prepare State Model Data ---
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