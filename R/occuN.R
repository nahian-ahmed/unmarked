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

    # --- Definitive Data Preparation (Using the required 2-step process) ---

    # 1. Prepare data for the Detection Model
    sc <- umf@siteCovs
    if(nrow(sc) > 0) sc <- sc[rep(1:M, each = J), , drop = FALSE]
    oc <- as.data.frame(lapply(umf@obsCovs, as.vector))
    det_data <- cbind(sc, oc)
    
    # 2. Use model.frame() FIRST to create a clean model frame. This is the crucial step.
    mf_det <- model.frame(det_formula, det_data, na.action = na.pass)
    
    # 3. THEN use model.matrix() on the resulting model frame.
    V_design <- model.matrix(det_formula, mf_det)

    # 4. Repeat the 2-step process for the State Model
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
    X <- designMats$X
    V <- designMats$V
    y <- designMats$y
    w <- data@w

    cat("SUCCESS! occuN function is running with the correct custom getDesign method.\n")
    cat("State design matrix (X) has", nrow(X), "rows (cells).\n")
    cat("Detection design matrix (V) has", nrow(V), "rows (sites x obs).\n")
}