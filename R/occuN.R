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

    # --- Definitive Formula Parsing (This is the final fix) ---
    # This is the standard, correct way to parse two-part formulas in unmarked
    det_formula <- as.formula(formula[[2]])
    state_formula <- as.formula(paste("~", formula[3], sep=""))

    # --- Data Preparation ---
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
    X <- designMats$X
    V <- designMats$V
    y <- designMats$y
    w <- data@w

    cat("SUCCESS! occuN function is running with the correct custom getDesign method.\n")
    cat("State design matrix (X) has", nrow(X), "rows (cells).\n")
    cat("Detection design matrix (V) has", nrow(V), "rows (sites x obs).\n")
}