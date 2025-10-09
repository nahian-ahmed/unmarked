# #' @export
# setClass("unmarkedFrameOccuN",
#          slots = c(w = "matrix",
#                    cellCovs = "data.frame"),
#          contains = "unmarkedFrameOccu")


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


# # Custom getDesign method for the unmarkedFrameOccuN class
# setMethod("getDesign", "unmarkedFrameOccuN",
#     function(umf, formula, na.rm = TRUE) {

#     M <- numSites(umf)
#     J <- obsNum(umf)

#     # --- Definitive Formula Parsing (This is the final fix) ---
#     # This is the standard, correct way to parse two-part formulas in unmarked
#     det_formula <- as.formula(formula[[2]])
#     state_formula <- as.formula(paste("~", formula[3], sep=""))

#     # --- Data Preparation ---
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
#     X <- designMats$X
#     V <- designMats$V
#     y <- designMats$y
#     w <- data@w

#     cat("SUCCESS! occuN function is running with the correct custom getDesign method.\n")
#     cat("State design matrix (X) has", nrow(X), "rows (cells).\n")
#     cat("Detection design matrix (V) has", nrow(V), "rows (sites x obs).\n")
# }


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

    # 1. Prepare data and parameters for TMB
    designMats <- getDesign(data, formula)
    
    # Add "model" to the data list to tell the dispatcher which model to run
    tmb_data <- list(model = "tmb_occuN",
                     y = designMats$y,
                     X = designMats$X,
                     V = designMats$V,
                     w = data@w)
    
    n_alpha <- ncol(designMats$V)
    n_beta <- ncol(designMats$X)

    if(missing(starts)) {
        starts <- rep(0, n_alpha + n_beta)
    }
    tmb_params <- list(alpha = starts[1:n_alpha],
                       beta = starts[(n_alpha + 1):(n_alpha + n_beta)])

    # 2. Call the TMB optimization engine
    # The DLL is now the main "unmarked" library, not a specific exports file.
    obj <- TMB::MakeADFun(data = tmb_data,
                          parameters = tmb_params,
                          DLL = "unmarked",
                          silent = TRUE)

    opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)

    # 3. Format and return the results (simplified for now)
    cat("SUCCESS! Model has been fit with TMB using the dispatcher.\n")
    print(opt)
    
}