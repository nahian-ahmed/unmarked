# /R/occuN.R

#' @export
setClass("unmarkedFrameOccuN",
         slots = c(w = "matrix"),
         contains = "unmarkedFrameOccu")

#' @export
unmarkedFrameOccuN <- function(y, siteCovs=NULL, obsCovs=NULL, w, mapInfo) {

    # (We will fill in the validation and processing code later)

    umf <- new("unmarkedFrameOccuN", y=y, siteCovs=siteCovs,
               obsCovs=obsCovs, w=w, mapInfo=mapInfo,
               numPrimary=1)

    return(umf)
}

#' @export
occuN <- function(formula, data,
                  starts, method = "BFGS", control = list(), se = TRUE) {

    # 1. Check for the correct data object type
    if(!is(data, "unmarkedFrameOccuN")) {
        stop("Data is not an object of class unmarkedFrameOccuN.")
    }

    # 2. Process the formulas (e.g., ~detform ~stateform)
    # This extracts the formulas for the detection and state (abundance) processes
    formula <- as.formula(formula)
    designMats <- getDesign(data, formula)
    X <- designMats$X # State covariates
    V <- designMats$V # Detection covariates
    y <- designMats$y # Response variable

    # 3. Extract the weights matrix from our custom unmarkedFrame
    w <- data@w

    # (Future steps will go here: prepare data for C++, call the optimizer)

    # For now, let's just print a success message
    cat("occuN function skeleton is running!\n")
    cat("Successfully extracted design matrices and weights matrix.\n")
}
