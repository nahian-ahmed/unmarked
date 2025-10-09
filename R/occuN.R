# /R/occuN.R

#' @export
setClass("unmarkedFrameOccuN",
         slots = c(w = "matrix"),
         contains = "unmarkedFrameOccu")


#' @export
unmarkedFrameOccuN <- function(y, siteCovs=NULL, obsCovs=NULL, w, mapInfo) {

    # First, call the parent constructor to create the base object.
    # This handles all the complex validation of y, siteCovs, etc.
    parentFrame <- unmarkedFrameOccu(y = y, siteCovs = siteCovs,
                                     obsCovs = obsCovs, mapInfo = mapInfo)

    # Now, create the new unmarkedFrameOccuN object.
    # We pass the parent object to new(), which copies all the inherited slots,
    # and then we provide the value for our new 'w' slot.
    umf <- new("unmarkedFrameOccuN", parentFrame, w = w)

    return(umf)
}


#' @export
occuN <- function(formula, data,
                  starts, method = "BFGS", control = list(), se = TRUE) {

    if(!is(data, "unmarkedFrameOccuN")) {
        stop("Data is not an object of class unmarkedFrameOccuN.")
    }

    formula <- as.formula(formula)
    designMats <- getDesign(data, formula)
    X <- designMats$X
    V <- designMats$V
    y <- designMats$y
    w <- data@w

    cat("occuN function skeleton is running!\n")
    cat("Successfully extracted design matrices and weights matrix.\n")
}