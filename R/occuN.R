# /R/occuN.R

#' @export
# Add a new slot 'pixelCovs' to hold the pixel-level data
setClass("unmarkedFrameOccuN",
         slots = c(w = "matrix",
                   pixelCovs = "data.frame"),
         contains = "unmarkedFrameOccu")


#' @export
# The constructor now takes 'pixelCovs' as a distinct argument
unmarkedFrameOccuN <- function(y, siteCovs = NULL, obsCovs = NULL,
                               pixelCovs, w, mapInfo = NULL) {

    # 1. Create a placeholder siteCovs if one isn't provided.
    #    This ensures the parent object always has a valid siteCovs table.
    if(is.null(siteCovs)) {
        siteCovs <- data.frame(site = 1:nrow(y))
    }

    # 2. Call the parent constructor with the SITE-LEVEL data.
    parentFrame <- unmarkedFrameOccu(y = y, siteCovs = siteCovs,
                                     obsCovs = obsCovs, mapInfo = mapInfo)

    # 3. Create the final object, passing the parent and the NEW data.
    umf <- new("unmarkedFrameOccuN", parentFrame,
               w = w,
               pixelCovs = pixelCovs)

    return(umf)
}


#' @export
# The main function is unchanged for now
occuN <- function(formula, data,
                  starts, method = "BFGS", control = list(), se = TRUE) {

    if(!is(data, "unmarkedFrameOccuN")) {
        stop("Data is not an object of class unmarkedFrameOccuN.")
    }

    # NOTE: The getDesign function will eventually need to be modified to
    #       use data@pixelCovs for the state formula.
    #       For now, the skeleton is fine.

    formula <- as.formula(formula)
    designMats <- getDesign(data, formula)
    X <- designMats$X
    V <- designMats$V
    y <- designMats$y
    w <- data@w

    cat("occuN function skeleton is running!\n")
    cat("Successfully extracted design matrices and weights matrix.\n")
}
