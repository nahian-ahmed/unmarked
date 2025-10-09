# /R/occuN.R

#' @export
# Defines the S4 class, now with the 'cellCovs' slot
setClass("unmarkedFrameOccuN",
         slots = c(w = "matrix",
                   cellCovs = "data.frame"),
         contains = "unmarkedFrameOccu")


#' @export
# The constructor function, now taking 'cellCovs' as an argument
unmarkedFrameOccuN <- function(y, siteCovs = NULL, obsCovs = NULL,
                               cellCovs, w, mapInfo = NULL) {

    # Create a placeholder for siteCovs if none is provided
    if(is.null(siteCovs)) {
        siteCovs <- data.frame(site = seq_len(nrow(y)))
    }

    # Call the parent constructor with site-level data
    parentFrame <- unmarkedFrameOccu(y = y, siteCovs = siteCovs,
                                     obsCovs = obsCovs, mapInfo = mapInfo)

    # Create the final object with the new 'w' and 'cellCovs' data
    umf <- new("unmarkedFrameOccuN", parentFrame,
               w = w,
               cellCovs = cellCovs)

    return(umf)
}


# Custom getDesign method for the unmarkedFrameOccuN class
setMethod("getDesign", "unmarkedFrameOccuN",
    function(umf, formula, na.rm = TRUE) {

    # Separate the formula into state (abundance) and detection parts
    formula <- as.formula(formula)
    form_parts <- unmarked:::split_formula(formula) # Use unmarked internal function
    det_formula <- form_parts$det
    state_formula <- form_parts$state

    # Build the detection design matrix (V) from obsCovs and siteCovs
    # This allows for formulas like ~ wind_speed + site_access
    det_data <- cbind(umf@obsCovs, umf@siteCovs)
    V_design <- model.matrix(det_formula, det_data)

    # Build the state design matrix (X) from our new 'cellCovs' slot
    # This is the core of the occuN model's functionality
    X_design <- model.matrix(state_formula, umf@cellCovs)

    y <- getY(umf)

    return(list(y = y, X = X_design, V = V_design))
})


#' @export
# The main user-facing function for fitting the occuN model
occuN <- function(formula, data,
                  starts, method = "BFGS", control = list(), se = TRUE) {

    # Check for the correct data object type
    if(!is(data, "unmarkedFrameOccuN")) {
        stop("Data is not an object of class unmarkedFrameOccuN.")
    }

    # This call now dispatches to our custom getDesign method
    designMats <- getDesign(data, formula)
    X <- designMats$X
    V <- designMats$V
    y <- designMats$y
    w <- data@w

    cat("SUCCESS! occuN function is running with the correct custom getDesign method.\n")
    cat("State design matrix (X) has", nrow(X), "rows (cells).\n")
    cat("Detection design matrix (V) has", nrow(V), "rows (sites x obs).\n")
}
