## Set generic transition class
Transition <- R6::R6Class("Transition",public=list(transitionLogProb=matrix(),nstates=numeric(),hmm=NULL),inherit=parameterObject)

Transition$set("public","initialize",function(nstates,params=list(),lowerBound=list(), upperBound=list(),fixed=list(),
                                              invariants=list(),chainSpecificParameters=NULL){
    self$checkParamConstraints(params,lowerBound,upperBound,fixed)
    ## Set number of states
    self$nstates=nstates
    ## Invariant check with no dependencies
    if(is.list(invariants)){
        self$invariants=invariants
    } else stop("Invariants must be in a list.")
    ## Perform parameter validity check and set parameters
    self$checkParamValidity(params)

    ## Build parameter index
    self$buildParamIndex(params,chainSpecificParameters)
    
    ## Set parameter vector
    self$setParamVector(params)

    ## Reshape parameter constraints to vectors and set self$lowerBound, self$upperBound, self$fixed, must be after params are set
    self$setParamConstraints(lowerBound,upperBound,fixed)

    ## Now perform additional checks for invariant validity
    self$checkInvariantValidity(invariants)
    self$checkTransitionValidity()
})

Transition$set("public","checkTransitionValidity", function(){
    errors=character()
    if(sum(!is.numeric(self$params) | is.na(self$params))>0){
        errors=c(errors,"Non-numeric value in parameter list.")
    }
    if (length(errors) == 0) TRUE else errors
})

## Updates log transition probabilities
Transition$set("public","updateTransitionProbabilities", function(){
        stop("Must implement a updateEmissionProbabilities function.")
})

## An optional function that allows certain transition updates to always happen, regardless of whether any transition parameters have been updated
## Really only do this if you have a really weird HMM model where the transition and emission models share some parameters so the HMM
## class cannot correctly track when to update the appropriate probabilities
Transition$set("public","forcedTransitionUpdates", function(){
    ## warning("No forcedTransitionUpdates function set")
    return(NULL)
})
