
## @params is a vector of paramters that can be optimized
## @paramType is a vector of labels for the parameters that can be used to access them by name (same order as params)
## @fixed a vector of logical values that dertermines if a parameter is to be optimized (same order as params)
## @lowerBound is a vector of values that sets a lower bound on the parameter for box optimization (same order as params)
## @upperBound is a vector of values that sets an upper bound on the parameter for box optimization (same order as params)
## @invariants is a list to hold any parameters that are constant by the nature of the problem and can never be optimized. For example
##             the enumeration of the tree states, or scaling values for the data matricies. The list elements can take any form for
##             ease of use.
## @chainSpecificParameters is a two layer list where the first matches a parameter name and the second groups together chains where the
##                          value of that parameter is tied


parameterObject <- R6Class("parameterObject",public=list(params=numeric(),paramType=character(),paramIndex=data.table(),fixed=logical(),
                                   lowerBound=numeric(),upperBound=numeric(),invariants=list()))

parameterObject$set("public","initialize",function(data,params=list(),lowerBound=list(), upperBound=list(),fixed=list(),
                                              invariants=list(),chainSpecificParameters=NULL){
    stop("parameterObject is an abstract class.")
})

#####
## Methods for checking constraints
#####

## Checks that parameter constraints are specified in a valid fashion
parameterObject$set("public","checkParamConstraints", function(params,lowerBound=list(),upperBound=list(),fixed=list()){
    ## For each parameter check that it's constraint is supplied in a viable fashion
    ## Must be either the same length as the supplied parameter vector, only one element, or null
    for(p in names(params)){
        ## browser()
        if(length(lowerBound[[p]])!=length(params[[p]]) & length(lowerBound[[p]])!=1 & !is.null(lowerBound[[p]])){
            stop(paste("Lower bound constraint of parameter",p,"is specified in an invalid fashion"))
        }
        if(length(upperBound[[p]])!=length(params[[p]]) & length(upperBound[[p]])!=1 & !is.null(upperBound[[p]])){
            stop(paste("Upper bound constraint of parameter",p,"is specified in an invalid fashion"))
        }
        if(length(fixed[[p]])!=length(params[[p]]) & length(fixed[[p]])!=1 & !is.null(fixed[[p]])){
            stop(paste("Fixed constraint of parameter",p,"is specified in an invalid fashion"))
        }
    }
    ## Make sure there are no constraints specified for non-existant parameters
    if(sum(!names(lowerBound) %in% names(params))>0){
        stop(paste("A lower bound constraint was specified for the following non-existant parameters:",names(lowerBound)[!names(lowerBound) %in% names(params)]))
    }
    if(sum(!names(upperBound) %in% names(params))>0){
        stop(paste("A upper bound constraint was specified for the following non-existant parameters:",names(upperBound)[!names(upperBound) %in% names(params)]))
    }
    if(sum(!names(fixed) %in% names(params))>0){
        stop(paste("A fixed constraint was specified for the following non-existant parameters:",names(fixed)[!names(fixed) %in% names(params)]))
    }
})

## A function that checks the validity of chainSpecificParameters
parameterObject$set("public","checkChainSpecificParameters", function(params,chainSpecificParameters){
    ## Make sure there is no parameter tying specified for non-existant parameters
    if(sum(!names(chainSpecificParameters) %in% names(params))>0){
        stop(paste("A chainSpecificParameter was specified for the following non-existant parameters:",
                   names(chainSpecificParameters)[!names(chainSpecificParameters) %in% names(params)]))
    }
    for(p in names(chainSpecificParameters)){
        ## Make sure that a chain is not in more than one group
        if(sum(duplicated(unlist(chainSpecificParameters[[p]])))>0){
            stop(paste("Parameter",p,"has a chain that appears in more than one group of chain specific parameters."))
        }
        ## Make sure that all chains are in at least one group
        if(sum(! 1:length(self$data) %in% unlist(chainSpecificParameters[[p]]) )>0){
            stop(paste("Parameter",p," has chains that have not been assigned to a group in chainSpecificParameters"))
        }
        ## Make sure that the list does not contain an out of range chain
        if(sum(unlist(chainSpecificParameters[[p]])>length(self$data))>0){
            stop(paste("Parameter",p,"is tied to chain that does not exist. Please verfy indicies of chainSpecificParameters"))
        }
    }
})

###
## Setter functions
###

## A function that takes the list of parameters a converts it to a vector which is used to set
## self$params and builds lookup table
## NOTE: self$data must be set for this to be used
parameterObject$set("public","setParamVector", function(params,chainSpecificParameters){
    ## Do validity check on chain specific parameters if they exist
    if(length(chainSpecificParameters)>0){
        self$checkChainSpecificParameters(params,chainSpecificParameters)
    }
    ## Iterate over all parameters
    ind=1
    pList=list()
    for(p in names(params)){
        ## If the parameter is a chain specific parameter, create multiple copies of it and
        ## label with a .chain.indexNum suffix
        if(p %in% names(chainSpecificParameters)){
            self$params=c(self$params,rep(params[[p]],length(chainSpecificParameters[[p]])))
            ## Build index for parameters
            start=seq(ind,ind+length(params[[p]])*length(chainSpecificParameters[[p]]),by=length(params[[p]]))
            end=start+length(params[[p]])-1
            pList[[p]]=rbindlist(lapply(as.list(1:length(chainSpecificParameters[[p]])),
                     function(y){
                         data.table(chains=chainSpecificParameters[[p]][[y]],start=start[y],end=end[y])
                     }),idcol="group")
            ind=ind+length(params[[p]])*length(chainSpecificParameters[[p]])
        } else {
            self$params=c(self$params,params[[p]])
            ## Use -1 to denote non-specific parameters
            pList[[p]]=data.table(chains=-1,group=-1,start=ind,end=ind+length(params[[p]])-1)
            ind=ind+length(params[[p]])
        }
    }
    self$paramIndex=rbindlist(pList,idcol="paramType")
    setkeyv(self$paramIndex,cols=c("paramType","chains"))
})

## A function that takes the list of parameter constraints and uses it to set constraint variables in self
parameterObject$set("public","setParamConstraints", function(params,lowerBound=list(),upperBound=list(),fixed=list()){
    self$lowerBound=numeric(length(self$params))
    self$upperBound=numeric(length(self$params))
    self$fixed=numeric(length(self$params))
    ## Iterate over all parameters
    for(p in names(params)){
        p.ind=self$getParamIndicies(p,NULL)
        if(is.null(lowerBound[[p]])){
            lowerBound[[p]]=-Inf
        }
        if(is.null(upperBound[[p]])){
            upperBound[[p]]=Inf
        }
        if(is.null(fixed[[p]])){
            fixed[[p]]=FALSE
        }
        self$lowerBound[p.ind]=rep(lowerBound[[p]],length.out=length(p.ind))
        self$upperBound[p.ind]=rep(upperBound[[p]],length.out=length(p.ind))
        self$fixed[p.ind]=rep(fixed[[p]],length.out=length(p.ind))
    }
    self$fixed=as.logical(self$fixed)
})

####
## Parameter getter function
####
## Returns the indicies that a parameter exists in, in the parameter vector
## Note: Relies on the construction of the parameter vector where parameters of the same type and chain are contiguous
parameterObject$set("public","getParamIndicies", function(paramType,chain=NULL){
    ## If a chain was specified (but may be ignored)
    if(!is.null(chain)){
        ## If it is a chain linked variable
        if(nrow(self$paramIndex[paramType])>1){
            out=as.numeric(self$paramIndex[eval(.(paramType,chain)),.(start,end)])
            out=out[1]:out[2]
            if(is.na(out[1])){
                stop(paste0("Invalid chain in parameter query (",paramType,",",chain,")"))
            }
        } else {
            ## If it's not a chain linked variable
            if(!is.na(self$paramIndex[paramType]$group)){
                out=as.numeric(self$paramIndex[paramType,.(start,end)])
                out=out[1]:out[2]
            } else {
                stop(paste("Invalid query for parameter of type",paramType))
            }
        }
    } else {
        ## If the query was just for all indicies of a given parameter type
        if(!is.na(self$paramIndex[paramType]$group[1])){
            out=unlist(apply(self$paramIndex[paramType,.(start,end)],1, function(x) x[1]:x[2]))
        } else {
            stop(paste("Invalid query for parameter of type",paramType))
        }
    }
    return(as.numeric(out))
})

###
## Convenience functions
###

## Returns a table of parameters for the HMM object
parameterObject$set("public","getParameterTable",function(){
    pTab=foreach(i=1:nrow(self$paramIndex)) %do% {
        with(self$paramIndex,data.table(paramName=rep(paramType[i],end[i]-start[i]+1),chain=rep(chains[i],end[i]-start[i]+1),
                                        value=self$params[self$paramIndex$start[i]:self$paramIndex$end[i]]))
    }
    pTab=rbindlist(pTab)
    return(pTab)
})

####
## Functions that may be instantiated in a subclass
####

## Checks that parameter constraints are satisfied, for superclass emission nothing is checked
parameterObject$set("public","checkParamValidity", function(params){
    warning("No validity check for parameter implemented.")
})

## Checks that invariant constraints are satisfied, for superclass emission nothing is checked
parameterObject$set("public","checkInvariantValidity", function(invariants){
        warning("No validity check for invariants implemented.")
})
