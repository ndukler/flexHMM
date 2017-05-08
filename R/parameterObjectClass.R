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
                                   lowerBound=numeric(),upperBound=numeric(),invariants=list(),replicatesPerParameter=numeric(0)))

parameterObject$set("public","initialize",function(data,params=list(),lowerBound=list(), upperBound=list(),fixed=list(),
                                              invariants=list(),chainSpecificParameters=NULL,replicatesPerParameter=numeric(0)){
    stop("parameterObject is an abstract class.")
})

#####
## Methods for checking settings
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

## A function that checks the validity of chainSpecificParameters
parameterObject$set("public","checkReplicates", function(params,replicatesPerParameter){
    ## Make sure there is no parameter tying specified for non-existant parameters
    if(sum(!names(replicatesPerParameter) %in% names(params))>0){
        stop(paste("A chainSpecificParameter was specified for the following non-existant parameters:",
                   names(replicatesPerParameter)[!names(replicatesPerParameter) %in% names(params)]))
    }
})

###
## Index building
###

## Function to build parameter index
parameterObject$set("public","buildParamIndex", function(params,chainSpecificParameters=list(),replicatesPerParameter=character()){
    ## Do validity check on chain specific parameters if they exist
    if(length(chainSpecificParameters)>0){
        self$checkChainSpecificParameters(params,chainSpecificParameters)
    }
    ## Do validity checks on replicates if they exist
    if(length(replicatesPerParameter)>0){
        self$checkReplicates(params,replicatesPerParameter)
    }
    ## Iterate over all parameters
    pIndList=list()
    ind=1
    for(p in names(params)){
        ## Get the length of each parameter
        p.len=length(params[[p]])
        ## If the parameter is replicate specific get number of replicates, else 1
        if(!is.na(replicatesPerParameter[p])){
            rpls=1:replicatesPerParameter[p]
        } else{
            rpls=-1
        }
        ## If the parameter is chain specific get the number of chains
        if(!is.null(chainSpecificParameters[[p]])){
            pIndList[[p]]=rbindlist(lapply(chainSpecificParameters[[p]], function(x) data.table(chain=rep(x,each=length(rpls)),replicate=rpls)),idcol="group")
        } else{
            pIndList[[p]]=data.table(group=-1,chain=-1,replicate=rpls)
        }
        pIndList[[p]][,tmp:=as.numeric(factor(paste0(group,".",replicate)))-1]
        ## Set start and end indicies
        pIndList[[p]][,start:=ind+tmp*p.len]
        pIndList[[p]][,end:=start+p.len-1]
        ind=max(pIndList[[p]]$end)+1
        pIndList[[p]][,tmp:=NULL]
    }
    self$paramIndex=rbindlist(pIndList,idcol='paramType')
    setkeyv(self$paramIndex,cols=c("paramType","chain","replicate"))
})

###
## Setter functions
###

parameterObject$set("public","setParamValue", function(paramType,value,chain=NULL,replicate=NULL){
    ind=self$getParamIndicies(paramType,chain,replicate)
    if(length(ind) %% length(value) != 0){
        stop('The number of parameter indicies must be a multiple of the length of the value vector when setting parameter values.')
    }
    self$params[ind]=rep(value,length.out=length(ind))
})

parameterObject$set("public","setFixedConstraint", function(paramType,value,chain=NULL,replicate=NULL){
    ind=self$getParamIndicies(paramType,chain,replicate)
    if(length(ind) %% length(value) != 0){
        stop('The number of parameter indicies must be a multiple of the length of the value vector when setting parameter values.')
    }
    self$fixed[ind]=rep(value,length.out=length(ind))
})

parameterObject$set("public","setUpperBoundConstraint", function(paramType,value,chain=NULL,replicate=NULL){
    ind=self$getParamIndicies(paramType,chain,replicate)
    if(length(ind) %% length(value) != 0){
        stop('The number of parameter indicies must be a multiple of the length of the value vector when setting parameter values.')
    }
    self$upperBound[ind]=rep(value,length.out=length(ind))
})

parameterObject$set("public","setLowerBoundConstraint", function(paramType,value,chain=NULL,replicate=NULL){
    ind=self$getParamIndicies(paramType,chain,replicate)
    if(length(ind) %% length(value) != 0){
        stop('The number of parameter indicies must be a multiple of the length of the value vector when setting parameter values.')
    }
    self$lowerBound[ind]=rep(value,length.out=length(ind))
})

## Sets entire parameter vector
parameterObject$set("public","setParamVector", function(params){
    self$params=numeric(max(self$paramIndex$end))
    ## Iterate over all parameters
    for(p in names(params)){
        self$setParamValue(p,params[[p]])
    }
})

## A function that takes the list of parameter constraints and uses it to set constraint variables in self
parameterObject$set("public","setParamConstraints", function(lowerBound=list(),upperBound=list(),fixed=list()){
    self$lowerBound=numeric(max(self$paramIndex$end))
    self$upperBound=numeric(max(self$paramIndex$end))
    self$fixed=numeric(max(self$paramIndex$end))
    ## Iterate over all parameters
    for(p in unique(self$paramIndex$paramType)){
        if(is.null(lowerBound[[p]])){
            lowerBound[[p]]=-Inf
        }
        if(is.null(upperBound[[p]])){
            upperBound[[p]]=Inf
        }
        if(is.null(fixed[[p]])){
            fixed[[p]]=FALSE
        }
        self$setLowerBoundConstraint(p,lowerBound[[p]])
        self$setUpperBoundConstraint(p,upperBound[[p]])
        self$setFixedConstraint(p,fixed[[p]])
    }
    self$fixed=as.logical(self$fixed)
})

####
## Parameter getter function
####
## Returns the indicies that a parameter exists in, in the parameter vector
## Note: Relies on the construction of the parameter vector where parameters of the same type and chain are contiguous
parameterObject$set("public","getParamIndicies", function(paramType,chain=NULL,replicate=NULL){
    ## Check that the parameter type is valid
    if(is.na(self$paramIndex[paramType]$group[1])){
        stop(paste("Invalid query for parameter of type",paramType))
    }
    ## If a chain was specified (but may be ignored)
    if(!is.null(chain)){
        ## If it is not a chain specific parameter
        if(self$paramIndex[eval(.(paramType))]$chain[1]==-1){
            chain=-1
        } else if(is.na(self$paramIndex[eval(.(paramType,chain))]$group[1])){
            stop("Out of range chain query")
        }
    } 
    ## If a replicate is specified (but may be ignored)
    if(!is.null(replicate)){
        ## If it is not a replicate specific parameter
        if(self$paramIndex[eval(.(paramType,chain))]$replicate[1]==-1){
            replicate=-1
        } else if(is.na(self$paramIndex[eval(.(paramType,chain,replicate))]$group[1])){
            stop("Out of range replicate query")
        }
    }
    qOut=self$paramIndex[eval(.(paramType,chain,replicate)),.(start,end)]
    
    out=do.call("c",foreach(z=1:nrow(qOut)) %do% {
        return(qOut$start[z]:qOut$end[z])
    })
    return(as.numeric(out))
})

parameterObject$set("public","getChainParameterMatrix", function(paramType,chain,replicate=NULL){
    ## Check that the parameter type is valid
    if(is.na(self$paramIndex[paramType]$group[1])){
        stop(paste("Invalid query for parameter of type",paramType))
    }
    ## If parameter is chain specific
    if(self$paramIndex[eval(.(paramType))]$chain[1]==-1){
        chain=-1
    } else if(is.na(self$paramIndex[eval(.(paramType,chain))]$group[1])){
        stop("Out of range chain query")
    }    
    ## If a replicate is specified (but may be ignored)
    if(!is.null(replicate)){
        ## If it is not a replicate specific parameter
        if(self$paramIndex[eval(.(paramType,chain))]$replicate[1]==-1){
            replicate=-1
        } else if(is.na(self$paramIndex[eval(.(paramType,chain,replicate))]$group[1])){
            stop("Out of range replicate query")
        }
    }
    qOut=self$paramIndex[eval(.(paramType,chain,replicate)),.(replicate,start,end)]
    pMat=matrix(0,ncol=max(abs(qOut$replicate)),nrow=qOut$end[1]-qOut$start[1]+1)
    for(z in 1:nrow(qOut)){
        pMat[,abs(qOut$replicate[z])]=self$params[qOut$start[z]:qOut$end[z]]
    }
    return(pMat)
})


###
## Convenience functions
###

## Returns a table of parameters for the HMM object
parameterObject$set("public","getParameterTable",function(){
    pTab=foreach(i=1:nrow(self$paramIndex)) %do% {
        with(self$paramIndex,data.table(paramName=rep(paramType[i],end[i]-start[i]+1),chain=rep(chain[i],end[i]-start[i]+1),
                                        replicate=rep(replicate[i],end[i]-start[i]+1),
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
