###
## Create generic emission class
###

## @data is in the form of a list of a list of matricies. The first layer of the list represents the chain index, the second layer
##       represents the feature index. The rows of the matricies are positions along the chain, and columns are replicates.
## @emissionLogProb is a list, where each  matrix in the list is a chain,each column is a state, and each row is an index
## @nstates defines the number of states in the HMM
## @chainReplicates is a 

Emission <- R6::R6Class("Emission",public=list(data=list(),emissionLogProb=list(),nstates=numeric(),chainReplicates=list()),inherit=parameterObject)

## Add generic constructor for all emission type objects
Emission$set("public","initialize",function(data,nstates,params=list(),lowerBound=list(), upperBound=list(),fixed=list(),
                                              invariants=list(),chainSpecificParameters=NULL,replicatesPerParameter=numeric()){
    ## Set number of states
    self$nstates=nstates

    ## First check validity for items with generic validity checkers
    self$checkDataValidity(data)
    self$data=data
    self$checkParamConstraints(params,lowerBound,upperBound,fixed)

    ## Invariant check with no dependencies
    if(is.list(invariants)){
        self$invariants=invariants
    } else stop("Invariants must be in a list.")
   
    ## Perform parameter validity check and set parameters
    self$checkParamValidity(params)
    
    ## Set parameter vector
    self$buildParamIndex(params,chainSpecificParameters,replicatesPerParameter)
    self$setParamVector(params)

    ## Reshape parameter constraints to vectors and set self$lowerBound, self$upperBound, self$fixed, must be after params are set
    self$setParamConstraints(lowerBound,upperBound,fixed)
  
    ## Now perform additional checks for invariant validity
    self$checkInvariantValidity(invariants)
    self$checkEmissionValidity()
})

###
## Functions for checking validity of emission inputs
###

## A function that checks for errors during emission construction
Emission$set("public","checkEmissionValidity", function(){
    errors=character()
    ## Check for general construction errors
    if(sum(!is.logical(self$fixed))>0){
        errors=c(errors,"Non-logical values in fixed.")
    }
    if(!is.numeric(self$nstates)){
        errors=c(errors,"Non-numeric number of states")
    }
    if(sum(!is.numeric(self$params) | is.na(self$params))>0){
        errors=c(errors,"Non-numeric value in parameter list.")
    }
    if (length(errors) == 0) TRUE else errors
})

## Function that checks for chain replicate validity
Emission$set("public","checkChainReplicateValidity", function(data,chainReplicates){
    if(sum(duplicated(unlist(chainReplicates)))>0){
        stop(paste("Chain(s)",unique(unlist(chainReplicates)[duplicated(unlist(chainReplicates))]),"are duplicated"))
    }
    if(sum(!1:length(data) %in% unlist(chainReplicates))>0){
        stop("Not all chains are present in chainReplicates")
    }
    if(sum(unlist(chainReplicates) %in% !1:length(data))>0){
        stop("Not out of range chain specified in chainReplicates")
    }
    ## If some chains are replicates of each other, check that they have the same length
    if(!is.null(chainReplicates)){
        ## Get the lengths of chains that are supposed to be replicates and place each group of replicates into a matrix where the column
        ## is a feature and the row is a replicate. Each list element is one group of replicates. 
        replicate.lengths=lapply(chainReplicates, function(x){do.call("rbind",lapply(data[as.numeric(x)],function(y) matrix(unlist(lapply(y,nrow)),byrow=TRUE,nrow=1)))})
        for(r in 1:length(replicate.lengths)){
            non.ident=apply(replicate.lengths[[i]],2, function(y) sum(length(unique(y))!=1))
            if(sum(non.ident)>0){
                stop(paste("For replicate group",i,"feature(s)",which(non.ident!=0),"do not have the same number of rows."))
            }
        }
    }
})

## A function that is used to check data structure and content validity
Emission$set("public","checkDataValidity", function(data,chainReplicates){
    errors=character()
    ## Check that the data is formated correctly, critical errors that will stop execution immediately
    if(length(data)==0){
        stop("Data must be set to a list of non-zero length")
    }
    if(class(data)!="list"){
        stop("The  first layer of data must be a list, each element representing a seperate chain.")
    }
    if(sum(unlist(lapply(data,function(x) class(x)[1]!="list")))>0){
        stop("The second layer of data must be lists, each element representing a seperate feature.")
    }
    ## Data content error
    if(sum(unlist(lapply(data,function(x) lapply(x, function(y) !is.numeric(y)))))>0){
        stop("Non-numeric value in data matrix.")
    }
 })

## Update log emission probability
Emission$set("public","updateEmissionProbabilities", function(){
        stop("Must implement a updateEmissionProbabilities function.")
})
