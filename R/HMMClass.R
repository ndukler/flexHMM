HMM <- R6Class("HMM",public=list(emission="Emission",transition="Transition",logPrior="numeric",alphaTable="list",betaTable="list",logLiklihood="numeric"))

#' checkHMMValidity
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

## Check validity of HMM object
HMM$set("public","checkHMMValidity",function(){
    errors=character()
    if(!"Emission" %in% class(self$emission)){
        errors=c(errors,"The emission argument does not inherit from an object of type Emission")
    }
    if(!"Transition" %in% class(self$transition)){
        errors=c(errors,"The transition argument does not inherit from an object of type Transition")
    }
    ## If the object types are wrong, short circuit the rest of the checks
    if (length(errors) != 0) {return(errors)}
    if(self$emission$nstates!=self$transition$nstates){
        errors=c(errors,"The transition and emission models must have the same number of states.")
    }
    if (length(errors) == 0) TRUE else errors
})

## constructor function for HMM object
HMM$set("public","initialize", function(emission,transition){
    transition$updateTransitionProbabilities()
    ## Compute the prior distribution of states as the steady state distribution of the transition matrix
    leadingEigenV=eigen(exp(transition$transitionLogProb))$vectors[,1]
    steadyState=as.numeric(log(leadingEigenV/sum(leadingEigenV)))
    ## Set object values
    self$emission=emission
    self$transition=transition
    self$logPrior=steadyState
},overwrite=TRUE)

##
## Numeric tricks
##

## Log sum exp trick
logSumExp <- function(z){
    a=max(z)
    return(a+log(sum(exp(z-a))))
}

## Vectorized version  for matricies of log sum exp trick
logSumExpV <- function(x,byRow=FALSE){
    if(byRow){
        out=apply(x,1,logSumExp)
    } else {
        out=apply(x,2,logSumExp)
    }
    return(out)
}

###
## HMM algorithms
###

## Implement method for the forward algorithm
HMM$set("public","forwardAlgorithm", function(ncores=1){
    self$alphaTable=mclapply(self$emission$sumChainReplicates(),function(x){
        return(forwardAlgorithmCpp(x,self$transition$transitionLogProb,self$logPrior))
    },mc.cores=ncores)
},overwrite=TRUE)


## Implement method for the backward algorithm
HMM$set("public","backwardAlgorithm",function(ncores=1){
    self$betaTable=mclapply(self$emission$sumChainReplicates(),function(x){
        return(backwardAlgorithmCpp(x,self$transition$transitionLogProb))
    },mc.cores = ncores)
})

## Method to compute viterbi path
HMM$set("public","computeViterbiPath",function(){
    vPathList=lapply(self$alphaTable,function(x){
        vPath = numeric(nrow(x))
        vPath[length(vPath)] = which.max(x[length(vPath),])
        if(nrow(x)>1){
            for(j in (nrow(x)-1):1){
                vPath[j]=which.max(x[j,]+self$transition$transitionLogProb[,vPath[j+1]])
            }
        }
        return(vPath)
    })
    return(vPathList)
})

## Method to compute liklihood from alpha table
HMM$set("public","computeLogLiklihood",function(){
    self$logLiklihood=sum(unlist(lapply(self$alphaTable,function(x) logSumExp(x[nrow(x),])),use.names=FALSE))
    ## using beta table
    ## logSumExp(hmm@betaTable[[1]][1,]+hmm@emission@emissionLogProb[[1]][1,]+hmm@logPrior)
})

## Method to compute marginal probability table
HMM$set("public","computeMarginalProb",function(){
    mProb=list()
    for(i in 1:length(self$alphaTable)){
        mProb[[paste0(as.character(names(self$alphaTable)[i]),i)]]=exp((self$alphaTable[[i]]+self$betaTable[[i]])-logSumExpV(self$alphaTable[[i]]+self$betaTable[[i]],TRUE))
    }
    return(mProb)
})

## Method to compute marginal probability table
HMM$set("public","getParameterTable",function(){
    return(rbind(self$emission$getParameterTable()[,pType:="emission"],self$transition$getParameterTable()[,pType:="transition"]))
})

## Function to pass to optim, that takes a set of values, and an HMM, and returns a negative log liklihood
updateAllParams <- function(x,hmmObj,nthreads){
    if(length(x) != sum(!hmmObj$transition$fixed)+sum(!hmmObj$emission$fixed) ){
        stop("Length of x does not match the number of expected parameters")
    }
    ## Take in parameters
    hmmObj$emission$params[!hmmObj$emission$fixed]=x[1:sum(!hmmObj$emission$fixed)]
    hmmObj$transition$params[!hmmObj$transition$fixed]=x[1:sum(!hmmObj$transition$fixed)+sum(!hmmObj$emission$fixed)]
    ## Update emission and transition probabilities
    hmmObj$emission$updateEmissionProbabilities()
    hmmObj$transition$updateTransitionProbabilities()
    ## Run forward algorithm
    hmmObj$forwardAlgorithm(nthreads)
    hmmObj$computeLogLiklihood()
    ## return(-hmmObj$logLiklihood)
    return(hmmObj$logLiklihood)
}

## Method to fit HMM
fitHMM <- function(hmm,nthreads=1){
    ## Pass non-fixed parameters for optimization
    ## optimx(par=c(hmm$emission$params[!hmm$emission$fixed],hmm$transition$params[!hmm$transition$fixed]),fn=updateAllParams,
    ##       lower=c(hmm$emission$lowerBound[!hmm$emission$fixed],hmm$transition$lowerBound[!hmm$transition$fixed]),
    ##       upper=c(hmm$emission$upperBound[!hmm$emission$fixed],hmm$transition$upperBound[!hmm$transition$fixed]),
    ##       method=c("L-BFGS-B"),itnmax=100,control=list(trace=1),hmmObj=hmm,nthreads=nthreads)
    rphast::optim.rphast(func=updateAllParams, params=c(hmm$emission$params[!hmm$emission$fixed],hmm$transition$params[!hmm$transition$fixed]),
                 lower = c(hmm$emission$lowerBound[!hmm$emission$fixed],hmm$transition$lowerBound[!hmm$transition$fixed]),
                 upper = c(hmm$emission$upperBound[!hmm$emission$fixed],hmm$transition$upperBound[!hmm$transition$fixed]),
                 logfile="foo.log",hmmObj=hmm,nthreads=nthreads)
}

setGeneric("plot.hmm",function(hmm=NULL,viterbi=NULL,marginal=NULL,truePath=NULL,start=NA_real_,end=NA_real_,chain=1){ standardGeneric("plot.hmm") })
## Now some methods to plot HMM
setMethod("plot.hmm",signature=c(hmm="ANY",viterbi="ANY",marginal="ANY",truePath="ANY",start="numeric",end="numeric",chain="ANY"),
          definition=function(hmm,viterbi=NULL,marginal=NULL,truePath=NULL,start=NA_real_,end=NA_real_,chain=1){
              ## If start and end not set, use full data length as default
              if(is.na(start))
                  start=1
              if(is.na(end))
                  end=nrow(hmm$emission$emissionLogProb[[chain]])
              if(start>end){
                  stop("Start must be less than end")
              }
              ## Add elements one layer at a time
              g=ggplot()
              if(!is.null(hmm$emission$emissionLogProb)){
                  dat=rbindlist(lapply(hmm$emission$data[[chain]],function(x){
                      x=as.data.table(x)
                      x[,index:=1:nrow(x)]
                      melt(x[index>=start & index <=end],id.vars="index")
                  }),idcol=TRUE)

                  dat[,variable:=NULL]
                  dat[,type:="Raw data"]
                  g=g+geom_point(data=dat,aes(x=index,y=value,color=.id),alpha=1/3,inherit.aes=FALSE)
              }
              if(!is.null(viterbi)){
                  vit=data.table(index=start:end,value=viterbi[[chain]][start:end],type="Viterbi")
                  tic=data.table(x=start,x1=end,y=1:max(hmm$transition$nstates),type="Viterbi")
                  g=g+geom_segment(data=tic,aes(x=x,y=y,xend=end,yend=y),linetype=2,color="gray",inherit.aes=FALSE)+
                      geom_line(data=vit,aes(x=index,y=value,linetype="Viterbi"),color="black",inherit.aes=FALSE)
              }
              if(!is.null(truePath)){
                  tru=data.table(index=start:end,value=truePath[[chain]][start:end],type="Viterbi")
                  g=g+geom_line(data=tru,aes(x=index,y=value,linetype="True Path"),color="red",inherit.aes=FALSE)+
                      guides(linetype=guide_legend(title="Path Type"))+
                      scale_linetype_manual(values=c(3,1))

              }
              if(!is.null(marginal)){
                  mar=data.table(index=start:end,value=marginal[[chain]][start:end,])
                  mar=melt(mar,id.vars=c("index"))
                  mar[,variable:=gsub("value.V","class.",variable)]
                  mar[,type:="Marginal"]
                  g=g+geom_bar(data=mar,aes(x=index,y=value,fill=variable),stat="identity",inherit.aes=FALSE)
              }
              g=g+facet_wrap(~type,ncol=1,scales="free_y")+
                  theme_bw()
              return(g)
})

