HMM <- R6::R6Class("HMM",public=list(emission="Emission",transition="Transition",logPrior="numeric",alphaTable="list",betaTable="list",logLiklihood="numeric"))

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
    ## Check which transitions are permissible that end in state X
    permTrans=lapply(split(self$transition$transitionLogProb,seq(ncol(self$transition$transitionLogProb))),function(x) which(!is.infinite(x))-1)
    permTransV=unlist(permTrans,use.names=FALSE)
    ## If matrix is more than 50% sparse use sparse 
    if(length(permTransV) <= 0.5 * self$emission$nstates^2){
        tLen=unlist(lapply(permTrans,length),use.names=FALSE)
        ## Compute forward table in parallel for each chain
        self$alphaTable=mclapply(self$emission$emissionLogProb,function(x){
            return(forwardAlgorithmSparseCpp(x,self$transition$transitionLogProb,self$logPrior,permTransV,tLen))
        },mc.cores=ncores)
    } else {
        self$alphaTable=mclapply(self$emission$emissionLogProb,function(x){
            return(forwardAlgorithmCpp(x,self$transition$transitionLogProb,self$logPrior))
        },mc.cores=ncores)
    }
},overwrite=TRUE)

## Implement method for the backward algorithm
HMM$set("public","backwardAlgorithm",function(ncores=1){
    ## Check which transitions are permissible that end in state X (reverse)
    permTrans=lapply(split(self$transition$transitionLogProb,seq(nrow(self$transition$transitionLogProb))),function(x) which(!is.infinite(x))-1)
    permTransV=unlist(permTrans,use.names=FALSE)
    ## If matrix is more than 25% sparse use sparse 
    if(length(permTransV) <= 0.75 * self$emission$nstates^2){
        tLen=unlist(lapply(permTrans,length),use.names=FALSE)
        ## Compute forward table in parallel for each chain
        self$betaTable=mclapply(self$emission$emissionLogProb,function(x){
            return(backwardAlgorithmSparseCpp(x,self$transition$transitionLogProb,permTransV,tLen))
        },mc.cores=ncores)
    } else {
        self$betaTable=mclapply(self$emission$emissionLogProb,function(x){
            return(backwardAlgorithmCpp(x,self$transition$transitionLogProb))
        },mc.cores=ncores)
    }
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
    if(!sum(hmmObj$emission$fixed)>0 || length(hmmObj$emission$emissionLogProb)==0){
        hmmObj$emission$updateEmissionProbabilities()
    }
    if(!sum(hmmObj$transition$fixed)>0 || length(hmmObj$transition$transitionLogProb)==0 ){
        hmmObj$transition$updateTransitionProbabilities()
    }
    ## Run forward algorithm
    hmmObj$forwardAlgorithm(nthreads)
    hmmObj$computeLogLiklihood()
    ## return(-hmmObj$logLiklihood)
    return(hmmObj$logLiklihood)
}

## Method to fit HMM
fitHMM <- function(hmm,nthreads=1){
    ## Pass non-fixed parameters for optimization
    rphast::optim.rphast(func=updateAllParams, params=c(hmm$emission$params[!hmm$emission$fixed],hmm$transition$params[!hmm$transition$fixed]),
                 lower = c(hmm$emission$lowerBound[!hmm$emission$fixed],hmm$transition$lowerBound[!hmm$transition$fixed]),
                 upper = c(hmm$emission$upperBound[!hmm$emission$fixed],hmm$transition$upperBound[!hmm$transition$fixed]),
                 logfile="foo.log",hmmObj=hmm,nthreads=nthreads)
}

setGeneric("plot.hmm",function(hmm=NULL,viterbi=NULL,marginal=NULL,truePath=NULL,misc=NULL,start=NA_real_,end=NA_real_,chain=1){ standardGeneric("plot.hmm") })
## Now some methods to plot HMM
setMethod("plot.hmm",signature=c(hmm="ANY",viterbi="ANY",marginal="ANY",truePath="ANY",start="numeric",end="numeric",chain="ANY"),
          definition=function(hmm,viterbi=NULL,marginal=NULL,truePath=NULL,start=NA_real_,end=NA_real_,chain=1){
              ## If start and end not set, use full data length as default
              if(chain < 1 | !is.integer(chain)){
                  stop("Chain value must be an integer greater than 0.")
              }
              if(is.na(start))
                  start=1
              if(is.na(end))
                  end=nrow(hmm$emission$emissionLogProb[[chain]])
              if(start>end){
                  stop("Start must be less than end")
              }
              ## Initialize the plots
              g.dat=ggplot2::ggplot()
              g.path=ggplot2::ggplot()
              g.mar=ggplot2::ggplot()
              ## Add elements one layer at a time
              if(!is.null(hmm$emission$emissionLogProb)){
                  dat=data.table::rbindlist(lapply(hmm$emission$data[[chain]],function(x){
                      x=as.data.table(x)
                      x[,index:=1:nrow(x)]
                      melt(x[index>=start & index <=end],id.vars="index")
                  }),idcol=TRUE)
                  dat[,.id:=factor(.id,levels=names(hmm$emission$data[[chain]]))]
                  dat[,variable:=NULL]
                  g.dat=g.dat+
                      ggplot2::geom_raster(data=dat,aes(x=index,y=.id,fill=value),alpha=1/3,inherit.aes=FALSE)                
              }
              if(!is.null(viterbi)){
                  tic=data.table::data.table(x=start,x1=end,y=1:max(hmm$transition$nstates))
                  vit=data.table::data.table(index=start:end,value=viterbi[[chain]][start:end])
                  g.path=g.path+
                      ggplot2::geom_segment(data=tic,aes(x=x,y=y,xend=end,yend=y),linetype=2,color="gray",inherit.aes=FALSE)+
                      ggplot2::geom_line(data=vit,aes(x=index,y=value,linetype="Viterbi"),color="black",inherit.aes=FALSE)
              }
              if(!is.null(truePath)){
                  tru=data.table::data.table(index=start:end,value=truePath[[chain]][start:end],type="Viterbi")
                  g.path=g.path+
                      ggplot2::geom_line(data=tru,aes(x=index,y=value,linetype="True Path"),color="red",inherit.aes=FALSE)+
                      ggplot2::guides(linetype=guide_legend(title="Path Type"))+
                      ggplot2::scale_linetype_manual(values=c(3,1))

              }
              if(!is.null(marginal)){
                  mar=data.table(index=start:end,value=marginal[[chain]][start:end,])
                  mar=melt(mar,id.vars=c("index"))
                  mar[,variable:=gsub("value.V","class.",variable)]
                  g.mar=g.mar+
                      ggplot2::geom_bar(data=mar,aes(x=index,y=value,fill=variable),stat="identity",inherit.aes=FALSE)
              }
              ## Allow for the addition of misc. information if it is a matrix or vector of the same length as the
              ## data. If it is a matrix it will be melted
              if(!is.null(misc)){
                  if(length(misc)==nrow(hmm$emission$data[[chain]][[1]])){
                      if(is.numeric(misc) | ncol(misc)==1){
                          misc.sub=data.table::data.table(index=start:end,variable=as.factor(1),value=misc[start:end])
                      }else if(is.matrix(misc)){
                          misc.sub=data.table::data.table(index=start:end,reshape2::melt(misc[start:end,])[,2:3])
                          data.table::setnames(misc.sub,colnames(misc.sub)[2],"variable")
                          misc.sub[,variable:=as.factor(variable)]
                      }
                  }
                  g.misc=ggplot()+
                      ggplot2::geom_raster(data=misc.sub,aes(x=index,y=variable,fill=value),inherit.aes=FALSE)

              }
              plots=paste(c("g.dat","g.mar","g.path","g.misc")[!c(is.null(TRUE),is.null(marginal),is.null(viterbi),is.null(misc))],collapse=",")
              eval(parse(text=paste0("cowplot::plot_grid(",plots,",ncol=1, align = 'v')")))
              return(g)
})

