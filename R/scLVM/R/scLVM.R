scLVM <- setClass(

#  scLVM class
#  slots:
#       geneID:                 G vector of geneIDs
#       Y:                      gene expression matrix [N, G]
#       tech_noise:             G vector of tech_noise
 
  "scLVM",
  
  slots = c(
    Y = "matrix",
    geneID = "character",
    tech_noise= "numeric",
    Ycorr = "matrix",
    var = "matrix",
    conv = "logical"
    )
  
)

setGeneric(name = "init",
           def = function(Object,Y,tech_noise=NULL){
             standardGeneric("init")
           })

setMethod("init","scLVM",function(Object,Y,tech_noise=NULL){ 
  require(rPython)
  #set data
  objName = deparse(substitute(Object)) #to make sure scLVM objects have the same name in R and python
  geneID = colnames(Y)
  if(is.null(geneID)){stop('Provide gene IDs as colnames of Y')}

  scLVM_py(objName,Y, geneID, tech_noise) #push it to python and call python constructor
  Object@Y = as.matrix(Y)
  Object@geneID = as.character(geneID)
  Object@tech_noise = as.numeric(tech_noise)
  Object
}
)


setGeneric(name = "fitFactor",
           def = function(Object,idx=NULL,geneSet = NULL, XKnown = NULL, k=1,standardize=FALSE, use_ard=FALSE,interaction=TRUE, initMethod='fast'){
             standardGeneric("fitFactor")
           })



setMethod(f = "fitFactor",
          signature = "scLVM",
          definition = function(Object,idx=NULL,geneSet = NULL, XKnown = NULL, k=1,standardize=FALSE, use_ard=FALSE,interaction=TRUE, initMethod='fast'){
              objName = deparse(substitute(Object))
            if(is.null(idx) + is.null(geneSet) != 1)stop("Provide either gene identifiers (geneSet) OR indices (idx)")
            
            if(!is.null(geneSet)){
              idx_geneSet <- na.omit(match(ens_ids_cc, geneID))
              idx = intersect(which(colMeans(Object@Y)>0),idx_geneSet)
              if(length(idx)==0){stop("Couldn't find any matches between geneSet and geneIDs. Make sure they are of the same type (ENSEMBL if possible)")}
            }
            
              res = fitLatent_py(objName,idx,XKnown,k,standardize, use_ard,interaction=interaction, initMethod=initMethod)
              return(res)  

          }
)

setGeneric(name = "setTechnicalNoise",
           def = function(obj,tech_noise=NULL){
             standardGeneric("setTechnicalNoise")
           }
)
setMethod(f = "setTechnicalNoise",
          signature = "scLVM",
          definition = function(obj,tech_noise=NULL){
            objName = deparse(substitute(obj))
            python.assign("tech_noise",tech_noise)
            python.exec(paste(objName,".set_tech_noise(SP.array(tech_noise))",sep=""))
            obj@tech_noise = tech_noise
            return(obj)              
          }
)


setGeneric(name = "varianceDecomposition",
           def = function(Object,K=NULL,idx = NULL){
             standardGeneric("varianceDecomposition")
           }
           )
setMethod(f = "varianceDecomposition",
          signature = "scLVM",
          definition = function(Object,K=NULL,idx = NULL){
            objName = deparse(substitute(Object))
            varianceDecomposition_py(objName,K,idx)
            res <- getVarianceComponents_py(objName)
            Object@var = res$var
            Object@conv = res$conv
            return(Object)  
            
          }
)

setGeneric(name = "setVarianceComponents",
           def = function(obj,normalize=TRUE){
             standardGeneric("setVarianceComponents")
           }
)
setMethod(f = "setVarianceComponents",
          signature = "scLVM",
          definition = function(obj,normalize=TRUE){  
          objName = deparse(substitute(obj))
          res <- getVarianceComponents_py(objName,normalize)
          obj@var = res$var
          obj@conv = res$conv
          return(obj)              
          }
)

setGeneric(name = "getVarianceComponents",
           def = function(obj){
             standardGeneric("getVarianceComponents")
           }
)
setMethod(f = "getVarianceComponents",
          signature = "scLVM",
          definition = function(obj){            
          res = list()
          res$var = obj@var
          res$conv = obj@conv
          return(res)              
          }
)

setGeneric(name = "setCorrectedExpression",
           def = function(obj,rand_eff_ids=NULL){
             standardGeneric("setCorrectedExpression")
           }
)
setMethod(f = "setCorrectedExpression",
          signature = "scLVM",
          definition = function(obj,rand_eff_ids=NULL){
            objName = deparse(substitute(obj))
            Ycorr = getCorrectedExpression_py(objName,rand_eff_ids)
            obj@Ycorr = Ycorr
            return(obj)              
          }
)

setGeneric(name = "getCorrectedExpression",
           def = function(obj,rand_eff_ids=NULL){
             standardGeneric("getCorrectedExpression")
           }
)
setMethod(f = "getCorrectedExpression",
          signature = "scLVM",
          definition = function(obj,rand_eff_ids=NULL){
            objName = deparse(substitute(obj))
            Ycorr = getCorrectedExpression_py(objName,rand_eff_ids)
            obj@Ycorr = Ycorr
            return(Ycorr)              
          }
)

setGeneric(name = "LMM",
           def = function(obj, expr = NULL, K=NULL, idx=NULL,geneSet=NULL, verbose=TRUE, recalc=TRUE){
             standardGeneric("LMM")
           }
)
setMethod(f = "LMM",
          signature = "scLVM",
          definition = function(obj,expr = NULL, K=NULL,idx=NULL,geneSet=NULL,verbose=TRUE,recalc=TRUE){   
            if(is.null(idx) + is.null(geneSet) != 1)stop("Provide either gene identifiers (geneSet) OR indices (idx)")
            
            if(!is.null(geneSet)){
              idx_geneSet <- na.omit(match(ens_ids_cc, geneID))
              idx = intersect(which(colMeans(Object@Y)>0),idx_geneSet)
              if(length(idx)==0){stop("Couldn't find any matches between geneSet and geneIDs. Make sure they are of the same type (ENSEMBL if possible)")}
            }
            
          objName = deparse(substitute(obj))
          res= fitLMM_py(objName,expr,K,idx,verbose,recalc)
            return(res)              
          }
)