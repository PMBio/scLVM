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

setMethod("initialize","scLVM",function(.Object,Y,geneID=NULL,tech_noise=NULL){  
  #set data
  scLVM_py(Y, geneID, tech_noise) #push it to python and call python constructor
  .Object@Y = as.matrix(Y)
  .Object@geneID = as.character(geneID)
  .Object@tech_noise = as.numeric(tech_noise)
  .Object
}
)


setGeneric(name = "fitGPLVM",
           def = function(obj,...){
             standardGeneric("fitGPLVM")
           })



setMethod(f = "fitGPLVM",
          signature = "scLVM",
          definition = function(obj,idx=NULL,k=1,standardize=FALSE,out_dir='./cache',
                          file_name=NULL,recalc=FALSE, use_ard=FALSE){
            
  res = fitGPLVM_py(idx,k,standardize,out_dir,file_name,recalc, use_ard)
  return(res)  

          }
)

setGeneric(name = "varianceDecomposition",
           def = function(obj,...){
             standardGeneric("varianceDecomposition")
           }
           )
setMethod(f = "varianceDecomposition",
          signature = "scLVM",
          definition = function(obj,K=NULL,i0=1,i1=1){            
            varianceDecomposition_py(K,i0,i1)
            obj <- setVarianceComponents(obj, normalize=TRUE)
            obj <- setCorrectedExpression(obj)
            return(obj)              
          }
)

setGeneric(name = "setVarianceComponents",
           def = function(obj,...){
             standardGeneric("setVarianceComponents")
           }
)
setMethod(f = "setVarianceComponents",
          signature = "scLVM",
          definition = function(obj,normalize=TRUE){            
          res <- getVarianceComponents_py(normalize)
          obj@var = res$var
          obj@conv = res$conv
          return(obj)              
          }
)

setGeneric(name = "getVarianceComponents",
           def = function(obj,...){
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
           def = function(obj,...){
             standardGeneric("setCorrectedExpression")
           }
)
setMethod(f = "setCorrectedExpression",
          signature = "scLVM",
          definition = function(obj){            
            Ycorr = getCorrectedExpression_py()
            obj@Ycorr = Ycorr
            return(obj)              
          }
)

setGeneric(name = "getCorrectedExpression",
           def = function(obj,...){
             standardGeneric("getCorrectedExpression")
           }
)
setMethod(f = "getCorrectedExpression",
          signature = "scLVM",
          definition = function(obj){            
            Ycorr = obj@Ycorr
            return(Ycorr)              
          }
)

setGeneric(name = "fitLMM",
           def = function(obj,...){
             standardGeneric("fitLMM")
           }
)
setMethod(f = "fitLMM",
          signature = "scLVM",
          definition = function(obj,K=NULL,i0=i0,i1=i1,verbose=TRUE, geneID=NULL){            
          res= fitLMM_py(K,i0,i1,verbose, geneID)
            return(res)              
          }
)