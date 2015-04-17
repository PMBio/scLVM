gpCLVM <- setClass(
  
  #  gpCLVM class

  
  "gpCLVM",
  
  slots = c(
    Y = "matrix",
    X= "matrix",
    K = "matrix",
    Ki = "matrix",
    var = "matrix",
    interaction = "logical"
  )
  
)

setMethod("initialize","gpCLVM",function(.Object,Y,X0=NULL,k=1,standardize=FALSE,interaction=TRUE){ 
  #set data
   gpCLVM_py(Y=Y,X0=X0,k=k,standardize=standardize,interaction=interaction) #push it to python and call python constructor
   if(standardize==TRUE){
   .Object@Y = as.matrix(scale(Y))}else{.Object@Y = as.matrix(Y)}
   .Object@interaction = as.logical(interaction)
  .Object
}
)


setGeneric(name = "optimize",
           def = function(obj,...){
             standardGeneric("optimize")
           })
            

setMethod(f = "optimize",
          signature = "gpCLVM",
          definition = function(obj){     
            res <- optimize_py(obj)
            obj@K <- as.matrix(res$K)
            obj@X <- as.matrix(res$X)
            #obj@var <- as.matrix(res$var)
            if(obj@interaction==TRUE){obj@Ki <- as.matrix(res$Ki)}
            
            return(obj)  
            
          }
)




setGeneric(name = "getK",
           def = function(obj,...){
             standardGeneric("getK")
           }
)

setMethod(f = "getK",
          signature = "gpCLVM",
          definition = function(obj){            
            K = obj@K
            return(K)              
          }
)

setGeneric(name = "getX",
           def = function(obj,...){
             standardGeneric("getX")
           }
)
setMethod(f = "getX",
          signature = "gpCLVM",
          definition = function(obj){            
            X = obj@X
            return(X)              
          }
)



setGeneric(name = "getKi",
           def = function(obj,...){
             standardGeneric("getKi")
           }
)
setMethod(f = "getKi",
          signature = "gpCLVM",
          definition = function(obj){
            if(obj@interaction==TRUE){
              Ki = obj@Ki
              return(Ki)}else{return(NULL)}              
          }
)


# setGeneric(name = "getVarianceComps",
#            def = function(obj,...){
#              standardGeneric("getVarianceComps")
#            }
# )
# setMethod(f = "getVarianceComps",
#           signature = "gpCLVM",
#           definition = function(obj){    
#             var = obj@var
#             return(var)              
#           }
# )
