scLVM_py = function(objName,Y=NULL,geneID=NULL,tech_noise=NULL){
  
  row.names(Y)=c() #strip Y of column and row names
  colnames(Y)=c()
  python.assign("Y",Y) #outfile is the hdf file generated in the previous setting
  python.exec("Y = SP.array(Y)") #matrices are passed as lists and need to be converted back to mtrices
  python.assign("tech_noise",tech_noise)
  if(!is.null(tech_noise)){
    python.exec("tech_noise = SP.array(tech_noise)")}
  python.assign("geneID",geneID)  
  
  python.exec(paste(objName," = scLVM(Y,geneID=geneID,tech_noise=tech_noise)", sep=""))
  
}


fitLatent_py = function(objName, idx=NULL, XKnown = NULL, k=1,standardize=FALSE, use_ard=FALSE,interaction=TRUE, initMethod='fast'){
  python.assign("idx",idx-1)#R indexing to python indeces
  python.assign("X0",XKnown)
  
  if(!is.null(XKnown)){
    python.exec("X0 = SP.array(X0)")
    if(is.null(dim(XKnown))){
      python.exec("X0 = SP.reshape(X0,(len(X0),1))")}
  }
  
  python.assign("k",k)
  python.assign("standardize",standardize)
  python.assign("use_ard",use_ard)
  python.assign("interaction",interaction)
  python.assign("initMethod",initMethod)
  python.exec(paste("X,K,Kint,varGPLVM_ARD = ",objName,".fitGPLVM(idx=idx, X0=X0, k=k,standardize=standardize, use_ard=use_ard, interaction=interaction, initMethod=initMethod)",sep=""))
  
  res=list()
  if(use_ard==TRUE){
    python.exec("Xard = list(varGPLVM_ARD['X_ARD'])")
    X_ard =  python.get("Xard")    
    res$X_ard = X_ard
  }
  
  if(k>1){
    python.exec("X = X.tolist()")
    X = do.call(rbind,python.get("X"))
    res$X = X
    
  }else{
    python.exec("X = X.tolist()")
    X =  python.get("X")    
    res$X = unlist(X)
  }
  
  python.exec("K = K.tolist()")
  K = do.call(rbind,python.get("K"))
  res$K = K
  
  if(interaction==TRUE & !is.null(XKnown)){
    python.exec("Kint = Kint.tolist()")
    Kint = do.call(rbind,python.get("Kint"))
    res$Kint = Kint
  }
  
  return(res)
}

varianceDecomposition_py = function(objName,K=NULL,idx=NULL){
  if(is.matrix(K)){
    K = as.matrix(K)
    row.names(K)=c() #strip Y of column and row names
    colnames(K)=c()
    python.assign("K",K)
    python.exec("K = SP.array(K)")   
  }else{
    python.exec("K = []")
    for(i in 1:length(K)){
      K_ = as.matrix(K[[i]])
      row.names(K_)=c() #strip Y of column and row names
      colnames(K_)=c()
      python.assign("K_",K_) 
      python.exec("K.append(SP.array(K_))")      
    }        
  }
  
  if(!is.null(idx)){idx = as.integer(idx-1)}
  
  python.assign("idx", idx)
  
  python.exec(paste(objName,".varianceDecomposition(K=K,idx=idx)", sep=""))
}


getVarianceComponents_py = function(objName,normalize=TRUE){
  python.assign("normalize", normalize)
  python.exec(paste('var, var_info = ',objName,'.getVarianceComponents(normalize=normalize)', sep=""))
  python.exec("conv = var_info['conv'].tolist()")
  gene_idx = python.get("var_info['gene_idx'].tolist()")
  geneID = python.get("geneID")
  res = list()
  res$conv = python.get("conv")
  res$var = do.call(rbind,python.get("var.tolist()"))
  colnames(res$var) = python.get("var_info['col_header'].tolist()")
  row.names(res$var) = geneID[gene_idx+1]
  python.assign(paste("geneID_vd",objName,sep="_"),geneID[gene_idx+1])
  
  return(res)
}


getCorrectedExpression_py = function(objName,rand_eff_ids=NULL){
  if(!is.null(rand_eff_ids))rand_eff_ids = rand_eff_ids-1
  python.assign("rand_eff_ids", rand_eff_ids)
  python.exec(paste("Ycorr = ",objName,".getCorrectedExpression(rand_eff_ids)", sep=""))
  Ycorr = do.call(rbind, python.get("Ycorr.tolist()"))
  colnames(Ycorr) <- python.get(paste("geneID_vd",objName,sep="_"))
  return(Ycorr)
}

fitLMM_py <- function(objName,expr = NULL, K=NULL,idx=NULL,verbose=TRUE, recalc=TRUE){
  
  if(!is.null(idx)){idx = as.integer(idx-1)}
  python.assign("idx", idx)
  
  row.names(expr)=c() 
  colnames(expr)=c()
  names(expr) = c()
  python.assign("expr",expr) 
  if(!is.null(expr)){
    python.exec("expr = SP.array(expr)")
    if(is.null(dim(expr))){
      python.exec("expr = SP.reshape(expr,(len(expr),1))")}
  }
  
  python.assign("verbose", verbose)
  python.assign("recalc", recalc)
  if(!is.null(K)){
    if(is.matrix(K)){
      K = as.matrix(K)
      row.names(K)=c() #strip Y of column and row names
      colnames(K)=c()
      python.assign("K",K)
      python.exec("K = SP.array(K)")   
    }else{
      python.exec("K = []")
      for(i in 1:length(K)){
        K_ = as.matrix(K[[i]])
        row.names(K_)=c() #strip Y of column and row names
        colnames(K_)=c()
        python.assign("K_",K_) 
        python.exec("K.append(SP.array(K_))")      
      }        
    }
  }else{python.assign("K",K)}
  
  python.exec(paste("pv,beta,info = ",objName,".fitLMM(expr = expr, K=K, idx=idx, verbose=verbose)", sep=""))
  python.exec("beta[SP.where(SP.isnan(beta))]=0")#nans can happen for K=None and they are trouble in JSON
  python.exec("pv[SP.where(SP.isnan(pv))]=1")
  pv = do.call(rbind,python.get("pv.tolist()"))
  beta = do.call(rbind,python.get("beta.tolist()"))
  
  gene_idx_row = python.get("info['gene_idx_row'].tolist()")+1
  
  row.names(beta) = geneID[gene_idx_row]
  row.names(pv) = geneID[gene_idx_row]
  
  if(is.null(expr)){
    colnames(beta) = geneID[gene_idx_row]
    colnames(pv) = geneID[gene_idx_row]
  }else{
    colnames(beta) = colnames(expr)
    colnames(pv) = colnames(expr)    
  }
  
  res = list()
  res$beta = beta
  res$pv = pv
  res$gene_idx = gene_idx_row
  return(res)
}



# #####gpCLVM#####
# 
# gpCLVM_py = function(Y=NULL,X0=NULL,k=1,standardize=FALSE,interaction=TRUE){
#   
#   row.names(Y)=c() #strip Y of column and row names
#   colnames(Y)=c()
#   python.assign("Y",Y) #outfile is the hdf file generated in the previous setting
#   python.exec("Y = SP.array(Y)") #matrices are passed as lists and need to be converted back to mtrices
#   python.assign("X0",X0)
#   if(!is.null(X0)){
#     python.exec("X0 = SP.array(X0)")}
#   if(is.null(dim(X0))){
#     python.exec("X0 = SP.reshape(X0,(len(X0),1))")}
#   python.assign("k",as.integer(k))
#   python.assign("standardize",standardize)
#   python.assign("interaction",interaction)
#   
#   python.exec("gp = gpCLVM(Y=Y,X0=X0,k=k,standardize=standardize,interaction=interaction)")
#   
# }
# 
# 
# 
# optimize_py = function(obj){
#   python.exec("params0 = gp.initParams()")
#   python.exec("conv = gp.optimize(params0)")
#   python.exec("X1 = gp.getX()")
#   python.exec("K1 = gp.getK()")
#   
#   
#   res=list()
#   python.exec("X1 = X1.tolist()")
#   X =  do.call(rbind,python.get("X1"))
#   res$X = X
#   
#   python.exec("K1 = K1.tolist()")
#   K = do.call(rbind,python.get("K1"))
#   res$K = K
#   
#   if(obj@interaction==TRUE){
#     python.exec("Ki = gp.getKi()")
#     python.exec("Ki = Ki.tolist()")
#     Ki = do.call(rbind,python.get("Ki"))
#     res$Ki = Ki
#   }
#   
#   #python.exec("var = gpCLVM.getVarianceComps()")
#   #res$var = var
#   #var = do.call(rbind,python.get("var.tolist()"))
#   
#   res
# }
