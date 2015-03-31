scLVM_py = function(Y=NULL,geneID=NULL,tech_noise=NULL){

  row.names(Y)=c() #strip Y of column and row names
  colnames(Y)=c()
  python.assign("Y",Y) #outfile is the hdf file generated in the previous setting
  python.exec("Y = SP.array(Y)") #matrices are passed as lists and need to be converted back to mtrices
  python.assign("tech_noise",tech_noise)
  if(!is.null(tech_noise)){
    python.exec("tech_noise = SP.array(tech_noise)")}
  python.assign("geneID",geneID)  
  
  python.exec("sclvm = scLVM(Y,geneID=geneID,tech_noise=tech_noise)")
  
}


fitGPLVM_py = function(idx=NULL,k=1,standardize=FALSE,out_dir='./cache',file_name=NULL,recalc=FALSE, use_ard=FALSE){
  python.assign("idx",idx-1)#R indexing to python indeces
  python.assign("k",k)
  python.assign("standardize",standardize)
  python.assign("out_dir",out_dir)
  python.assign("file_name",file_name)
  python.assign("recalc",recalc)
  python.assign("use_ard",use_ard)
  python.exec("X,Kcc_ARD,varGPLVM_ARD = sclvm.fitGPLVM(idx=idx,k=k,standardize=standardize,out_dir=out_dir,file_name=file_name,recalc=recalc, use_ard=use_ard)")
  
  res=list()
  if(use_ard==TRUE){
    python.exec("Xard = list(varGPLVM_ARD['X_ARD'])")
    X_ard =  python.get("Xard")    
    res$X_ard = X_ard
  }
  
  python.exec("X = X.tolist()")
  X =  python.get("X")    
  res$X = X
  
  python.exec("Kcc = Kcc_ARD.tolist()")
  Kcc = do.call(rbind,python.get("Kcc"))
  res$Kcc = Kcc
  
  return(res)
}

varianceDecomposition_py = function(K=NULL,i0=1,i1=1){
  K = as.matrix(K)
  row.names(K)=c() #strip Y of column and row names
  colnames(K)=c()
  python.assign("K",K)
  python.exec("K = SP.array(K)") 
  
  python.assign("i0", as.integer(i0-1))
  python.assign("i1", as.integer(i1))
  
  python.exec("sclvm.varianceDecomposition(K=K,i0=i0,i1=i1)")
}


getVarianceComponents_py = function(normalize=normalize){
  python.assign("normalize", normalize)
  python.exec('var, var_info = sclvm.getVarianceComponents(normalize=normalize)')
  python.exec("conv = var_info['conv'].tolist()")
  gene_idx = python.get("var_info['gene_idx'].tolist()")
  geneID = python.get("geneID")
  res = list()
  res$conv = python.get("conv")
  res$var = do.call(rbind,python.get("var.tolist()"))
  colnames(res$var) = python.get("var_info['col_header'].tolist()")
  row.names(res$var) = geneID[gene_idx+1]
  python.assign("geneID_vd",geneID[gene_idx+1])
  
  return(res)
}


getCorrectedExpression_py = function(){
  python.exec("Ycorr = sclvm.getCorrectedExpression()")
  Ycorr = do.call(rbind, python.get("Ycorr.tolist()"))
  colnames(Ycorr) <- python.get("geneID_vd")
  return(Ycorr)
}

fitLMM_py <- function(K=NULL,i0=i0,i1=i1,verbose=TRUE, geneID=NULL){

  python.assign("i0", as.integer(i0-1))
  python.assign("i1", as.integer(i1))
  python.assign("verbose", verbose)
  python.assign("K",K)
  if(!is.null(K))python.exec("K = SP.array(K)") 
  
  python.exec("pv,beta,info = sclvm.fitLMM(K=K,i0=i0,i1=i1,verbose=verbose)")
  python.exec("beta[SP.where(SP.isnan(beta))]=0")#nans can happen for K=None and they are trouble in JSON
  python.exec("pv[SP.where(SP.isnan(pv))]=1")
  pv = do.call(rbind,python.get("pv.tolist()"))
  beta = do.call(rbind,python.get("beta.tolist()"))
  gene_idx_row = python.get("info['gene_idx_row'].tolist()")+1
  
  row.names(beta) = geneID[gene_idx_row]
  row.names(pv) = geneID[gene_idx_row]
  colnames(beta) = geneID
  colnames(pv) = geneID
  
  res = list()
  res$beta = beta
  res$pv = pv
  res$gene_idx = gene_idx_row
  return(res)
}
