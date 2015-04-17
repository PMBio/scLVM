fitTechnicalNoise <- function(nCountsEndo,nCountsERCC=NULL,
                              use_ERCC = TRUE,fit_type = "counts",plot=TRUE, fit_opts=NULL){  
  if(is.null(nCountsERCC) & use_ERCC==TRUE){
    print("You didn't provide ERCC counts so I will set use_ERCC to FALSE")
    use_ERCC = FALSE  
  }
  
  if( use_ERCC==FALSE &  (fit_type %in% c('counts','logvar'))){
    warning("Without ERCCs 'fit_type' 'log' is recommedned")
  }
  
  if((fit_type %in% c('counts', 'log','logvar'))==F){stop("'fit_type' needs to be 'counts', 'log' or 'logvar'")}
  
  if(fit_type=="count" & use_ERCC==FALSE){
    print("Without ERCCs fit needs to be perfromed in log-space")
    use_ERCC = FALSE  
  }
  
  if(use_ERCC==TRUE){
    if(fit_type=="counts"){
      meansEndo <- rowMeans( nCountsEndo )
      varsEndo <- rowVars( nCountsEndo )
      cv2Endo <- varsEndo / meansEndo^2
      
      meansERCC <- rowMeans( nCountsERCC )
      varsERCC <- rowVars( nCountsERCC )
      cv2ERCC <- varsERCC / meansERCC^2
      
      #Do fitting of technical noise
      if(!is.null(fit_opts)){
        if("mincv2" %in% names(fit_opts)){mincv2 = fit_opts$mincv2}else{mincv2=.3}
        if("quan" %in% names(fit_opts)){quan = fit_opts$quan}else{quan=0.8}
      }else{
        mincv2 = 0.3
        quan = 0.8
      }
      
      #normalised counts (with size factor)
      minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > mincv2 ) ], quan ) )
      useForFitA <- meansERCC >= minMeanForFitA
      fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
                          cv2ERCC[useForFitA] )
      
      #4. Transform to log-space and propagate error
      eps=1
      LogNcountsEndo=log10(nCountsEndo+eps)
      dLogNcountsEndo=1/((meansEndo+eps)*log(10))
      var_techEndo=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansEndo)*meansEndo^2
      LogVar_techEndo=(dLogNcountsEndo*sqrt(var_techEndo))^2 #error propagation 
      
      if(plot==TRUE){
        #plot fit
        par(mfrow=c(1,2))
        plot( meansERCC, cv2ERCC, log="xy", col=1+2*useForFitA)
        xg <- 10^seq( -3, 5, length.out=100 )
        lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='blue' )
        segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
                  meansERCC[useForFitA], fitA$fitted.values, col="gray" )
        legend('bottomleft',c('Genes used for fit', 'Fit technical noise'),pch=c(1, NA),lty =c(NA,1),col=c('green','blue'),cex=0.8)
        title('Mean-CV2 fit using ERCCs')
        
        #plot fot with all genes
        plot( meansEndo, cv2Endo, log="xy", col=1, xlab = 'Means', ylab = 'CV2')
        points(meansERCC, cv2ERCC, col='blue', pch=15)
        xg <- 10^seq( -3, 5, length.out=100 )
        lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='blue',lwd=2 )
        legend('bottomleft',c('Endogenous genes','ERCCs', 'Fit technical noise'),pch=c(1,15, NA),lty =c(NA,NA,1),col=c('black','blue','blue'),cex=0.8)        
        title('Mean-CV2 relationship')
        par(mfrow=c(1,1))
      }
      res = list()
      res$fit = fitA
      res$techNoiseLog = LogVar_techEndo
      
    }else{#with ERCCs in log space
      if(fit_type=="log"){
        
        LCountsEndo <- log10(nCountsEndo+1)
        LmeansEndo <- rowMeans( LCountsEndo )
        LvarsEndo <- rowVars( LCountsEndo )
        Lcv2Endo <- LvarsEndo / LmeansEndo^2
        
        LCountsERCC = log10(nCountsERCC+1)
        LmeansERCC <- rowMeans( LCountsERCC )
        LvarsERCC <- rowVars( LCountsERCC )
        Lcv2ERCC <- LvarsERCC / LmeansERCC^2
        
        if(!is.null(fit_opts)){
          if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=2}
        }else{
          minmean = .5
        }
        
        LogNcountsList=list()
        useForFitL=LmeansERCC>minmean
        LogNcountsList$mean=LmeansERCC[useForFitL]
        LogNcountsList$cv2=Lcv2ERCC[useForFitL]
        fit_loglin=nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=20,k=1))
        LogVar_techEndo_logfit <- coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*LmeansEndo)*LmeansEndo^2
        
        if(plot==TRUE){
          plot( LmeansEndo, Lcv2Endo, log="y", col=1,ylim=c(1e-3,1e2),xlab='meansLogEndo',ylab='cv2LogEndo')
          xg <- seq( 0, 5.5, length.out=100 )
          lines( xg, coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*xg ),lwd=2,col='green' )
          points(LmeansERCC, Lcv2ERCC,col='blue',pch=15,cex=1.1)
          legend('topright',c('Endo','ERCC'),pch=c(1,1,15),col=c('black','blue'))
        }
        res = list()
        res$fit = fit_loglin
        res$techNoiseLog = LogVar_techEndo_logfit
      }else{#with ERCCs fit variance in log space with loess
        LCountsEndo <- log10(nCountsEndo+1)
        LmeansEndo <- rowMeans( LCountsEndo )
        LvarsEndo <- rowVars( LCountsEndo )
        Lcv2Endo <- LvarsEndo / LmeansEndo^2
        
        LCountsERCC = log10(nCountsERCC+1)
        LmeansERCC <- rowMeans( LCountsERCC )
        LvarsERCC <- rowVars( LCountsERCC )
        
        if("span" %in% names(fit_opts)){span = fit_opts$span}else{span=0.8}
        if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=0.5}
        
        useForFitA <- LmeansERCC >= minmean
        fit_var2 = loess(LvarsERCC[useForFitA] ~ LmeansERCC[useForFitA], span=span, control=loess.control(surface="direct"))
        xg <- seq( 0, 5.5, length.out=100 )
        Var_techEndo_logfit_loess <-  predict(fit_var2, xg)
        
        minVar_ERCC = min(LvarsERCC[LmeansERCC>3])
        
        if(any(xg>3 & (Var_techEndo_logfit_loess<0.8*minVar_ERCC))){
          idx_1 = which(xg>3 & (Var_techEndo_logfit_loess<0.8*minVar_ERCC))[1]
          idx_end = length(Var_techEndo_logfit_loess)
          Var_techEndo_logfit_loess[idx_1:idx_end] = 0.8*minVar_ERCC        
        }
        
        if(plot==TRUE){
          plot( LmeansEndo, LvarsEndo, col=1,ylim=c(1e-3,150.5),log="y",xlab='meansLogEndo',ylab='VarLogEndo')
          points(LmeansERCC, LvarsERCC,col='blue',pch=15,cex=1.1)
          lines(xg, Var_techEndo_logfit_loess,lwd=3,col='blue',lty=1)  
          legend('topright',c('Endo. genes','ERCC', 'Tech. noise fit'),pch=c(1,15,NA), lty = c(NA,NA,1),col=c('black','blue', 'blue'))
        }
        
        #use model for endogenous genes
        xg=LmeansEndo
        Var_techEndo_logfit_loess <-  predict(fit_var2, xg)      
        
        if(any(xg>3 & Var_techEndo_logfit_loess<0.8*minVar_ERCC)){
          idx_1 = which(xg>3 & Var_techEndo_logfit_loess<0.8*minVar_ERCC)[1]
          idx_end = length(Var_techEndo_logfit_loess)
          Var_techEndo_logfit_loess[idx_1:idx_end] = 0.8*minVar_ERCC       
        }          
        
        res = list()
        res$fit = fit_var2
        res$techNoiseLog = Var_techEndo_logfit_loess
        
      }
      
    }
  }else{#no ERCCs available
    if(fit_type=="log"){
      LCountsEndo <- log10(nCountsEndo+1)
      LmeansEndo <- rowMeans( LCountsEndo )
      LvarsEndo <- rowVars( LCountsEndo )
      Lcv2Endo <- LvarsEndo / LmeansEndo^2
      
      if(!is.null(fit_opts)){
        if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=2}
      }else{
        minmean = 0.5
      }
      
      LogNcountsList = list()
      useForFitL = LmeansEndo>0.3
      LogNcountsList$mean = LmeansEndo[useForFitL]
      LogNcountsList$cv2 = Lcv2Endo[useForFitL]
      fit_loglin = nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=10,k=2))
      LogVar_techEndo_logfit <- coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*LmeansEndo)*LmeansEndo^2
      
      if(plot==TRUE){
        plot( LmeansEndo, Lcv2Endo, log="y", col=1,ylim=c(1e-3,1e2),xlab='meansLogEndo',ylab='cv2LogEndo')
        xg <- seq( 0, 5.5, length.out=100 )
        lines( xg, coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*xg ),lwd=2,col='green' )
      }
      
      res = list()
      res$fit = fit_loglin
      res$techNoiseLog = LogVar_techEndo_logfit
      
    }
    if(fit_type=='counts'){
      meansEndo <- rowMeans( nCountsEndo )
      varsEndo <- rowVars( nCountsEndo )
      cv2Endo <- varsEndo / meansEndo^2
      
      #Do fitting of technical noise
      if(!is.null(fit_opts)){
        if("mincv2" %in% names(fit_opts)){mincv2 = fit_opts$mincv2}else{mincv2=.3}
        if("quan" %in% names(fit_opts)){quan = fit_opts$quan}else{quan=0.8}
      }else{
        mincv2 = 0.3
        quan=0.8
      }
      
      #normalised counts (with size factor)
      minMeanForFitA <- unname( quantile( meansEndo[ which( cv2Endo > mincv2 ) ], quan ) )
      useForFitA <- meansEndo >= minMeanForFitA
      fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansEndo[useForFitA] ),
                          cv2Endo[useForFitA] )
      
      #4. Transform to log-space and propagate error
      eps=1
      LogNcountsEndo=log10(nCountsEndo+eps)
      dLogNcountsEndo=1/((meansEndo+eps)*log(10))
      var_techEndo=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansEndo)*meansEndo^2
      LogVar_techEndo=(dLogNcountsEndo*sqrt(var_techEndo))^2 #error propagation 
      
      if(plot==TRUE){
        #plot fit        
        plot( meansEndo, cv2Endo, log="xy", col=1+2*useForFitA, xlab = 'Means', ylab = 'CV2')
        xg <- 10^seq( -3, 5, length.out=100 )
        lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='blue' )
        legend('bottomleft',c('Genes used for fit', 'Fit baseline variation'),pch=c(1, NA),lty =c(NA,1),col=c('green','blue'),cex=0.8)
        title('Mean-CV2 fit using endogeneous genes')
      }
      res = list()
      res$fit = fitA
      res$techNoiseLog = LogVar_techEndo
      
      
    }
    if(fit_type=='logvar'){
      LCountsEndo <- log10(nCountsEndo+1)
      LmeansEndo <- rowMeans( LCountsEndo )
      LvarsEndo <- rowVars( LCountsEndo )
      Lcv2Endo <- LvarsEndo / LmeansEndo^2
      
      
      
      if("span" %in% names(fit_opts)){span = fit_opts$span}else{span=0.8}
      if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=0.5}
      
      useForFitA <- LmeansEndo >= minmean
      fit_var2 = loess(LvarsEndo[useForFitA] ~ LmeansEndo[useForFitA], span=span, control=loess.control(surface="direct"))
      xg <- seq( 0, 5.5, length.out=100 )
      Var_techEndo_logfit_loess <-  predict(fit_var2, xg)
      
      minVar_ERCC = min(LvarsEndo[LmeansEndo>3])
      
      if(any(xg>2.5 & (Var_techEndo_logfit_loess<0.6*minVar_ERCC))){
        idx_1 = which(xg>2.5 & (Var_techEndo_logfit_loess<0.6*minVar_ERCC))[1]
        idx_end = length(Var_techEndo_logfit_loess)
        Var_techEndo_logfit_loess[idx_1:idx_end] = 0.6*minVar_ERCC        
      }
      
      if(plot==TRUE){
        plot( LmeansEndo, LvarsEndo, col=1,ylim=c(1e-3,150.5),log="y",xlab='meansLogEndo',ylab='VarLogEndo')
        lines(xg, Var_techEndo_logfit_loess,lwd=3,col='blue',lty=1)  
        legend('topright',c('Endo. genes', 'Tech. noise fit'),pch=c(1,NA), lty = c(NA,1),col=c('black', 'blue'))
      }
      
      #use model for endogenous genes
      xg=LmeansEndo
      Var_techEndo_logfit_loess <-  predict(fit_var2, xg)      
      
      if(any(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_ERCC)){
        idx_1 = which(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_ERCC)[1]
        idx_end = length(Var_techEndo_logfit_loess)
        Var_techEndo_logfit_loess[idx_1:idx_end] = 0.6*minVar_ERCC       
      }          
      
      res = list()
      res$fit = fit_var2
      res$techNoiseLog = Var_techEndo_logfit_loess
      
    }
    
  }
  res    
}



getVariableGenes <- function(nCountsEndo, fit, method = "fit", threshold = 0.1, fit_type=NULL,sfEndo=NULL, sfERCC=NULL, plot=T){
  if(!(method %in% c("fdr","fit"))){
    stop("'method' needs to be either 'fdr' or 'fit'")
  }
  
  
  if(is.null(fit_type)){
    print("No 'fit_type' specified. Trying to guess its from parameter names")
    if("a0" %in% names(coefficients(fit)) & "a1tilde" %in% names(coefficients(fit))){fit_type="counts"}else{
      if("a" %in% names(coefficients(fit)) & "k" %in% names(coefficients(fit))){fit_type="log"}else{
        if(is.call(fit$call)){fit_type="logvar"}    
      }
    }
    print(paste("Assuming 'fit_type' is ","'",fit_type,"'",sep=""))
  }
  
  
  if(is.null(fit_type)){stop("Couldn't guess fit_type. Please specify it or run the getTechincalNoise 
                           function to obtain the fit")}
  
  
  if(!(fit_type %in% c("counts","log", "logvar")) & !is.null(fit_type)){
    stop("'fit_type' needs to be either 'fdr' or 'fit'")
  }
  
  
  if(method=='fdr' & fit_type!="counts"){stop("method='fdr', can only be used with fit_type 'counts'")}
  if(method=='fdr' & (is.null(sfERCC) | is.null(sfEndo))){stop("Please specify sfERCC and sfEndo when using method='fdr'")}
  
  
  if(method=='fdr'){
    meansEndo <- rowMeans( nCountsEndo )
    varsEndo <- rowVars( nCountsEndo )
    cv2Endo <- varsEndo / meansEndo^2

      
    

    minBiolDisp <- .5^2
    xi <- mean( 1 / sfERCC )
    m <- ncol(nCountsEndo)
    psia1thetaA <- mean( 1 / sfERCC ) +
      ( coefficients(fit)["a1tilde"] - xi ) * mean( sfERCC / sfEndo )
    cv2thA <- coefficients(fit)["a0"] + minBiolDisp + coefficients(fit)["a0"] * minBiolDisp
    testDenomA <- ( meansEndo * psia1thetaA + meansEndo^2 * cv2thA ) / ( 1 + cv2thA/m )
    
    pA <- 1 - pchisq( varsEndo * (m-1) / testDenomA, m-1 )
    padjA <- p.adjust( pA, "BH" )
   print( table( padjA < .1 ))
    is_het =  padjA < threshold
    is_het[is.na(is_het)] = FALSE
    
    if(plot==TRUE){
      plot( meansEndo, cv2Endo, log="xy", col=1+is_het,ylim=c(0.1,95), xlab='Mean Counts', ylab='CV2 Counts')
      xg <- 10^seq( -3, 5, length.out=100 )
      lines( xg, coefficients(fit)[1] + coefficients(fit)[2]/xg,lwd=2,col='green' )      
      legend('bottomright',c('Endo. genes','Var. genes',"Fit"),pch=c(1,1,NA),lty = c(NA,NA,1),col=c('black','red', 'green'),cex=0.7)   
    }
    
  }
  if(method=='fit' & fit_type=='log'){
    LCountsEndo <- log10(nCountsEndo+1)
    LmeansEndo <- rowMeans( LCountsEndo )
    Lcv2Endo = rowVars(LCountsEndo)/LmeansEndo^2
    is_het = (coefficients(fit)["a"] *10^(-coefficients(fit)["k"]*LmeansEndo) < Lcv2Endo) &  LmeansEndo>0.3  
    
    if(plot==TRUE){
      plot( LmeansEndo, Lcv2Endo, log="y", col=1+is_het,ylim=c(1e-3,1e2),xlab='meansLogEndo',ylab='cv2LogEndo')
      xg <- seq( 0, 5.5, length.out=100 )
      lines( xg, coefficients(fit)[1] *10^(-coefficients(fit)[2]*xg ),lwd=2,col='green' )
      legend('bottomright',c('Endo. genes','Var. genes',"Fit"),pch=c(1,1,NA),lty = c(NA,NA,1),col=c('black','red', 'blue'),cex=0.7)   
      
    }
    
  }
  
  if(method=='fit' & fit_type=='counts'){
    meansEndo <- rowMeans( nCountsEndo )
    varsEndo <- rowVars( nCountsEndo )
    cv2Endo <- varsEndo/meansEndo^2
    is_het = (coefficients(fit)[[1]] + coefficients(fit)[[2]]/meansEndo) < cv2Endo &  meansEndo>2
    
    if(plot==TRUE){
      plot( meansEndo, cv2Endo, log="xy", col=1+is_het,ylim=c(0.1,95), xlab='Mean Counts', ylab='CV2 Counts')
      xg <- 10^seq( -3, 5, length.out=100 )
      lines( xg, coefficients(fit)[1] + coefficients(fit)[2]/xg,lwd=2,col='green' )      
      legend('bottomright',c('Endo. genes','Var. genes',"Fit"),pch=c(1,1,NA),lty = c(NA,NA,1),col=c('black','red', 'green'),cex=0.7)   
    }
    
    
  }
  
  if(method=='fit' & fit_type=='logvar'){
    LCountsEndo <- log10(nCountsEndo+1)
    LmeansEndo <- rowMeans( LCountsEndo )
    LVarsEndo <- rowVars( LCountsEndo )
    
    xg = LmeansEndo
    
    Var_techEndo_logfit_loess =  predict(fit, LmeansEndo)
    
    minVar_Endo = min(LVarsEndo[LmeansEndo>2.5])
    
    if(any(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)){
      idx = which(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)
      Var_techEndo_logfit_loess[idx] = 0.6*minVar_Endo       
    }      
    
    is_het = (Var_techEndo_logfit_loess < LVarsEndo) &  LmeansEndo>0.3
    print(sum(is_het))
    
    
    if(plot==TRUE){
      plot( LmeansEndo, LVarsEndo, log="y", col=1+is_het,,xlab='meansLogEndo',ylab='varsLogEndo')
      xg <- seq( 0, 5.5, length.out=100 )
      Var_techEndo_logfit_loess =  predict(fit, xg)
      if(any(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)){
        idx_1 = which(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)[1]
        idx_end = length(Var_techEndo_logfit_loess)
        Var_techEndo_logfit_loess[idx_1:idx_end] = 0.6*minVar_Endo       
      }      
      
      lines( xg, Var_techEndo_logfit_loess,lwd=2,col='green' )
      legend('bottomright',c('Endo. genes','Var. genes',"Fit"),pch=c(1,1,NA),lty = c(NA,NA,1),col=c('black','red', 'green'),cex=0.7)   
      
      #print("Plotting not available for this method")}
    }

  }
  
  is_het
}


getEnsembl <- function(term, species = 'mMus'){
  if(!(species %in%c('mMus','Hs'))){stop("'species' needs to be either 'mMus' or 'Hs'")}
  
  if(species=='mMus'){
    if(require(org.Mm.eg.db)){
    xxGO <- AnnotationDbi::as.list(org.Mm.egGO2EG)}else{
      stop("Install org.Mm.eg.db package for retrieving gene lists from GO")
    }
  }else{
    if(require(org.Hs.eg.db)){
  xxGO <- AnnotationDbi::as.list(org.Hs.egGO2EG)}else{
    stop("Install org.Hs.eg.db package for retrieving gene lists from GO")
  }
  }
  cell_cycleEG <-unlist(xxGO[term])
  #get ENSEMBLE ids
  x <- org.Mm.egENSEMBL
  mapped_genes <- mappedkeys(x)
  xxE <- as.list(x[mapped_genes])
  ens_ids_cc<-unlist(xxE[cell_cycleEG])  
  
  ens_ids_cc
}


getSymbols <- function(ensIds, species = 'mMus'){
  if(!(species %in%c('mMus','Hs'))){stop("'species' needs to be either 'mMus' or 'Hs'")}
  
  if(species=='mMus'){
    require(org.Mm.eg.db)
  x <- org.Mm.egSYMBOL
  xxenseg <- AnnotationDbi::as.list(org.Mm.egENSEMBL2EG)}else{
    require(org.Hs.eg.db)
    x <- org.Hs.egSYMBOL
    xxenseg <- AnnotationDbi::as.list(org.Hs.egENSEMBL2EG)        
  }
  # Get the gene symbol that are mapped to an entrez gene identifiers
  gene_names = ensIds
  
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  xx <- as.list(x[mapped_genes])  
  gene_syms=unlist(xx[unlist(xxenseg[gene_names])])
  gene_names_list<-(lapply(xxenseg[gene_names],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
  sym_names=unlist(lapply(xx[unlist(gene_names_list)],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
  sym_names[is.na(sym_names)]=gene_names[is.na(sym_names)]
  
  sym_names
}

## Function for arranging ggplots
## Code adapted from Stephen Turner: https://gist.github.com/stephenturner/3724991#file-arrange_ggplot2-r

vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  require(grid)
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}

  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row	
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}


configLimix <- function(limix_path){
  python.assign('limix_path', limix_path)
  python.exec("sys.path.append(limix_path)")
  python.load(system.file("py","init_data.py",package="scLVM"))  
}
