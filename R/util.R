fitTechnicalNoise <- function(nCountsEndo,nCountsERCC=NULL,
                              use_ERCC = FALSE,fit_type = "counts",plot=TRUE, fit_opts=NULL){  
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
        quan=0.8
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
      
    }else{
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
  }
  res    
}



getVariableGenes <- function(nCountsEndo, fit, method = "fit", threshold = 0.1, fit_type=NULL,sfEndo=NULL, sfERCC=NULL){
  if(!(method %in% c("fdr","fit"))){
    stop("'method' needs to be either 'fdr' or 'fit'")
  }
  if(is.null(fit_type)){
    print("No 'fit_type' specified. Trying to guess it from parameter names")
    if("a0" %in% names(coefficients(fit)) & "a1tilde" %in% names(coefficients(fit))){fit_type="counts"}else{
      if("a" %in% names(coefficients(fit)) & "k" %in% names(coefficients(fit))){fit_type="log"}else{
        if(is.call(techNoiseLog$fit$call)){fit_type="logvar"}    
      }
    }
    print(paste("Assuming 'fit_type' was ","'",fit_type,"'",sep=""))
  }
  
  if(is.null(fit_type)){stop("Couldn't guess fit_type. Please specify it or run getTechincalNoise 
                           function to obtain the fit")}
  if(method=='fdr' & fit_type!="counts"){stop("method='fdr', can only be used with fit_type 'counts'")}
  if(method=='fdr' & (is.null(sfERCC) | is.null(sfEndo))){stop("Please specify sfERCC and sfEndo when using method='fdr'")}
  
  
  if(method=='fdr'){
    meansEndo <- rowMeans( nCountsEndo )
    varsEndo <- rowVars( nCountsEndo )
    cv2Endo <- varsEndo/meansEndo^2
    
    minBiolDisp <- .5^2
    xi <- mean( 1 / sfERCC )
    m <- ncol(nCountsEndo)
    psia1thetaA <- mean( 1 / sfERCC ) +
      ( coefficients(fit)["a1tilde"] - xi ) * mean( sfERCC / sfEndo )
    cv2thA <- coefficients(fit)["a0"] + minBiolDisp + coefficients(fit)["a0"] * minBiolDisp
    testDenomA <- ( meansEndo * psia1thetaA + meansEndo^2 * cv2thA ) / ( 1 + cv2thA/m )
    pA <- 1 - pchisq( varsEndo * (m-1) / testDenomA, m-1 )
    padjA <- p.adjust( pA, "BH" )
    is_het =  padjA < threshold
    is_het[is.na(is_het)] = FALSE
    
  }
  if(method=='fit' & fit_type=='log'){
    LCountsEndo <- log10(nCountsEndo+1)
    LmeansEndo <- rowMeans( LCountsEndo )
    Lcv2Endo = rowVars(LCountsEndo)/LmeansEndo^2
    is_het = (coefficients(fit)["a"] *10^(-coefficients(fit)["k"]*LmeansEndo) < Lcv2Endo) &  LmeansEndo>0.3  
  }
  
  if(method=='fit' & fit_type=='counts'){
    meansEndo <- rowMeans( nCountsEndo )
    varsEndo <- rowVars( nCountsEndo )
    cv2Endo <- varsEndo/meansEndo^2
    is_het = (coefficients(fit)[[1]] + coefficients(fit)[[2]]/meansEndo) < cv2Endo &  meansEndo>2
  }
  
  if(method=='fit' & fit_type=='logvar'){
    LCountsEndo <- log10(nCountsEndo+1)
    LmeansEndo <- rowMeans( LCountsEndo )
    LVarsEndo <- rowVars( LCountsEndo )
    is_het = predict(fit, LmeansEndo) < LVarsEndo &  LmeansEndo>0.3
  }
  
  is_het
}

