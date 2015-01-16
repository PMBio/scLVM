library(statmod)
library(gdata)
library(genefilter)
library(EBImage)
library(rhdf5)
library(DESeq)
library(statmod)
library(hom.Hs.inp.db)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(rPython)


#1. REad in data
load('./data_Tcells.Rdata')

dataMouse[ 1:10, 1:5 ]

geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
  substr( rownames(dataMouse), 1, 4 ) ] )

#2. calculate normalisation for counts
countsMmus <- dataMouse[ which( geneTypes=="ENSM" ), ]
countsERCC <- dataMouse[ which( geneTypes=="ERCC" ), ]
lengthsMmus <- dataMouse[ which( geneTypes=="ENSM" ), 1 ]
lengthsERCC <- dataMouse[ which( geneTypes=="ERCC" ), 1 ]


sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
#we do not want to normalize for differences in biological starting material... 
#as this can be due to different stages in cell cycle
sfMmus <- sfERCC 

nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsMmus <- t( t(countsMmus) / sfMmus )


meansERCC <- rowMeans( nCountsERCC )
varsERCC <- rowVars( nCountsERCC )
cv2ERCC <- varsERCC / meansERCC^2

#normalized counts (brennecke)
meansMmus <- rowMeans( nCountsMmus )
varsMmus <- rowVars( nCountsMmus )
cv2Mmus <- varsMmus / meansMmus^2

meansERCC <- rowMeans( nCountsERCC )
varsERCC <- rowVars( nCountsERCC )
cv2ERCC <- varsERCC / meansERCC^2

#Do fitting of technical noise

#normalised counts (with size factor)
minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .8 ) )
useForFitA <- meansERCC >= minMeanForFitA
fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
                    cv2ERCC[useForFitA] )

#plot fit
plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA)
xg <- 10^seq( -3, 5, length.out=100 )
lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg )
segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
          meansERCC[useForFitA], fitA$fitted.values, col="gray" )

#perfrom statistical test
minBiolDisp <- .5^2
xi <- mean( 1 / sfERCC )
m <- ncol(countsMmus)
psia1thetaA <- mean( 1 / sfERCC ) +
  ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfERCC / sfMmus )
cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
testDenomA <- ( meansMmus * psia1thetaA + meansMmus^2 * cv2thA ) / ( 1 + cv2thA/m )
pA <- 1 - pchisq( varsMmus * (m-1) / testDenomA, m-1 )
padjA <- p.adjust( pA, "BH" )
table( padjA < .1 )


#plot mean/cv2 relationship and 
plot( meansMmus, cv2Mmus, log="xy", col=1+(padjA<0.1),ylim=c(0.1,95), xlab='Mean Counts', ylab='CV2 Counts')
xg <- 10^seq( -3, 5, length.out=100 )
lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg,lwd=2,col='blue' )
points(meansERCC, cv2ERCC,col='blue',pch=15,cex=1.1)
#points(meansMmus[cc_gene_indices], cv2Mmus[cc_gene_indices],col=rgb(0,255,0,100,maxColorValue=255),pch=2,cex=0.75)
#points(meansMmus[ccCBall_gene_indices], cv2Mmus[ccCBall_gene_indices],col=rgb(0,255,0,20,maxColorValue=255),pch=2,cex=0.8)
legend('bottomleft',c('T-cells (padj >= 0.1)','T-cells (padj<0.1)','ERCC','Cell Cycle genes'),pch=c(1,1,15),col=c('black','red','blue','green'),cex=0.7)

#4. Transform to log-space and propagate error
eps=1
LogNcountsMmus=log10(nCountsMmus+eps)
dLogNcountsMmus=1/(meansMmus+eps)
var_techMmus=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansMmus)*meansMmus^2
LogVar_techMmus=(dLogNcountsMmus*sqrt(var_techMmus))^2 #error propagation 


#gene names in the T-cell data.
gene_names=rownames(nCountsMmus)
gene_names_het=gene_names[which(padjA<0.1)]

#all Cycle base genes homologs (top 600 genes)
hu2musAll=inpIDMapper(dataCB[1:600,3],'HOMSA','MUSMU',srcIDType='ENSEMBL',destIDType='ENSEMBL')
ccCBall_gene_indices=match(unlist(hu2musAll),rownames(nCountsMmus))

#get cell cycle genes from GO 
xxGO <- as.list(org.Mm.egGO2EG)
cell_cycleEG <-unlist(xxGO['GO:0007049'])
#get ENSEMBLE ids
x <- org.Mm.egENSEMBL
mapped_genes <- mappedkeys(x)
xxE <- as.list(x[mapped_genes])
ens_ids_cc<-unlist(xxE[cell_cycleEG])
cc_gene_indices <- na.omit(match(ens_ids_cc, rownames(dataMouse)))


#ensemble IDs to gene symbols
x <- org.Mm.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
xxenseg <- as.list(org.Mm.egENSEMBL2EG)
gene_syms=unlist(xx[unlist(xxenseg[gene_names])])
gene_names_list<-(lapply(xxenseg[gene_names],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
sym_names=unlist(lapply(xx[unlist(gene_names_list)],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
sym_names[is.na(sym_names)]=gene_names[is.na(sym_names)]

sym_names_het=sym_names[which(padjA<0.1)] #gene symbols of variable genes


#rename a few variables...
cellcyclegenes <- ens_ids_cc
cellcyclegenes_filter <- cc_gene_indices
cell_names <- colnames(nCountsMmus)
Y <- nCountsMmus
genes_heterogen <- (padjA<0.1)*1
countsERCC_mat=as.matrix(countsERCC * 1)
countsMmus_mat = as.matrix(countsMmus * 1)

h5save(ccCBall_gene_indices,gene_names,sym_names,sym_names_het,cellcyclegenes_filter,cellcyclegenes,cell_names,nCountsMmus,genes_heterogen,LogVar_techMmus,LogNcountsMmus,countsMmus_mat,sfERCC,countsERCC_mat,file='data_Tcells_normCounts.h5f')



