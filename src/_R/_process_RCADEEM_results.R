library(gplots)
library(glmnet)
library(pROC)
library(data.table)

############ These functions implement non-negative logistic regression

# Sigmoid function
sigmoid <- function(z) 1/(1+exp(-z))

# The logistic cost Function. Note that except for theta[1], which represents the intercept,
# the rest of the values are actually log(b), where b is the coefficient. This enforces 'b'
# to be positive.
logistic_cost <- function( theta, X, Y )
{
  m <- nrow(X)
  g <- sigmoid( X %*% exp(theta[2:length(theta)]) + theta[1] )
  J <- (1/m)*sum((-Y*log(g)) - ((1-Y)*log(1-g)))
  return(J)
}

# Predict the probability based on a fit logistic function
logistic_predict <- function( theta, X )
{
  g <- sigmoid( X %*% theta[2:length(theta)] + theta[1] )
  return(g)
}

# The actual function for non-negative logistic regression
nnlogistic <- function( X, Y )
{
  fit.glm <- glm(Y~X,family="binomial") # first fit a simple logistic
  initial_theta <- coef( fit.glm ) # store the coefficients
  ncoef <- length(initial_theta)
  initial_theta[ 2:ncoef ][ initial_theta[ 2:ncoef ] <= 0 ] <- 0.00001
  initial_theta[ 2:ncoef ] <- log( initial_theta[ 2:ncoef ] )
  
  theta_optim <- optim( par=initial_theta, fn=logistic_cost, X=X, Y=Y )
  coefs <- theta_optim$par
  coefs[ 2:ncoef ] <- exp( coefs[ 2:ncoef ] )
  return( coefs )
}


############ Here's where the actual action happens

# Read the arguments, including the job label
args <- commandArgs(trailingOnly = TRUE)
jobLabel <- args[1]

input <- paste0("./out/",jobLabel,"/results.PFM.scores.txt")


# Read the HMM scores
cat("Reading the PFM (HMM) scores  ...\n")
hmm_scores <- try( fread( input, sep="\t" ) )
if( class(hmm_scores)[1]=="try-error" ) quit()

ncolumns <- ncol( hmm_scores )
nmotifs <- ncolumns-2
cat(paste0(nmotifs," PFM scores were read for ", nrow(hmm_scores), " sequences, including ", sum(hmm_scores$Status), " positive sequences.\n"))
if( nmotifs <= 0 )
{
  cat("ERROR: No motifs found in the input\n")
  quit()
}

# Display a heatmap of the HMM scores
if( nmotifs > 1 )
{
  jpeg(file=paste0("./out/",jobLabel,"/graphs.PFM_scores.jpg"),width=500,height=800)
  heatmap.2( as.matrix( hmm_scores[ hmm_scores$Status==1, 3:ncolumns, with=F ] ), hclustfun = function(x) hclust(x,method =  "complete"), margins=c(12,2), Colv = F, dendrogram = "row", trace="none", key=T, breaks=seq(0,1,length.out=256),col=colorRampPalette( c(rgb(0.95,0.95,0.95),"yellow", "orange", "red","dark red"))(255), labRow = F, ylab = "Peaks", key.title = "", key.xlab = "Binding score", key.ylab = "", density.info="none", lhei = c(0.3,2) )
  new_device <- dev.off()
}

cat("Performing non-negative logistic regression for distinguishing positive from negative sequences ...\n")
# Perform non-negative logistic regression to create a predictor of DNA-binding
coefs <- nnlogistic( as.matrix( hmm_scores[,3:ncolumns,with=F] ), hmm_scores$Status )
write.table( coefs, paste0("./out/",jobLabel,"/logistic_regression.coefficients.txt"), sep="\t", quote=F, col.names = F )

# Display a ROC curve of the predictions
predicted <- logistic_predict( coefs, as.matrix( hmm_scores[,3:ncolumns,with=F] ) )
myroc <- roc( hmm_scores$Status, c(predicted ) )
par(mar=c(10,10,10,10))
jpeg(file=paste0("./out/",jobLabel,"/graphs.logistic.ROC.jpg"),width=300,height=300)
plot(myroc)
new_device <- dev.off()

cat("Calculating weighted PFM scores ...\n")
# Use the coefficients to create a weighted HMM score matrix
weighted_hmm_scores <- t( t( as.matrix( hmm_scores[ , 3:ncolumns, with=F ] ) ) * coefs[2:(nmotifs+1)] )

# Display a heatmap of the weighted scores
if( nmotifs > 1 )
{
  jpeg(file=paste0("./out/",jobLabel,"/graphs.weighted_PFM_scores.jpg"),width=500,height=800)
  heatmap.2( weighted_hmm_scores[ hmm_scores$Status==1, ]  + coefs[1]/nmotifs, hclustfun = function(x) hclust(x,method =  "complete"), margins=c(12,2), Colv = F, dendrogram = "row", trace="none", key=T, breaks=seq(0,3,length.out=256),col=colorRampPalette( c(rgb(0.95,0.95,0.95),"yellow", "orange", "red","dark red"))(255), labRow = F, ylab = "Peaks", key.title = "", key.xlab = "Binding score", key.ylab = "", density.info="none", lhei = c(0.3,2) )
  new_device <- dev.off()
}
write.table( cbind( hmm_scores[,1:2,with=F], weighted_hmm_scores), paste0("./out/",jobLabel,"/graphs.weighted_PFM_scores.txt"), sep="\t", quote=F, row.names = F )


cat("Calculating ZF binding scores ...\n")
# Use the motif names to figure out what are the stat and end ZFs of each motif
# and then create a matrix that represents the relationship between motifs and the ZFs
motif_names <- colnames( weighted_hmm_scores )
zfs <- matrix( rep(0,nmotifs*100), nrow = nmotifs  )
max_zf <- 0
for( i in 1:nmotifs )
{
  range <- as.numeric( strsplit( strsplit( strsplit( motif_names[i], ":" )[[1]][2], "\\|"  )[[1]][1], "-" )[[1]] )
  max_zf <- max( max_zf, range )
  zfs[ i, range[1]:range[2] ] <- 1/(range[2]-range[1]+1)
}

# Create a profile representing the weighted probability of binding of each sequence to each ZF
meta_pfm_profile <- weighted_hmm_scores %*% zfs[ , 1:max_zf ]
meta_pfm_profile <- meta_pfm_profile + coefs[1]/max_zf

# Display the profile in a heatmap
jpeg(file=paste0("./out/",jobLabel,"/graphs.ZF_binding_scores.jpg"),width=500,height=800)
heatmap.2( meta_pfm_profile[ hmm_scores$Status==1, ], hclustfun = function(x) hclust(x,method =  "complete"), margins=c(4,2), Colv = F, dendrogram = "row", trace="none", key=T, breaks=seq(0, 0.6, length.out=256),col=colorRampPalette( c(rgb(0.95,0.95,0.95),"yellow", "orange", "red","dark red"))(255), labRow = F, xlab = "Zinc finger", ylab = "Peaks", key.title = "", key.xlab = "Binding score", key.ylab = "", density.info="none", lhei = c(0.3,2) )
write.table( cbind( hmm_scores[,1:2,with=F], meta_pfm_profile), paste0("./out/",jobLabel,"/graphs.ZF_binding_scores.txt"), sep="\t", quote=F, row.names = F )
new_device <- dev.off()
cat("Done!\n")
