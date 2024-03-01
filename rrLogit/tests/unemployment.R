# tests using the data from survey with randomized response published
# by van den Hout and van der Heijden (2007)

library(rrLogit)

# enter data in wide format
dfWide <- data.frame( Yes.Yes=68, Yes.No=52, No.Yes=103, No.No=189 )
# create perturbation matrix
Tmat <- matrix( c( .8, .2, .2, .8 ), 2, 2 )
bigT <- kronecker(Tmat, Tmat)
rownames(bigT) <- colnames(bigT) <- c("Yes.Yes", "Yes.No", "No.Yes", "No.No")
# fit intercept-only model with DAP
fitModeWide <- rrLogit( cbind(Yes.Yes, Yes.No, No.Yes, No.No) ~ 1,
   data=dfWide, pertMat=bigT,
   prior="DAP", priorFreqTot=2, priorAlloc=rep(1/4,4) )
summary(fitModeWide)


# enter the data in narrow format with frequencies
Y1.Y2 <- as.factor( 1:4 )
levels(Y1.Y2) <- c("Yes.Yes", "Yes.No", "No.Yes", "No.No" )
Freq <- c(68, 52, 103, 189)
dfNarrow <- data.frame( Y1.Y2, Freq )
# fit the intercept-only model using the same DAP as before
fitModeNarrow <- rrLogit( Y1.Y2 ~ 1, data=dfNarrow, freq=Freq, pertMat=bigT,
   prior="DAP", priorFreqTot=2, priorAlloc=rep(1/4, 4) )
summary(fitModeNarrow)

# check for equality of coefficients, fitted values, and residuals
all.equal( coef(fitModeWide), coef(fitModeNarrow) )
all.equal( fitted(fitModeWide), fitted(fitModeNarrow) )
all.equal( residuals(fitModeWide), residuals(fitModeNarrow) )

# run saturated model on Wide data
fitModeWideSat <- rrLogit( cbind(Yes.Yes, Yes.No, No.Yes, No.No) ~ 1,
   data=dfWide, pertMat=bigT,
   prior="DAP", priorFreqTot=2, priorAlloc=rep(1/4,4),
   saturated=TRUE )
   summary(fitModeWideSat)
all.equal( fitted(fitModeWide), fitted(fitModeWideSat) )
all.equal( residuals(fitModeWide), residuals(fitModeWideSat) )

# run saturated model on narrow data
fitModeNarrowSat <- rrLogit( Y1.Y2 ~ 1, data=dfNarrow, freq=Freq, pertMat=bigT,
   prior="DAP", priorFreqTot=2, priorAlloc=rep(1/4, 4),
   saturated=TRUE )
all.equal( fitted(fitModeWide), fitted(fitModeNarrowSat) )
all.equal( residuals(fitModeWide), residuals(fitModeNarrowSat) )

# run MCMC on wide data, producing ten imputations
set.seed(5812)
fitMCMCWide <- rrLogit(fitModeWide, method="MCMC",
   control=list( tuneRWM=c(100,1), imputeEvery=500))
summary(fitMCMCWide)
impList(fitMCMCWide)

# run MCMC on narrow data, producing ten imputations
set.seed(2231)
fitMCMCNarrow <- rrLogit(fitModeNarrow, method="MCMC",
   control=list( tuneRWM=c(100,1), imputeEvery=500))
summary(fitMCMCNarrow)
impList(fitMCMCNarrow)

# run MCMC on wide data, saturated model
set.seed(598)
fitMCMCWideSat <- rrLogit(fitModeWideSat, method="MCMC",
   control=list( tuneRWM=c(100,1), imputeEvery=500) )
summary(fitMCMCWideSat)
impList(fitMCMCWideSat)


