library(rrLogit)
data(cig2019)

# get race/ethnicity variable, shorten the category names
# for better displays
y <- cig2019$hispallp_a
levels(y) <- c( "Hispanic", "NH White", "NH Black", "NH Asian",
 "NH AIAN", "NH AIAN+", "Other" )

# display frequencies and percentages for the true response
yFreqTable <- data.frame( table(y) )
yFreqTable$Pct <- round( 100*table(y)/length(y), 2 )
yFreqTable

# perturbation matrix described by Wang, Wu and Hu (2016)
pertMat <- rrPertMat( privLoss=1.5, nLevels=7, method="WWH" )
round(pertMat, 4)

# infuse with noise
set.seed(721)
ystar <- rrPerturbResponse( y, pertMat )

# display frequencies and percentages for the noisy response
ystarFreqTable <- data.frame( table(ystar) )
ystarFreqTable$Pct <- round( 100*table(ystar)/length(ystar), 2 )
ystarFreqTable

# method-of-moments estimates for true percentages
ystarFreqTable$MOM <- round( 100 * solve(pertMat) %*%
   table(ystar)/length(ystar), 2 )

# ML estimates for true percentages
fit <- rrLogit(ystar ~ 1, pertMat=pertMat)
fit <- rrLogit(fit, method="NR")
ystarFreqTable$ML <- round( 100 * as.vector(fitted(fit)), 2 )

# add prior information
fit <- rrLogit(ystar ~ 1, pertMat=pertMat,
   prior="DAP", priorFreqTot=10, priorAlloc=rep(1/7,7) )
   
# put y and ystar into data frame, remove incomplete observations
metro <- cig2019$urbrrl
region <- cig2019$region
age <- cig2019$agep_a
catAge <- cut( age, breaks=c(18,25,35,45,55,65,75,85), include.lowest=TRUE)
sex <- cig2019$sex_a
dF <- data.frame(y, ystar, metro, region, catAge, sex)
# identify rows with missing values and eliminate them
incomplete <- apply( is.na(dF), 1, any )
dF <- dF[ ! incomplete, ]
# see how many rows are left
NROW(dF)

# fit some models
fitA <- rrLogit( ystar ~ catAge + sex,
   pertMat=pertMat, data=dF,
   prior="DAP", priorFreqTot=7, priorAlloc=rep(1/7,7) )
fitA <- rrLogit(fitA, method="NR")
round( coef(fitA), 5 )
round( head( fitted(fitA) ), 6 )

fitB <- rrLogit( ystar ~ catAge + sex + region,
   pertMat=pertMat, data=dF,
   prior="DAP", priorFreqTot=7, priorAlloc=rep(1/7,7) )
fitB <- rrLogit(fitB, method="NR")
round( coef(fitB), 5 )
round( head( fitted(fitB) ), 6 )
