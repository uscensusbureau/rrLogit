# all non-exported functions -- accessible to other package functions, but
# not invoked directly by users -- have names that begin with .

rrPertMat <- function( privLoss=Inf, nLevels=2L, method=c("WWH") ) {
   # double-constant perturbation matrix of Wang, Wu and Hu
   privLoss <- as.double( privLoss )[1L]
   stopifnot( privLoss > 0 )
   nLevels <- as.integer( nLevels )[1L]
   stopifnot( nLevels >= 2L )
   method <- match.arg(method)
   p <- 1 / ( 1 + ( nLevels - 1 ) * exp( -privLoss ) )
   pertMat <- matrix( numeric(), nLevels, nLevels )
   pertMat[ row(pertMat) == col(pertMat) ] <- p
   pertMat[ row(pertMat) != col(pertMat) ] <- ( 1 - p ) / ( nLevels - 1 )
   return( pertMat )
}

.isPertMat <- function( object, order, tol=1e-08 ) {
   if( ! is.matrix(object) ) return(FALSE)
   if( NROW(object) != NCOL(object) ) return(FALSE)
   if( NROW(object) != order ) return(FALSE)
   if( any( object < 0 ) || any( object > 1 ) ) return(FALSE)
   if( any( abs( apply( object, 2, sum ) - 1 ) > tol ) ) return(FALSE)
   return(TRUE)
}

rrPerturbResponse <- function(y, pertMat) {
   # Function for infusing noise into a factor (categorical variable)
   #        y : factor variable
   #  pertMat : perturbation matrix
   # Returns a noise-infused version of y
   stopifnot( is.factor(y) )
   stopifnot( .isPertMat( pertMat, nlevels(y) ) )
   ystar <- y
   for( j in 1:nlevels(y) ) {
      w <- ( y == levels(y)[j] )
      ystar[w] <- sample( levels(y), size=sum(w), prob=pertMat[,j],
         replace=TRUE )
   }
   return(ystar)
}

.unstackCoef <- function( coefVec, predNames, responseBaseLevels,
   baselineLev ) {
   p <- length(predNames)
   r <- length(responseBaseLevels)
   stopifnot( length(coefVec) == p*(r-1L) )
   stopifnot( baselineLev %in% responseBaseLevels )
   stopifnot( sum( match(responseBaseLevels, baselineLev,
      nomatch=0L) ) == 1L )
   result <- matrix( numeric(1L), p, r )
   rownames(result) <- predNames
   colnames(result) <- responseBaseLevels
   posn <- 0L
   for(k in 1:r) {
      if( responseBaseLevels[k] != baselineLev ) {
         for(j in 1:p) {
            posn <- posn + 1L
	    result[j,k] <- coefVec[posn]
         }
      }
   }
   return(result)
}

.fixRHS <- function(obj) {
   # replaces RHS of formula by a list of predictor variables
   # separated by '+', as expected for a saturated model
   form <- obj
   form[[2]] <- form[[3]]
   form[[3]] <- NULL
   form <- if( length( all.vars(form) ) == 0 ) as.formula( "~ 1" ) else
      as.formula( paste( "~", paste( all.vars(form), collapse=" + ") ) )
   obj[[3]] <- form[[2]]
   obj
}

.isWide <- function( obj ) {
   # checks formula, returns TRUE if wide format, FALSE if narrow
   if(  ( ! inherits(obj, "formula") ) ||
      ( length(obj) != 3L ) ) stop( gettext(
      "Argument 'obj' is not a two-sided formula" ), domain = NA ) 
   if( is.call( obj[[2L]] ) && identical( obj[[2L]][[1L]], quote(cbind) ) ) {
      if( length( obj[[2L]] ) < 3L ) stop( gettext(
      "Need at least two variables enclosed by 'cbind()'" ), domain = NA )
      return( TRUE )
   } else if( is.symbol( obj[[2L]] ) ) {
      return( FALSE )
   }
   stop( gettext(
      "Left-hand side of 'obj' must be a factor or a call to cbind()" ),
      domain = NA )
}

.processWideData <- function(obj, mF ) {
   fMat <- model.response(mF)
   # note: fMat could be integer or double, might be a long array
   if( ! is.numeric(fMat) ) stop( gettext(
      "Response-frequency matrix is non-numeric" ), domain = NA )
   if( any( fMat < 0 ) ) stop( gettext(
      "Response-frequency matrix contains negative values" ), domain = NA )
   if( any( is.na(fMat) ) ) stop( gettext(
      "Response-frequency matrix contains missing values" ), domain = NA )
   freqRowInputData <- apply( fMat, 1, sum )
   names( freqRowInputData ) <- NULL
   storage.mode(freqRowInputData) <- "double"
   weightRowInputData <- freqRowInputData
   # collapse to unique covariate patterns
   x <- model.matrix( obj, mF )
   if( any( is.na(x) ) ) stop( gettext(
      "Model matrix contains missing values" ), domain = NA )
   xdf <- as.data.frame(x)
   cx <- do.call( "paste", c( structure( xdf[,,drop=FALSE], names=NULL ),
      sep="\r" ) )
   uniqueRowNumbers.x <- (1:NROW(x))[ ! duplicated(cx) ]
   cxUnique <- cx[ uniqueRowNumbers.x ]
   nCovPatt <- length(cxUnique)   # will be integer
   covPatt <- match(cx, cxUnique) # will be integer
   names(covPatt) <- NULL
   modelMatrix <- x[uniqueRowNumbers.x,,drop=FALSE] # will be double
   # aggregate frequencies within unique covariate patterns
   aggForm <- obj
   aggForm[[3L]] <- quote(`(covPatt)`)
   dFtmp <- as.data.frame( fMat )
   dFtmp$`(covPatt)` <- covPatt
   aggF <- aggregate( aggForm, data=dFtmp, FUN=sum )
   aggF <- aggF[ order(aggF$`(covPatt)`), ]
   aggF$`(covPatt)` <- NULL
   aggF <- data.matrix(aggF)
   weightForCovPatt <- freqForCovPatt <- apply(aggF, 1, sum)
   names(weightForCovPatt) <- names(freqForCovPatt) <- NULL
   storage.mode(weightForCovPatt) <- storage.mode(freqForCovPatt) <- "double"
   # data patterns
   weightForDataPatt <- freqForDataPatt <- as.vector( t(aggF) )
   covPattForDataPatt <- as.vector( t( row(aggF) ) )
   responseForDataPatt <- as.vector( t( col(aggF) ) )
   names(weightForDataPatt) <- names(freqForDataPatt) <-
      names(covPattForDataPatt) <- names(responseForDataPatt) <- NULL
   nDataPatt <- length(weightForDataPatt)
   storage.mode(weightForDataPatt) <- storage.mode(freqForDataPatt) <- "double"
   # empty covariate patterns
   emptyCovPatt <- ( freqForCovPatt == 0 )
   # count parameters
   nParamThisModel <- NCOL(modelMatrix) * ( NCOL(fMat) - 1L )
   nParamSatModel <- sum( ! emptyCovPatt ) * ( NCOL(fMat) - 1L )
   # objects to return
   dimVec <- c( nrowInputData = NROW(mF), nLevels = NCOL(fMat),
      nCovPatt = nCovPatt, nDataPatt = nDataPatt,
      ncolMM = NCOL(modelMatrix),
      nParamThisModel = nParamThisModel, nParamSatModel = nParamSatModel,
      wideFormatInt = 1L )
   storage.mode(dimVec) <- storage.mode(covPatt) <-
      storage.mode(covPattForDataPatt) <-
      storage.mode(responseForDataPatt) <- "integer" 
   storage.mode(modelMatrix) <- storage.mode(freqRowInputData) <- 
      storage.mode(weightRowInputData) <-
      storage.mode(freqForCovPatt) <- storage.mode(weightForCovPatt) <-
      storage.mode(freqForDataPatt) <- storage.mode(weightForDataPatt) <-
          "double"
   list(
      dimVec = dimVec,
      rowNamesInputData = rownames(mF),
      responseLevels = colnames(fMat),
      freqRowInputData = freqRowInputData,
      weightRowInputData = weightRowInputData,
      covPatt = covPatt,
      modelMatrix = modelMatrix,
      freqForCovPatt = freqForCovPatt,
      weightForCovPatt = weightForCovPatt,
      freqForDataPatt = freqForDataPatt,
      weightForDataPatt = weightForDataPatt,
      covPattForDataPatt = covPattForDataPatt,
      responseForDataPatt = responseForDataPatt,
      freqSupplied = TRUE,
      weightSupplied = FALSE,
      wideFormat = TRUE,
      dataPatt = NULL, # b/c rows of input data hold multiple data patterns
      emptyCovPatt = emptyCovPatt,
      fMat = fMat,
      responseRowInputData = NULL, # b/c rows of input data have all responses
      microData = FALSE
      )
}

.processNarrowData <- function(obj, mF) {
   y <- model.response(mF)
   if( ! is.factor(y) ) stop( gettext(
      "Response variable is not a factor" ), domain = NA )
   nLevels <- nlevels(y)
   responseLevels <- levels(y)
   y <- unclass(y)
   # collapse to unique covariate patterns
   x <- model.matrix( obj, mF )
   if( any( is.na(x) ) ) stop( gettext(
      "Model matrix contains missing values" ), domain = NA )
   xdf <- as.data.frame(x)
   cx <- do.call( "paste", c( structure( xdf[,,drop=FALSE], names=NULL ),
      sep="\r" ) )
   uniqueRowNumbers.x <- (1:NROW(x))[ ! duplicated(cx) ]
   cxUnique <- cx[ uniqueRowNumbers.x ]
   nCovPatt <- length(cxUnique)
   covPatt <- match(cx, cxUnique)
   names(covPatt) <- NULL
   modelMatrix <- x[uniqueRowNumbers.x,,drop=FALSE]
   # identify unique data patterns
   cxy <- do.call( "paste", list( covPatt, y, sep="\r" ) )
   uniqueRowNumbers.xy <- (1:NROW(x))[ ! duplicated(cxy) ]
   cxyUnique <- cxy[ uniqueRowNumbers.xy ]
   o <- order( covPatt[ uniqueRowNumbers.xy ], y[ uniqueRowNumbers.xy ] )
   cxyUnique <- cxyUnique[o]
   nDataPatt <- length(cxyUnique)
   dataPatt <- match( cxy, cxyUnique )
   names(dataPatt) <- NULL
   covPattForDataPatt <- covPatt[ uniqueRowNumbers.xy ][o]
   responseForDataPatt <- y[ uniqueRowNumbers.xy ][o]
   names(covPattForDataPatt) <- names(responseForDataPatt) <- NULL
   # aggregate frequencies / weights within patterns
   weightSupplied <- freqSupplied <- FALSE
   freqRowInputData <- weightRowInputData <- rep( 1, NROW(mF) )
   microData <- TRUE
   if( ! is.null( mF$`(freq)` ) ) {
      freqSupplied <- TRUE
      microData <- FALSE
      if( ! is.numeric( mF$`(freq)` ) )  stop( gettext(
         "Supplied 'freq' is not numeric" ), domain = NA )
      if( any( is.na( mF$`(freq)` ) ) )  stop( gettext(
         "Missing values in 'freq' are not allowed" ), domain = NA )
      if( any( mF$`(freq)` < 0 ) )  stop( gettext(
         "Negative values in 'freq' are not allowed" ), domain = NA )
      weightRowInputData <- freqRowInputData <- mF$`(freq)`
      names(freqRowInputData) <- NULL
   }
   if( ! is.null( mF$`(weight)` ) ) {
      weightSupplied <- TRUE
      if( ! is.numeric( mF$`(weight)` ) )  stop( gettext(
         "Supplied 'weight' is not numeric" ), domain = NA )
      if( any( is.na( mF$`(weight)` ) ) )  stop( gettext(
         "Missing values in 'weight' are not allowed" ), domain = NA )
      if( any( mF$`(weight)` < 0 ) )  stop( gettext(
         "Negative values in 'weight' are not allowed" ), domain = NA )
      weightRowInputData <- mF$`(weight)`
      names(weightRowInputData) <- NULL
   }
   # make fMat
   fMat <- matrix( 0, NROW(mF), length(responseLevels) )
   storage.mode(fMat) <- storage.mode(freqRowInputData)
   for( j in 1:length(responseLevels) ) fMat[y==j,j] <- freqRowInputData[y==j]
   rownames(fMat) <- rownames(mF)
   colnames(fMat) <- responseLevels
   # 
   storage.mode(freqRowInputData) <-
      storage.mode(weightRowInputData) <- "double"
   dFtmp <- data.frame( freqRowInputData, weightRowInputData,
      covPatt, dataPatt )
   names(dFtmp) <- c( "freq", "weight", "covPatt", "dataPatt" )
   aggF <- aggregate( cbind(freq, weight) ~ covPatt, data=dFtmp, FUN=sum )
   aggF <- aggF[ order(aggF$covPatt), ]
   freqForCovPatt <- aggF$freq
   weightForCovPatt <- aggF$weight
   names(freqForCovPatt) <- names(weightForCovPatt) <- NULL
   storage.mode(weightForCovPatt) <- storage.mode(freqForCovPatt) <- "double"
   aggF <- aggregate( cbind(freq, weight) ~ dataPatt, data=dFtmp, FUN=sum )
   aggF <- aggF[ order(aggF$dataPatt), ]
   freqForDataPatt <- aggF$freq
   weightForDataPatt <- aggF$weight
   names(freqForDataPatt) <- names(weightForDataPatt) <- NULL
   storage.mode(weightForDataPatt) <- storage.mode(freqForDataPatt) <- "double"
   # empty covariate patterns
   emptyCovPatt <- ( weightForCovPatt == 0 )
   # count parameters
   nParamThisModel <- NCOL(modelMatrix) * ( nLevels - 1L )
   nParamSatModel <- sum( ! emptyCovPatt ) * ( nLevels - 1L )
   # objects to return
   dimVec <- c( nrowInputData = NROW(mF), nLevels = nLevels,
      nCovPatt = nCovPatt, nDataPatt = nDataPatt,
      ncolMM = NCOL(modelMatrix),
      nParamThisModel = nParamThisModel, nParamSatModel = nParamSatModel,
      wideFormatInt = 0L )
   storage.mode(dimVec) <- storage.mode(covPatt) <- storage.mode(dataPatt) <-
      storage.mode(covPattForDataPatt) <-
      storage.mode(responseForDataPatt) <- "integer" 
   storage.mode(modelMatrix) <- storage.mode(freqRowInputData) <- 
      storage.mode(weightRowInputData) <-
      storage.mode(freqForCovPatt) <- storage.mode(weightForCovPatt) <-
      storage.mode(freqForDataPatt) <- storage.mode(weightForDataPatt) <-
          "double"
   list(
      dimVec = dimVec,
      rowNamesInputData = rownames(mF),
      responseLevels = responseLevels,
      freqRowInputData = freqRowInputData,
      weightRowInputData = weightRowInputData,
      covPatt = covPatt,
      modelMatrix = modelMatrix,
      freqForCovPatt = freqForCovPatt,
      weightForCovPatt = weightForCovPatt,
      freqForDataPatt = freqForDataPatt,
      weightForDataPatt = weightForDataPatt,
      covPattForDataPatt = covPattForDataPatt,
      responseForDataPatt = responseForDataPatt,
      freqSupplied = freqSupplied,
      weightSupplied = weightSupplied,
      wideFormat = FALSE,
      dataPatt = dataPatt,
      emptyCovPatt = emptyCovPatt,
      fMat = fMat,
      responseRowInputData = y,
      microData = microData
   )
}

rrLogitControl <- function( iterMaxNR = 50L, iterMaxFS = 200L,
   iterMaxEM = 1000L, iterMaxMstep = 25L, critConverge = 1e-08,
   critBoundary = 1e-08, iterApproxBayes = 1L, imputeApproxBayes = FALSE,
   iterMCMC = 5000L, burnMCMC = 0L, thinMCMC = 1L, imputeEvery = 0L,
   typeMCMC = c("RWM","DA"), tuneDA = c(100,.2,.5),
   tuneRWM = c(1000,.2), stuckLimit = 25L, startValJitter = 0 ) {
   stopifnot( ( iterMaxNR <- as.integer(iterMaxNR)[1L] ) >= 0L )
   stopifnot( ( iterMaxFS <- as.integer(iterMaxFS)[1L] ) >= 0L )
   stopifnot( ( iterMaxEM <- as.integer(iterMaxEM)[1L] ) >= 0L )
   stopifnot( ( iterMaxMstep <- as.integer(iterMaxMstep)[1L] ) >= 0L )
   stopifnot( ( critConverge <- as.double(critConverge)[1L] ) > 0 )
   stopifnot( ( critBoundary <- as.double(critBoundary)[1L] ) > 0 )
   stopifnot( ( iterApproxBayes <- as.integer(iterApproxBayes)[1L] ) >= 0L )
   imputeApproxBayes <- as.logical( imputeApproxBayes )[1L]
   stopifnot( ( iterMCMC <- as.integer(iterMCMC)[1L] ) >= 0L )
   stopifnot( ( burnMCMC <- as.integer(burnMCMC)[1L] ) >= 0L )
   stopifnot( ( thinMCMC <- as.integer(thinMCMC)[1L] ) >= 1L )
   stopifnot( ( imputeEvery <- as.integer(imputeEvery)[1L] ) >= 0L )
   typeMCMC <- match.arg( typeMCMC )
   tuneDA <- as.numeric( tuneDA )[1:3]
   stopifnot( tuneDA[1L] > 0 )
   stopifnot( tuneDA[3L] > 0 )
   tuneRWM <- as.numeric(tuneRWM)[1:2]
   stopifnot( tuneRWM[1L] > 0 )
   stopifnot( tuneRWM[2L] > 0 )
   stopifnot( ( stuckLimit <- as.integer(stuckLimit)[1L] ) > 0 )
   stopifnot( ( startValJitter <- as.double(startValJitter)[1L] ) >= 0 )
   list(
      iterMaxNR = iterMaxNR,
      iterMaxFS = iterMaxFS,
      iterMaxEM = iterMaxEM,
      iterMaxMstep = iterMaxMstep,
      critConverge = critConverge,
      critBoundary = critBoundary,
      iterApproxBayes = iterApproxBayes,
      imputeApproxBayes = imputeApproxBayes,
      iterMCMC = iterMCMC,
      burnMCMC = burnMCMC,
      thinMCMC = thinMCMC,
      imputeEvery = imputeEvery,
      typeMCMC = typeMCMC,
      tuneDA = tuneDA,
      tuneRWM = tuneRWM,
      stuckLimit = stuckLimit,
      startValJitter = startValJitter )
}


rrLogit <- function( obj, ...) {
   # S3 generic function
   UseMethod("rrLogit")
}
   
rrLogit.default <- function( obj, ...) {
   stop( gettext(
      'First argument must be an object of class "formula" or "rrLogit"'),
      domain = NA )
}

.responseVarName <- function(obj) {
   stopifnot( inherits(obj, "formula")  )
   stopifnot( length(obj) > 2 )
   if( .isWide(obj) ) NULL else as.character( obj[[2]] )
}

rrLogit.formula <- function( obj, data, baseline=1L, freq, weight, 
   pertMat=NULL, privLoss=NULL, saturated=FALSE,
   prior = c("none", "DAP"), priorFreqTot = NULL, priorAlloc = NULL,
   method = c("EM", "NR", "FS", "MCMC", "approxBayes"), 
   startVal = NULL, control = list(), ...) { 
   #----------------------------------------------
   saturated <- as.logical( saturated )[1L]
   saturatedInt <- as.integer(saturated)
   if( saturated ) obj <- .fixRHS( obj )
   #--------------------------------------------
   # create model frame, process data
   wideFormat <- .isWide(obj)
   if( wideFormat ) {
      if( ! missing(freq) ) warning( gettext(
         "Data in wide format, argument 'freq' was ignored" ),
         domain = NA ) 
      if( ! missing(weight) ) warning( gettext(
         "Data in wide format, argument 'weight' was ignored" ),
	 domain = NA ) 
   } else {
      if( ( ! missing(freq) ) && ( ! missing(weight) ) ) stop( gettext(
         "You cannot supply both 'freq' and 'weight'" ), domain = NA ) 
   }
   mc <- match.call( expand.dots=FALSE )
   mc[[1L]] <- quote(stats::model.frame)
   mc$baseline <- mc$privLoss <- mc$pertMat <- mc$saturated <- mc$prior <-
      mc$priorFreqTot <- mc$priorAlloc <- mc$method <- mc$startVal <-
      mc$control <- NULL
   if( wideFormat ) {
      mc$freq <- mc$weight <- NULL
      m <- match( c("obj", "data"), names(mc), nomatch=0L )
   } else {
      m <- match( c("obj", "data", "freq", "weight"), names(mc), nomatch=0L )
   }
   mc <- mc[ c(1L,m) ]
   names(mc)[2] <- "formula"
   mc$na.action <- as.name("na.fail")  # maybe change this later
   mc$drop.unused.levels <- TRUE       # maybe change this later
   mF <- eval( mc, parent.frame() )
   pD <- if( wideFormat ) .processWideData(obj, mF) else
      .processNarrowData(obj, mF)
   surveyMode <- pD$weightSupplied
   if( surveyMode ) {
      weightScaleFac <-
         sum( pD$freqRowInputData ) / sum( pD$weightRowInputData )
      fitweightRowInputData <- pD$weightRowInputData * weightScaleFac
      fitweightCovPatt <- pD$weightForCovPatt * weightScaleFac
      fitweightDataPatt <- pD$weightForDataPatt * weightScaleFac
   } else {
      fitweightRowInputData <- pD$freqRowInputData
      fitweightCovPatt <- pD$freqForCovPatt
      fitweightDataPatt <- pD$freqForDataPatt
   }
   #--------------------------------------------
   method <- match.arg( method )
   methodInt <- match( method, c("EM", "NR", "FS", "MCMC", "approxBayes") )
   storage.mode(methodInt) <- "integer"
   if( ( ! ( method %in% c("EM","NR","FS") ) ) & surveyMode ) stop( gettextf(
      "Method '%s' not defined with survey weights", method), domain = NA )
   #--------------------------------------------
   # ensure modelMatrix has full rank
   if( ! saturated ) {
      if( qr(pD$modelMatrix)$rank < NCOL(pD$modelMatrix) ) stop( gettext(
         "modelMatrix does not have full rank" ), domain=NA )
      emptyCovPatt <- ( fitweightCovPatt == 0 )
      if( qr(pD$modelMatrix[(!emptyCovPatt),,drop=FALSE])$rank <
         NCOL(pD$modelMatrix) ) stop( gettext(
       "modelMatrix does not have full rank due to empty covariate patterns" ), 
         domain=NA )
   }      
   #--------------------------------------------
   prior <- match.arg( prior )
   priorFreqTotSupplied <- FALSE
   priorAllocSupplied <- FALSE
   if( prior == "none" ) {
      priorInt <- 1L
      priorFreqTot <- 0
      priorAlloc <- numeric( pD$dimVec["nLevels"] )
   } else if( prior == "DAP" ) {
      priorInt <- 2L
      if( is.null( priorFreqTot ) ) {
         priorFreqTot <- as.double( pD$dimVec["nParamThisModel"] )
      } else {
         priorFreqTot <- as.double( priorFreqTot )[1L]
         stopifnot( priorFreqTot >= 0 )
	 priorFreqTotSupplied <- TRUE
      }
      if( is.null( priorAlloc ) ) {
         priorAlloc <- numeric( pD$dimVec["nLevels"] )
      } else {
         stopifnot( is.numeric( priorAlloc ) )
         priorAlloc <- as.vector( priorAlloc )
	 stopifnot( length( priorAlloc ) == pD$dimVec["nLevels"] )
         priorAllocSupplied <- TRUE
      }
   }
   storage.mode(priorInt) <- "integer"
   storage.mode(priorAlloc) <- storage.mode(priorFreqTot) <- "double"
   names(priorAlloc) <- pD$responseLevels
   #--------------------------------------------
   if( is.numeric( baseline ) ) {
      baselineInt <- as.integer( baseline[[1L]] )
      if( ( baselineInt < 0 ) || ( baselineInt > pD$dimVec["nLevels"] ) ) 
         stop( gettext(
         "Integer provided for 'baseline' is out of range" ), domain = NA )
      baselineLev <- pD$responseLevels[ baselineInt ]
   } else if( is.character( baseline ) ) {
      baselineLev <- baseline[[1L]]
      if( ! ( baselineLev %in% pD$responseLevels ) ) stop( gettext(
        "Character string provided for 'baseline' is not a response category"),
         domain = NA )
      baselineInt <- match( baselineLev, pD$responseLevels )
   } else stop( gettext(
      "Argument 'baseline' must be an integer or character string" ),
      domain = NA )
   #--------------------------------------------
   if( ( ! is.null( pertMat ) ) & ( ! is.null( privLoss ) ) ) stop( gettext(
      "You cannot specify both 'pertMat' and 'privLoss'" ), domain = NA )
   if( ! is.null( pertMat ) ) {
      if( ! .isPertMat( pertMat, pD$dimVec["nLevels"] ) ) stop( gettextf(
         "'pertMat' is not a perturbation matrix of order %i",
          pD$dimVec["nLevels"] ), domain = NA )
   } else if( ! is.null( privLoss ) ) {
      privLoss <- as.double( privLoss )[1L]
      stopifnot( privLoss >= 0 )
      pertMat <- rrPertMat( privLoss, pD$dimVec["nLevels"] )
   } else {
      pertMat <- diag( pD$dimVec["nLevels"] )
   }
   rownames(pertMat) <- colnames(pertMat) <- pD$responseLevels
   if( qr( pertMat )$rank < NROW(pertMat) ) stop( gettext(
      "Perturbation matrix 'pertMat' is singular" ), domain = NA )
   pertMatInv <- solve(pertMat)  
   storage.mode(pertMatInv) <- storage.mode(pertMat) <- "double"
   #--------------------------------------------
   ident <- all( diag( pD$dimVec["nLevels"] ) == pertMat )
   if( ident & ( method %in% c("EM","FS") ) & ( !saturated ) ) {
      method <- "NR"
      methodInt <- match( method, c("EM", "NR", "FS", "MCMC", "approxBayes") )
      storage.mode(methodInt) <- "integer"
   }
   #--------------------------------------------
   if( saturated ) {
      if( method %in% c("NR","FS") ) method <- "EM"
      if( method == "approxBayes" ) stop( gettext(
         'method="approxBayes" not available when saturated=TRUE'),
	 domain = NA )
   }
   #--------------------------------------------
   # control parameters
   control <- do.call( "rrLogitControl", control )
   if( saturated && method == "MCMC" && control$typeMCMC == "RWM" )
      control$typeMCMC <- "DA"
   ctrlInt <- c( iterMaxNR = control$iterMaxNR,
      iterMaxFS = control$iterMaxFS,
      iterMaxEM = control$iterMaxEM,
      iterMaxMstep = control$iterMaxMstep,
      iterApproxBayes = control$iterApproxBayes,
      imputeApproxBayesInt = as.integer( control$imputeApproxBayes ),
      iterMCMC = control$iterMCMC,
      burnMCMC = control$burnMCMC,
      thinMCMC = control$thinMCMC,
      imputeEvery = control$imputeEvery,
      typeMCMCInt = match( control$typeMCMC, c("DA","RWM"), nomatch=0L ),
      stuckLimit = control$stuckLimit )
   storage.mode(ctrlInt) <- "integer"
   ctrlReal <- c( critConverge = control$critConverge,
      critBoundary = control$critBoundary,   # 2
      dfDA = control$tuneDA[1L],             # 3
      stepSizeDA = control$tuneDA[2L],       # 4
      scaleFacDA = control$tuneDA[3L],       # 5
      dfRWM = control$tuneRWM[1L],           # 6
      scaleFacRWM = control$tuneRWM[2L],     # 7
      startValJitter = control$startValJitter ) # 8
   storage.mode( ctrlReal ) <- "double"
   #--------------------------------------------
   mvcode <- .Machine$double.xmax
   specialCodes <- c( mvcode=mvcode, nancode=NaN, infcode=Inf, neginfcode=-Inf)
   #--------------------------------------------
   # arrays for holding results from Fortran call
   outputsInt <- integer(6L)
   names(outputsInt) <- c("iter", "convergedInt", "boundaryInt",
      "abortedInt", "vhatFailedInt", "DAPFailedInt" )
   outputsReal <- numeric(2L)
   names(outputsReal) <- c("loglik", "logprior")
   storage.mode(outputsReal) <- "double"
   coefficients <- matrix( numeric(1L), pD$dimVec["ncolMM"],
      pD$dimVec["nLevels"] )
   rownames( coefficients ) <- colnames( pD$modelMatrix )
   colnames( coefficients ) <- pD$responseLevels
   coefVec <- numeric( pD$dimVec["nParamThisModel"] )
   names(coefVec) <- paste(
      rep( colnames( pD$modelMatrix ), pD$dimVec["nLevels"] - 1L ),
      rep( pD$responseLevels[ -baselineInt ], each=pD$dimVec["ncolMM"] ),
      sep="." )
   score <- coefVec
   vhatCoef <- matrix( numeric(1L), pD$dimVec["nParamThisModel"],
      pD$dimVec["nParamThisModel"] )
   rownames(vhatCoef) <- colnames(vhatCoef) <- names(coefVec)
   hessA <- hess <- vhatCoef
   fittedPi <- matrix( 1 / pD$dimVec["nLevels"], pD$dimVec["nCovPatt"],
      pD$dimVec["nLevels"] )  # uniform probabilities
   rownames(fittedPi) <- NULL
   colnames(fittedPi) <- pD$responseLevels
   storage.mode(coefficients) <- storage.mode(coefVec) <-
      storage.mode(vhatCoef) <- storage.mode(fittedPi) <- "double"
   fitweightMat <- fstarMat <- fhatMat <- fittedPistar <- fittedPi
   margProportionsYstar <- numeric( pD$dimVec["nLevels"] )
   names(margProportionsYstar) <- pD$responseLevels
   storage.mode(margProportionsYstar) <- "double"
   margProportionsY <- margProportionsYstar
   #--------------------------------------------
   if( saturated ) {
      if( !is.null( startVal ) ) {
         if( ( ! is.numeric(startVal) ) || ( ! is.matrix(startVal) ) ||
            ( NROW(startVal) != pD$dimVec["nCovPatt"] ) || 
            ( NCOL(startVal) != pD$dimVec["nLevels"] ) ||
            any( startVal[ ! pD$emptyCovPatt, ] < 0 ) ||
            any( startVal[ ! pD$emptyCovPatt, ] > 1 ) || 
            any( abs( apply( startVal[ ! pD$emptyCovPatt, ], 1, sum ) 
	        - 1 ) > 1e-08 ) ) stop( gettext(
                "'startVal' is not a valid matrix of probabilities" ),
                domain = NA )
         fittedPi[,] <- startVal[,]
         startValSource <- "supplied by user"
      } else {
         startValSource <- "default"
      }
   } else {
      if( !is.null( startVal ) ) {
         if( ( ! is.numeric(startVal) ) || ( ! is.matrix(startVal) ) ||
            ( NROW(startVal) != pD$dimVec["ncolMM"] ) || 
            ( NCOL(startVal) != pD$dimVec["nLevels"] ) ||
            any( startVal[,baselineInt] != 0 ) ) stop( gettext(
                "'startVal' is not a valid matrix of coefficients" ),
                domain = NA )
         coefficients[,] <- startVal[,]
         startValSource <- "supplied by user"
      } else {
         startValSource <- "default"
      }
   }
   startValSourceInt <- match( startValSource,
      c("default", "supplied by user", "from rrLogit object"),
      nomatch=0L)
   if( startValSourceInt == 0L ) stop( gettextf(
      "startValSource '%s' not recognized", startValSource ), domain = NA )
   #--------------------------------------------
   if( method == "MCMC" && control$typeMCMC == "RWM" ) stop( gettext(
     "Cannot run random-walk Metropolis without first getting a mode"),
      domain = NA )
   if( method == "approxBayes" ) stop( gettext(
     "Cannot run approximate Bayes without first getting a mode"),
      domain = NA )
   #--------------------------------------------
   # needed for MCMC or approxBayes
   iterMCMC <- as.integer( ceiling( control$iterMCMC / control$thinMCMC ) *
      control$thinMCMC )
   if( method == "MCMC" ) {
      seriesLength <- iterMCMC / control$thinMCMC
      nImpute <- if( control$imputeEvery == 0L ) 0L else
         floor( iterMCMC / control$imputeEvery )
   } else if( method == "approxBayes" ) { 
      seriesLength <- control$iterApproxBayes
      nImpute <- if( control$imputeApproxBayes )
         control$iterApproxBayes else 0L
   } else {
      seriesLength <- 0L
      nImpute <- 0L
   }
   dimVecMCMC <- c(
      seriesLength = seriesLength,
      nImpute = nImpute )
   storage.mode( dimVecMCMC ) <- "integer"
   if( saturated ) {
      vhatCoefRWM <- matrix( numeric(), 0L, 0L )
      coefMCMC <- matrix( numeric(), 0L, 0L )
      coefVecMCMC <- numeric()
      vhatCoefMCMC <- matrix( numeric(), 0L, 0L )
      coefVecSeries <- matrix( numeric(), seriesLength, 0L )
   } else {
      vhatCoefRWM <- vhatCoef
      coefMCMC <- coefficients
      coefVecMCMC <- coefVec
      vhatCoefMCMC <- vhatCoef
      coefVecSeries <- matrix( numeric(), seriesLength,
         pD$dimVec["nParamThisModel"] )
      colnames(coefVecSeries) <- names(coefVec)
   }
   approxBayesLogImpRatios <- numeric( seriesLength )
   fittedPistarMCMC <- fittedPiMCMC <- fittedPi
   nActual <- c( nIterActual = 0L, nSampleActual = 0L, nImpActual = 0L )
   startLogP <- mhAcceptRate <- numeric(1L)
   logPSeries <- numeric(seriesLength)
   impMatSeries <- array( integer(), c( pD$dimVec["nrowInputData"],
      pD$dimVec["nLevels"], nImpute ) )
   dimnames(impMatSeries) <- list( rownames(mF), pD$responseLevels, NULL )
   impVecSeries <- matrix( integer(), pD$dimVec["nrowInputData"], nImpute )
   rownames(impVecSeries) <- rownames(mF)
   piMargMCMC <- numeric( pD$dimVec["nLevels"] )
   names(piMargMCMC) <- pD$responseLevels
   piMargSeries <- matrix( numeric(1L), seriesLength, pD$dimVec["nLevels"] )
   colnames(piMargSeries) <- pD$responseLevels
   #--------------------------------------------
   # integerized frequencies needed for random imputation within MCMC
   freqRowInputDataInt <- as.integer(pD$freqRowInputData)
   if( any( freqRowInputDataInt != pD$freqRowInputData ) &
      ( ! ( method %in% c("EM","NR","FS") ) ) ) stop( gettextf(
      "Method '%s' not defined for non-integer frequencies", method),
       domain = NA ) 
   fMatInt <- as.integer(pD$fMat)
   if( any( fMatInt != pD$fMatInt ) &
      ( ! ( method %in% c("EM","NR","FS") ) ) ) stop( gettextf(
      "Method '%s' not defined for non-integer frequencies", method),
       domain = NA ) 
   responseRowInputData <- if( wideFormat ) integer() else
      pD$responseRowInputData
   freqForDataPattInt <- as.integer( pD$freqForDataPatt )
   if( any( freqForDataPattInt != pD$freqForDataPatt ) &
      ( ! ( method %in% c("EM","NR","FS") ) ) ) stop( gettextf(
      "Method '%s' not defined for non-integer frequencies", method),
       domain = NA ) 
   #--------------------------------------------
   # needed for dotCall64::.C64
   SIGNATURE <- c(
      # inputs
      dimVec = "integer",
      modelMatrix = "double",
      fitweightRowInputData = "double",
      fitweightCovPatt = "double",
      fitweightDataPatt = "double",
      freqForDataPatt = "double",
      covPatt = "integer",
      dataPatt = "integer",
      covPattForDataPatt = "integer",
      responseForDataPatt = "integer",
      surveyModeInt = "integer",
      baselineInt = "integer",
      pertMat = "double",
      pertMatInv = "double",
      priorInt = "integer",
      priorFreqTot = "double",
      priorAllocSuppliedInt = "integer",
      priorAlloc = "double",
      saturatedInt = "integer",
      methodInt = "integer",
      ctrlInt = "integer",
      ctrlReal = "double",
      specialCodes = "double",
      # MCMC inputs
      dimVecMCMC = "integer",
      vhatCoefRWM = "double",
      microDataInt = "integer",
      freqRowInputDataInt = "integer",
      fMatInt = "integer",
      responseRowInputData = "integer",
      freqForDataPattInt = "integer",
      # starting values indicator
      startValSourceInt = "integer",
      # outputs
      outputsInt = "integer",
      outputsReal = "double",
      score = "double",
      hess = "double",
      coefficients = "double",
      coefVec = "double",
      vhatCoef = "double",
      fittedPi = "double",
      fittedPistar = "double",
      margProportionsYstar = "double",
      margProportionsY = "double",
      fhatMat = "double",
      hessA = "double",
      fstarMat = "double",
      fitweightMat = "double",
      # MCMC outputs
      coefMCMC = "double",
      coefVecMCMC = "double",
      vhatCoefMCMC = "double",
      fittedPiMCMC = "double",
      fittedPistarMCMC = "double",
      nActual = "integer",
      mhAcceptRate = "double",
      startLogP = "double",
      logPSeries = "double",
      coefVecSeries = "double",
      impMatSeries = "integer",
      impVecSeries = "integer",
      piMargMCMC = "double",
      piMargSeries = "double",
      approxBayesLogImpRatios = "double",
      # messaging
      status = "integer",
      msgLenMax = "integer",
      msgCodes = "integer",
      msgLenActual = "integer" )
   INTENT <- c(
      # inputs
      dimVec = "r",
      modelMatrix = "r",
      fitweightRowInputData = "r",
      fitweightCovPatt = "r",
      fitweightDataPatt = "r",
      freqForDataPatt = "r",
      covPatt = "r",
      dataPatt = "r",
      covPattForDataPatt = "r",
      responseForDataPatt = "r",
      surveyModeInt = "r",
      baselineInt = "r",
      pertMat = "r",
      pertMatInv = "r",
      priorInt = "r",
      priorFreqTot = "r",
      priorAllocSuppliedInt = "r",
      priorAlloc = "r",
      saturatedInt = "r",
      methodInt = "r",
      ctrlInt = "r",
      ctrlReal = "r",
      specialCodes = "r",
      # MCMC inputs
      dimVecMCMC = "r",
      vhatCoefRWM = "r",
      microDataInt = "r",
      freqRowInputDataInt = "r",
      fMatInt = "r",
      responseRowInputData = "r",
      freqForDataPattInt = "r",
      # starting values indicator
      startValSourceInt = "r",
      # outputs
      outputsInt = "w",
      outputsReal = "w",
      score = "w",
      hess = "w",
      coefficients = "rw",
      coefVec = "w",
      vhatCoef = "w",
      fittedPi = "rw",
      fittedPistar = "w",
      margProportionsYstar = "w",
      margProportionsY = "w",
      fhatMat = "w",
      hessA = "w",
      fstarMat = "w",
      fitweightMat = "w",
      # MCMC outputs
      coefMCMC = "w",
      coefVecMCMC = "w",
      vhatCoefMCMC = "w",
      fittedPiMCMC = "w",
      fittedPistarMCMC = "w",
      nActual = "w",
      mhAcceptRate = "w",
      startLogP = "w",
      logPSeries = "w",
      coefVecSeries = "w",
      impMatSeries = "w",
      impVecSeries = "w",
      piMargMCMC = "w",
      piMargSeries = "w",
      approxBayesLogImpRatios = "w",
      # messaging
      status = "w",
      msgLenMax = "r",
      msgCodes = "w",
      msgLenActual = "w" )
   #--------------------------------------------
   # create a matrix for holding message codes
   msgLenMax <- 40L
   msgCodes <- matrix( 0L, msgLenMax, 17L )
   #--------------------------------------------
   result <- dotCall64::.C64("rrlogit",
      SIGNATURE = SIGNATURE,
      # inputs
      dimVec = pD$dimVec,
      modelMatrix = pD$modelMatrix,
      fitweightRowInputData = fitweightRowInputData,
      fitweightCovPatt = fitweightCovPatt,
      fitweightDataPatt = fitweightDataPatt,
      freqForDataPatt = pD$freqForDataPatt,
      covPatt = pD$covPatt,
      dataPatt = if( pD$wideFormat ) integer() else pD$dataPatt,
      covPattForDataPatt = pD$covPattForDataPatt,
      responseForDataPatt = pD$responseForDataPatt,
      surveyModeInt = as.integer(surveyMode),
      baselineInt = as.integer(baselineInt),
      pertMat = pertMat,
      pertMatInv = pertMatInv,
      priorInt = priorInt,
      priorFreqTot = priorFreqTot,
      priorAllocSuppliedInt = as.integer(priorAllocSupplied),
      priorAlloc = priorAlloc,
      saturatedInt = saturatedInt,
      methodInt = methodInt,
      ctrlInt = ctrlInt,
      ctrlReal = ctrlReal,
      specialCodes = specialCodes,
      # MCMC inputs
      dimVecMCMC = dimVecMCMC,
      vhatCoefRWM = vhatCoefRWM,
      microDataInt = as.integer(pD$microData),
      freqRowInputDataInt = freqRowInputDataInt,
      fMatInt = fMatInt,
      responseRowInputData = responseRowInputData,
      freqForDataPattInt = freqForDataPattInt,
      # starting values indicator
      startValSourceInt = startValSourceInt,
      # outputs
      outputsInt = outputsInt,
      outputsReal = outputsReal,
      score = score,
      hess = hess,
      coefficients = coefficients,
      coefVec = coefVec,
      vhatCoef = vhatCoef,
      fittedPi = fittedPi,
      fittedPistar = fittedPistar,
      margProportionsYstar = margProportionsYstar,
      margProportionsY = margProportionsY,
      fhatMat = fhatMat,
      hessA = hessA,
      fstarMat = fstarMat,
      fitweightMat = fitweightMat,
      # MCMC outputs
      coefMCMC = coefMCMC,
      coefVecMCMC = coefVecMCMC,
      vhatCoefMCMC = vhatCoefMCMC,
      fittedPiMCMC = fittedPiMCMC,
      fittedPistarMCMC = fittedPistarMCMC,
      nActual = nActual,
      mhAcceptRate = mhAcceptRate,
      startLogP = startLogP,
      logPSeries = logPSeries,
      coefVecSeries = coefVecSeries,
      impMatSeries = impMatSeries,
      impVecSeries = impVecSeries,
      piMargMCMC = piMargMCMC,
      piMargSeries = piMargSeries,
      approxBayesLogImpRatios = approxBayesLogImpRatios,
      # messaging
      status = integer(1L),
      msgLenMax = msgLenMax,
      msgCodes = msgCodes,
      msgLenActual = integer(1L),
      # other args to .C64
      NAOK = TRUE,
#      INTENT = INTENT,
      PACKAGE = "rrLogit" )
   #--------------------------------------------
   # display message from Fortran, if present
   msg.lines <- .msg( result$msgCodes, result$msgLenActual )
   if( is.null( msg.lines ) ){
      msg <- "OK"
   } else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( result$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   if( saturated ) {
      result$fittedPi[ result$fittedPi == specialCodes["mvcode"] ] <- NA
      result$fittedPistar[ result$fittedPistar == specialCodes["mvcode"] ] <- NA
      result$fhatMat[ result$fhatMat == specialCodes["mvcode"] ] <- NA
   }
   #--------------------------------------------
   emptyCovPatt <- ( fitweightCovPatt == 0 )
   tmp <- result$fittedPi * fitweightCovPatt
   tmp <- apply( tmp[ ! emptyCovPatt,,drop=FALSE ], 2, sum )
   result$fittedPiMarg <- tmp / sum(tmp)
   #--------------------------------------------
   result$call <- match.call()
   result$mF <- mF
   result$formula <- obj
   result$formulaStr <- deparse(obj)
   result$responseVarName <- .responseVarName(obj)
   result$contrasts <- attr(pD$modelMatrix, "contrasts")
   result$xlevels <- .getXlevels( attr(mF, "terms"), mF )
   result$prior <- prior
   result$wideFormat <- wideFormat
   result$surveyMode <- surveyMode
   result$freqSupplied <- pD$freqSupplied
   result$weightSupplied <- pD$weightSupplied
   result$dataFreqTot <- sum( pD$freqForDataPatt )
   if( surveyMode ) {
      result$weightTot <- sum( pD$weightRowInputData )
      result$weightScaleFac <- weightScaleFac
   }
   result$responseLevels <- pD$responseLevels
   result$freqRowInputData <- pD$freqRowInputData
   result$weightRowInputData <- pD$weightRowInputData
   result$freqForCovPatt <- pD$freqForCovPatt
   result$weightForCovPatt <- pD$weightForCovPatt
   result$freqForDataPatt <- pD$freqForDataPatt
   result$weightForDataPatt <- pD$weightForDataPatt
   result$dataPatt <- pD$dataPatt
   result$startValSource <- startValSource
   result$priorFreqTotSupplied <- priorFreqTotSupplied
   result$priorAllocSupplied <- priorAllocSupplied
   result$saturated <- saturated
   result$emptyCovPatt <- pD$emptyCovPatt
   result$rowNamesInputData <- pD$rowNamesInputData
   result$fMat <- pD$fMat
   result$microData <- pD$microData
   #--------------------------------------------
   result$method <- method
   result$control <- control
   result$iter <- result$outputsInt[ "iter" ]
   result$converged <- as.logical( result$outputsInt[ "convergedInt" ] )
   result$boundary <- as.logical( result$outputsInt[ "boundaryInt" ] )
   result$aborted <- as.logical( result$outputsInt[ "abortedInt" ] )
   result$vhatFailed <- as.logical( result$outputsInt[ "vhatFailedInt" ] )
   result$DAPFailed <- as.logical( result$outputsInt[ "DAPFailedInt" ] )
   result$loglik <- result$outputsReal[ "loglik" ]
   result$logprior <- result$outputsReal[ "logprior" ]
   result$logP <- result$loglik + result$logprior
   names(result$loglik) <- names(result$logprior) <-
      names(result$logP) <- NULL
   #--------------------------------------------
   result$beganAtMode <- FALSE
   result$atMode <- ( method %in% c("EM","NR","FS" ) ) & result$converged
   result$vhatCoefRWMStore <- if( ( !saturated ) & result$atMode &
      ( !result$vhatFailed ) ) result$vhatCoef else NULL
   #--------------------------------------------
   result$nIterActual <- result$nActual[1L]
   result$nSampleActual <- result$nActual[2L]
   result$nImpActual <- result$nActual[3L]
   #--------------------------------------------
   result$logPSeries[ result$logPSeries == mvcode ] <- NA
   result$coefVecSeries[ result$coefVecSeries == mvcode ] <- NA
   result$impMatSeries[ result$impMatSeries == -1L ] <- NA
   result$impVecSeries[ result$impVecSeries == -1L ] <- NA
   result$approxBayesLogImpRatios[
      result$approxBayesLogImpRatios == mvcode ] <- NA
   #--------------------------------------------
   result$SIGNATURE <- SIGNATURE
   result$INTENT <- INTENT
   #--------------------------------------------
   if( result$DAPFailed ) {
      warning( gettext(
         "Model fitting not attempted" ), domain = NA )
   } else {
      if( method %in% c("EM","NR","FS") ) {
         if( ! result$converged ) warning( gettextf(
            "Procedure failed to converge by iteration %i", result$iter ),
            domain = NA )
         if( result$boundary ) {
            if( ! result$saturated ) {
               message( gettext(
            "Estimate at or near boundary; standard errors may be unreliable" ),
               domain = NA )
            }
         }
      }
   }
   #--------------------------------------------
   structure( result, class = c("rrLogit", "list" ) )
}

rrLogit.rrLogit <- function( obj, method = obj$method,
   saturated = obj$saturated, control = NULL, startVal = NULL, ...) {
   #----------------------------------
   stopifnot( inherits(obj, "rrLogit") )
   #----------------------------------
   argsSupplied <- names( as.list( match.call() )[-1L] )
   badArgs <- setdiff( argsSupplied,
      c("obj", "method", "saturated", "control", "startVal") )
   for( i in seq_along(badArgs) ) warning( gettextf(
      "Argument '%s' was ignored", badArgs[i] ), domain=NA ) 
   #----------------------------------
   stopifnot( method %in% c("EM", "NR", "FS", "MCMC", "approxBayes") )
   methodInt <- match( method, c("EM", "NR", "FS", "MCMC", "approxBayes") )
   if( ( ! (method %in% c("EM","NR","FS")) ) & obj$surveyMode ) stop( gettextf(
      'Method "%s" not defined with survey weights', method), domain = NA )
   saturated <- as.logical(saturated)[1L]
   saturatedInt <- as.integer(saturated)
   #----------------------------------
   # ensure modelMatrix has full rank
   if( ! saturated ) {
      if( qr(obj$modelMatrix)$rank < NCOL(obj$modelMatrix) ) stop( gettext(
         "modelMatrix does not have full rank" ), domain=NA )
      emptyCovPatt <- ( obj$fitweightCovPatt == 0 )
      if( qr(obj$modelMatrix[(!emptyCovPatt),,drop=FALSE])$rank <
         NCOL(obj$modelMatrix) ) stop( gettext(
       "modelMatrix does not have full rank due to empty covariate patterns" ), 
         domain=NA )
   }      
   #--------------------------------------------
   if( ( obj$prior == "DAP" ) & ( obj$DAPFailed ) ) stop( gettext(
      "Default prior allocation for DAP has failed" ), domain = NA )
   #----------------------------------
   if( is.null(control) ) {
      control <- obj$control
   } else {
      stopifnot( is.list(control) )
      ctrl <- obj$control
      for( i in seq_along(control) ) ctrl[[ names(control)[i] ]] <- 
         control[[i]]
      control <- do.call( "rrLogitControl", ctrl )
   }
   if( saturated && method == "MCMC" && control$typeMCMC == "RWM" )
      control$typeMCMC <- "DA"
   ctrlInt <- c( iterMaxNR = control$iterMaxNR,
      iterMaxFS = control$iterMaxFS,
      iterMaxEM = control$iterMaxEM,
      iterMaxMstep = control$iterMaxMstep,
      iterApproxBayes = control$iterApproxBayes,
      imputeApproxBayesInt = as.integer( control$imputeApproxBayes ),
      iterMCMC = control$iterMCMC,
      burnMCMC = control$burnMCMC,
      thinMCMC = control$thinMCMC,
      imputeEvery = control$imputeEvery,
      typeMCMCInt = match( control$typeMCMC, c("DA","RWM"), nomatch=0L ),
      stuckLimit = control$stuckLimit )
   storage.mode(ctrlInt) <- "integer"
   ctrlReal <- c( critConverge = control$critConverge,
      critBoundary = control$critBoundary,   # 2
      dfDA = control$tuneDA[1L],             # 3
      stepSizeDA = control$tuneDA[2L],       # 4
      scaleFacDA = control$tuneDA[3L],       # 5
      dfRWM = control$tuneRWM[1L],           # 6
      scaleFacRWM = control$tuneRWM[2L],     # 7
      startValJitter = control$startValJitter ) # 8
   storage.mode( ctrlReal ) <- "double"
   #----------------------------------
   coefficients <- obj$coefficients
   fittedPi <- obj$fittedPi
   startValSource <- "from rrLogit object"
   if( saturated ) {
      if( ! is.null(startVal) ) {
         if( ( ! is.numeric(startVal) ) || ( ! is.matrix(startVal) ) ||
            ( NROW(startVal) != obj$dimVec["nCovPatt"] ) || 
            ( NCOL(startVal) != obj$dimVec["nLevels"] ) ||
            any( startVal[ ! obj$emptyCovPatt, ] < 0 ) ||
            any( startVal[ ! obj$emptyCovPatt, ] > 1 ) || 
            any( abs( apply( startVal[ ! obj$emptyCovPatt, ], 1, sum ) 
	        - 1 ) > 1e-08 ) ) stop( gettext(
                "'startVal' is not a valid matrix of probabilities" ),
                domain = NA )
         fittedPi[,] <- startVal[,]
	 startValSource <- "supplied by user"
     }
   } else {
      if( !is.null( startVal ) ) {
         if( ( ! is.numeric(startVal) ) || ( ! is.matrix(startVal) ) ||
            ( NROW(startVal) != obj$dimVec["ncolMM"] ) || 
            ( NCOL(startVal) != obj$dimVec["nLevels"] ) ||
            any( startVal[,obj$baselineInt] != 0 ) ) stop( gettext(
                "'startVal' is not a valid matrix of coefficients" ),
                domain = NA )
         coefficients[,] <- startVal[,]
	 startValSource <- "supplied by user"
      } else {
         if( obj$saturated ) {
	    coefficients[,] <- 0
	    warning( gettext(
	       "Default starting values used because obj$saturated is TRUE" ),
	       domain = NA )
	    startValSource <- "default"
         }
      }
   }
   startValSourceInt <- match( startValSource,
      c("default", "supplied by user", "from rrLogit object"),
      nomatch=0L)
   if( startValSourceInt == 0L ) stop( gettextf(
      "startValSource '%s' not recognized", startValSource ), domain = NA )
   vhatCoef <- obj$vhatCoef
   coefVec <- obj$coefVec
   mvcode <- obj$specialCodes["mvcode"]
   fittedPi[ is.na(fittedPi) ] <- mvcode
   #----------------------------------
   doingRWM <-  ( method == "MCMC" ) & ( control$typeMCMC == "RWM" )
   doingApproxBayes <- ( method == "approxBayes" )
   if( saturated & doingRWM ) stop( gettext(
      "Random-walk Metropolis not defined for saturated model" ),
      domain = NA )
   if( saturated & doingApproxBayes ) stop( gettext(
      "Approx Bayes not defined for saturated model" ),
      domain = NA )
   atMode <- obj$atMode
   if( ( ! saturated ) & obj$saturated ) atMode <- FALSE
   if( saturated  & ( ! obj$saturated ) ) atMode <- FALSE
   if( doingRWM & ( is.null(obj$vhatCoefRWMStore) | ( ! atMode ) ) ) 
      stop( gettext(
      "Cannot run random-walk Metropolis without starting at a mode"),
          domain = NA )
   if( doingApproxBayes & ( is.null(obj$vhatCoefRWMStore) | ( ! atMode ) ) ) 
      stop( gettext(
      "Cannot run approx Bayes without starting at a mode"),
          domain = NA )
   #--------------------------------------------
   # needed for MCMC or approxBayes
   iterMCMC <- as.integer( ceiling( control$iterMCMC / control$thinMCMC ) *
      control$thinMCMC )
   if( method == "MCMC" ) {
      seriesLength <- iterMCMC / control$thinMCMC
      nImpute <- if( control$imputeEvery == 0L ) 0L else
         floor( iterMCMC / control$imputeEvery )
   } else if( method == "approxBayes" ) { 
      seriesLength <- control$iterApproxBayes
      nImpute <- if( control$imputeApproxBayes )
         control$iterApproxBayes else 0L
   } else {
      seriesLength <- 0L
      nImpute <- 0L
   }
   dimVecMCMC <- c(
      seriesLength = seriesLength,
      nImpute = nImpute )
   storage.mode( dimVecMCMC ) <- "integer"
   #
   if( saturated ) {
      vhatCoefRWM <- matrix( numeric(), 0L, 0L )
      coefMCMC <- matrix( numeric(), 0L, 0L )
      coefVecMCMC <- numeric()
      vhatCoefMCMC <- matrix( numeric(), 0L, 0L )
      coefVecSeries <- matrix( numeric(), seriesLength, 0L )
   } else {
      vhatCoefRWM <-  if( doingRWM | doingApproxBayes )
         obj$vhatCoefRWMStore else vhatCoef
      coefMCMC <- coefficients
      coefVecMCMC <- coefVec
      vhatCoefMCMC <- vhatCoef
      coefVecSeries <- matrix( numeric(), seriesLength,
         obj$dimVec["nParamThisModel"] )
      colnames(coefVecSeries) <- names(coefVec)
   }
   approxBayesLogImpRatios <- numeric( seriesLength )
   fittedPistarMCMC <- fittedPiMCMC <- fittedPi
   nActual <- c( nIterActual = 0L, nSampleActual = 0L, nImpActual = 0L )
   startLogP <- mhAcceptRate <- numeric(1L)
   logPSeries <- numeric(seriesLength)
   impMatSeries <- array( integer(), c( obj$dimVec["nrowInputData"],
      obj$dimVec["nLevels"], nImpute ) )
   dimnames(impMatSeries) <- list( rownames(obj$mF), obj$responseLevels, NULL )
   impVecSeries <- matrix( integer(), obj$dimVec["nrowInputData"], nImpute )
   rownames(impVecSeries) <- rownames(obj$mF)
   piMargMCMC <- numeric( obj$dimVec["nLevels"] )
   names(piMargMCMC) <- obj$responseLevels
   piMargSeries <- matrix( numeric(1L), seriesLength, obj$dimVec["nLevels"] )
   colnames(piMargSeries) <- obj$responseLevels
   #----------------------------------
   # create a matrix for holding message codes
   msgLenMax <- 40L
   msgCodes <- matrix( 0L, msgLenMax, 17L )
   #--------------------------------------------
   result <- dotCall64::.C64("rrlogit",
      SIGNATURE = obj$SIGNATURE,
      # inputs
      dimVec = obj$dimVec,
      modelMatrix = obj$modelMatrix,
      fitweightRowInputData = obj$fitweightRowInputData,
      fitweightCovPatt = obj$fitweightCovPatt,
      fitweightDataPatt = obj$fitweightDataPatt,
      freqForDataPatt = obj$freqForDataPatt,
      covPatt = obj$covPatt,
      dataPatt = if( obj$wideFormat ) integer() else obj$dataPatt,
      covPattForDataPatt = obj$covPattForDataPatt,
      responseForDataPatt = obj$responseForDataPatt,
      surveyModeInt = obj$surveyModeInt,
      baselineInt = obj$baselineInt,
      pertMat = obj$pertMat,
      pertMatInv = obj$pertMatInv,
      priorInt = obj$priorInt,
      priorFreqTot = obj$priorFreqTot,
      priorAllocSuppliedInt = if( obj$prior == "DAP" ) 1L else 0L,
      priorAlloc = obj$priorAlloc,
      saturatedInt = saturatedInt,
      methodInt = methodInt,
      ctrlInt = ctrlInt,
      ctrlReal = ctrlReal,
      specialCodes = obj$specialCodes,
      # MCMC inputs
      dimVecMCMC = dimVecMCMC,
      vhatCoefRWM = vhatCoefRWM,
      microDataInt = obj$microDataInt,
      freqRowInputDataInt = obj$freqRowInputDataInt,
      fMatInt = obj$fMatInt,
      responseRowInputData = obj$responseRowInputData,
      freqForDataPattInt = obj$freqForDataPattInt,
      # starting values indicator
      startValSourceInt = startValSourceInt,
      # outputs
      outputsInt = obj$outputsInt,
      outputsReal = obj$outputsReal,
      score = obj$score,
      hess = obj$hess,
      coefficients = coefficients,
      coefVec = coefVec,
      vhatCoef = vhatCoef,
      fittedPi = fittedPi,
      fittedPistar = obj$fittedPistar,
      margProportionsYstar = obj$margProportionsYstar,
      margProportionsY = obj$margProportionsY,
      fhatMat = obj$fhatMat,
      hessA = obj$hessA,
      fstarMat = obj$fstarMat,
      fitweightMat = obj$fitweightMat,
      # MCMC outputs
      coefMCMC = coefMCMC,
      coefVecMCMC = coefVecMCMC,
      vhatCoefMCMC = vhatCoefMCMC,
      fittedPiMCMC = fittedPiMCMC,
      fittedPistarMCMC = fittedPistarMCMC,
      nActual = nActual,
      mhAcceptRate = mhAcceptRate,
      startLogP = startLogP,
      logPSeries = logPSeries,
      coefVecSeries = coefVecSeries,
      impMatSeries = impMatSeries,
      impVecSeries = impVecSeries,
      piMargMCMC = piMargMCMC,
      piMargSeries = piMargSeries,
      approxBayesLogImpRatios = approxBayesLogImpRatios,
      # messaging
      status = integer(1L),
      msgLenMax = msgLenMax,
      msgCodes = msgCodes,
      msgLenActual = integer(1L),
      # other args to .C64
      NAOK = TRUE,
#      INTENT = obj$INTENT,
      PACKAGE = "rrLogit")
   #--------------------------------------------
   # display message from Fortran, if present
   msg.lines <- .msg( result$msgCodes, result$msgLenActual )
   if( is.null( msg.lines ) ){
      msg <- "OK"
   } else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( result$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   if( saturated ) {
      result$fittedPi[ result$fittedPi == mvcode ] <- NA
      result$fittedPistar[ result$fittedPistar == mvcode ] <- NA
      result$fhatMat[ result$fhatMat == mvcode ] <- NA
   }
   #--------------------------------------------
   emptyCovPatt <- ( obj$fitweightCovPatt == 0 )
   tmp <- result$fittedPi * obj$fitweightCovPatt
   tmp <- apply( tmp[ ! emptyCovPatt,,drop=FALSE ], 2, sum )
   result$fittedPiMarg <- tmp / sum(tmp)
   #--------------------------------------------
   result$call <- match.call()
   result$mF <- obj$mF
   result$formula <- obj$formula
   result$formulaStr <- obj$formulaStr
   result$responseVarName <- obj$responseVarName
   result$contrasts <- obj$contrasts
   result$xlevels <- obj$xlevels
   result$prior <- obj$prior
   result$wideFormat <- obj$wideFormat
   result$surveyMode <- obj$surveyMode
   result$freqSupplied <- obj$freqSupplied
   result$weightSupplied <- obj$weightSupplied
   result$dataFreqTot <- obj$dataFreqTot
   if( result$surveyMode ) {
      result$weightTot <- obj$weightTot
      result$weightScaleFac <- obj$weightScaleFac
   }
   result$responseLevels <- obj$responseLevels
   result$freqRowInputData <- obj$freqRowInputData
   result$weightRowInputData <- obj$weightRowInputData
   result$freqForCovPatt <- obj$freqForCovPatt
   result$weightForCovPatt <- obj$weightForCovPatt
   result$freqForDataPatt <- obj$freqForDataPatt
   result$weightForDataPatt <- obj$weightForDataPatt
   result$dataPatt <- obj$dataPatt
   result$startValSource <- startValSource
   result$priorFreqTotSupplied <- obj$priorFreqTotSupplied
   result$priorAllocSupplied <- obj$priorAllocSupplied
   result$saturated <- saturated
   result$emptyCovPatt <- obj$emptyCovPatt
   result$rowNamesInputData <- obj$rowNamesInputData
   result$fMat <- obj$fMat
   result$responseRowInputData <- obj$responseRowInputData
   result$microData <- obj$microData
   #--------------------------------------------
   result$method <- method
   result$control <- control
   result$iter <- result$outputsInt[ "iter" ]
   result$converged <- as.logical( result$outputsInt[ "convergedInt" ] )
   result$boundary <- as.logical( result$outputsInt[ "boundaryInt" ] )
   result$aborted <- as.logical( result$outputsInt[ "abortedInt" ] )
   result$vhatFailed <- as.logical( result$outputsInt[ "vhatFailedInt" ] )
   result$DAPFailed <- as.logical( result$outputsInt[ "DAPFailedInt" ] )
   result$loglik <- result$outputsReal[ "loglik" ]
   result$logprior <- result$outputsReal[ "logprior" ]
   result$logP <- result$loglik + result$logprior
   names(result$loglik) <- names(result$logprior) <-
      names(result$logP) <- NULL
   #--------------------------------------------
   result$beganAtMode <- atMode
   result$atMode <- ( method %in% c("EM","NR","FS" ) ) & result$converged
   result$vhatCoefRWMStore <- if( ( !saturated ) & result$atMode &
      ( !result$vhatFailed ) ) result$vhatCoef else obj$vhatCoefRWMStore
   #--------------------------------------------
   result$nIterActual <- result$nActual[1L]
   result$nSampleActual <- result$nActual[2L]
   result$nImpActual <- result$nActual[3L]
   #--------------------------------------------
   result$logPSeries[ result$logPSeries == mvcode ] <- NA
   result$coefVecSeries[ result$coefVecSeries == mvcode ] <- NA
   result$impMatSeries[ result$impMatSeries == -1L ] <- NA
   result$impVecSeries[ result$impVecSeries == -1L ] <- NA
   result$approxBayesLogImpRatios[
      result$approxBayesLogImpRatios == mvcode ] <- NA
   #--------------------------------------------
   result$SIGNATURE <- obj$SIGNATURE
   result$INTENT <- obj$INTENT
   #--------------------------------------------
   if( result$DAPFailed ) {
      warning( gettext(
         "Model fitting not attempted" ), domain = NA )
   } else {
      if( method %in% c("EM","NR","FS") ) {
         if( ! result$converged ) warning( gettextf(
            "Procedure failed to converge by iteration %i", result$iter ),
            domain = NA )
         if( result$boundary ) {
            if( ! result$saturated ) {
               message( gettext(
            "Estimate at or near boundary; standard errors may be unreliable" ),
               domain = NA )
            }
         }
      }
   }
   #--------------------------------------------
   structure( result, class = c("rrLogit", "list" ) )
}


summary.rrLogit <- function(object, showCoef=FALSE, digits=4L,
      dispersion=1.0, ...) {
   stopifnot( inherits(object, "rrLogit") )
   if( inherits(object, "summary.rrLogit") ) return(object)
   result <- list()
   #--------------------------------------------
   result$formula <- object$formula
   result$formulaStr <- object$formulaStr
   result$dataFormat <- if( object$wideFormat ) "wide" else "narrow"
   result$freqSupplied <- object$freqSupplied
   result$weightSupplied <- object$weightSupplied
   result$dataFreqTot <- object$dataFreqTot
   result$weightTot <- object$weightTot
   result$weightScaleFac <- object$weightScaleFac
   result$nrowInputData <- object$dimVec[ "nrowInputData" ]
   result$nLevels <- object$dimVec[ "nLevels" ]
   result$nCovPatt <- object$dimVec[ "nCovPatt" ]
   result$nCovPattEmpty <- sum( object$emptyCovPatt )
   result$nDataPatt <- object$dimVec[ "nDataPatt" ]
   result$nParamThisModel <- object$dimVec[ "nParamThisModel" ]
   result$nParamSatModel <- object$dimVec[ "nParamSatModel" ]
   result$nParamEstimated <- if( object$saturated ) 
      result$nParamSatModel else result$nParamThisModel
   result$df.residual <- if( object$saturated ) 0L else
      result$nParamSatModel - result$nParamThisModel
   names(result$df.residual) <- NULL
   dispersion <- as.double( dispersion )[1L]
   stopifnot( dispersion > 0 )
   result$dispersion <- dispersion
   #--------------------------------------------
   result$responseLevels <- object$responseLevels
   result$baselineInt <- object$baselineInt
   result$baselineLev <- result$responseLevels[ result$baselineInt ]
   result$predNames <- rownames( object$coefficients )
   result$ident <- all( diag( object$dimVec["nLevels"] ) == object$pertMat )
   #--------------------------------------------
   result$prior <- object$prior
   result$priorFreqTot <- object$priorFreqTot
   result$priorFreqTotSupplied <- object$priorFreqTotSupplied
   result$priorAlloc <- object$priorAlloc
   result$priorAllocSupplied <- object$priorAllocSupplied
   result$saturated <- object$saturated
   #--------------------------------------------
   result$startValSource <- object$startValSource
   #--------------------------------------------
   result$method <- object$method
   result$control <- object$control
   result$aborted <- object$aborted
   result$vhatFailed <- object$vhatFailed
   result$DAPFailed <- object$DAPFailed
   if( object$method %in% c("EM","NR","FS") ) {
      result$iter <- object$iter
      result$converged <- object$converged
      result$boundary <- object$boundary
      result$lenGrad <- if( object$saturated ) NULL else
         sqrt( sum( object$score^2 ) )
   } else if( object$method == "MCMC" ) {
      result$iter <- object$nIterActual
      result$discarded <- if( result$iter >= object$control$burnMCMC ) 
         object$control$burnMCMC else ( result$iter - object$control$burnMCMC )
      result$afterBurnIn <- max( result$iter - object$control$burnMCMC, 0 )
      result$thin <- object$control$thinMCMC
      result$imputeEvery <- object$control$imputeEvery
      result$nSampleActual <- object$nSampleActual
      result$nImpActual <- object$nImpActual
      result$mhAcceptRate <- if(object$saturated) NULL else
         object$mhAcceptRate
   } else if( object$method == "approxBayes" ) {
      result$iter <- object$nIterActual
      result$nSampleActual <- object$nSampleActual
      result$nImpActual <- object$nImpActual
   }
   #--------------------------------------------
   result$loglik <- object$loglik
   result$logprior <- object$logprior
   result$logP <- object$logP
   #--------------------------------------------
   result$coefficients <- object$coefficients
   result$coefVec <- object$coefVec
   result$score <- object$score
   result$hess <- object$hess
   result$vhatCoefVec <- object$vhatCoefVec
   #--------------------------------------------
   result$margProportionsYstar <- object$margProportionsYstar
   result$margProportionsY <- object$margProportionsY
   result$pertMat <- object$pertMat
   #--------------------------------------------
   result$showCoef <- as.logical( showCoef )[1L]
   result$digits <- as.integer( digits ) [1L]
   #--------------------------------------------
   if( ! object$saturated ) {
      if( object$method %in% c("EM","NR","FS") ) { 
         result$headerCoef <-
            gettext("Estimated coefficients with Hessian-based SEs" )
         result$headerCoef <- paste( result$headerCoef,
            gettextf( "  (dispersion parameter taken to be %f)", dispersion),
            sep="\n" )
         coefArray <- array( numeric(),
            c( NROW(object$coefficients), 4L, NCOL(object$coefficients) ) )
         dimnames(coefArray) <- list( rownames(object$coefficients),
               c("coef", "SE", "zstat", "pval" ), 
               paste( "Response =", colnames( object$coefficients ) ) )
         coefArray[,"coef",] <- object$coefficients
         if( ! object$vhatFailed ) {
            coefArray[,"SE",] <- .unstackCoef(
               sqrt( dispersion * diag(object$vhatCoef) ),
               result$predNames, result$responseLevels, result$baselineLev )
            zstat <- coefArray[,"coef",] / coefArray[,"SE",]
            pval <- 2 * pnorm( - abs(zstat) )
            coefArray[,"zstat",] <- round(zstat, 2)
            coefArray[,"pval",] <- round(pval, 4)
         }
      } else if( object$method %in% c("MCMC","approxBayes") ) {
         nDraws <- max( object$nIterActual -  object$control$burnMCMC, 0L )
         if( nDraws == 0L ) {
            result$headerCoef <- gettext(
               "No MCMC samples available for estimating coefficients" )
            coefArray <- NULL
	 } else {
	    result$headerCoef <- gettextf(
               "Direct estimates and SE's based on %i MCMC samples",
               nDraws ) 
            coefArray <- array( numeric(),
               c( NROW(object$coefficients), 4L, NCOL(object$coefficients) ) )
            dimnames(coefArray) <- list( rownames(object$coefficients),
                  c("coef", "SE", "zstat", "pval" ), 
                  paste( "Response =", colnames( object$coefficients ) ) )
            coefArray[,"coef",] <- object$coefMCMC
            coefArray[,"SE",] <- .unstackCoef(
               sqrt( diag(object$vhatCoefMCMC) ),
               result$predNames, result$responseLevels, result$baselineLev )
            zstat <- coefArray[,"coef",] / coefArray[,"SE",]
            pval <- 2 * pnorm( - abs(zstat) )
            coefArray[,"zstat",] <- round(zstat, 2)
            coefArray[,"pval",] <- round(pval, 4)
         }
      }
      result$coefficients <- coefArray
   } else result$coefficients <- result$headerCoef <- NULL
   #--------------------------------------------
   if( object$method %in% c("EM", "NR", "FS") ) {
      result$headerPiMarg <-
         result$headerPiMarg <- gettext(
            "Estimated marginal probs. for true response" )
         result$piMarg <- object$fittedPiMarg
   } else {
      nDraws <- max( object$nIterActual -  object$control$burnMCMC, 0L )
      if( nDraws == 0L ) {
         result$headerPiMarg <- gettext(
            "Marginal probs. for true response based on final parameters" )
         result$piMarg <- object$fittedPiMarg
      } else {
         result$headerPiMarg <- gettextf(
           "Marginal probs. for true response estimated from %i MCMC samples",
            nDraws )
         result$piMarg <- object$piMargMCMC
      }
   }
   #--------------------------------------------
   structure( result, class = "summary.rrLogit" )
}
   
print.rrLogit <- function(x, ...) {
   stopifnot( inherits(x, "rrLogit") )
   print( summary(x, ...) )
   invisible()
}

print.summary.rrLogit <- function(x, digits=x$digits, ...) {
   stopifnot( inherits(x, "summary.rrLogit") )
   #--------------------------------------------
   cat(x$formulaStr, sep="\n")
   strA <- format( c("Saturated option:", "Prior:"), justify="right")
   strB <- format( c( as.character(x$saturated), x$prior), justify="left")
   cat( paste(strA, strB), sep="\n" )
   cat("\n") 
   #--------------------------------------------
   strA <- format( c( "Data format:", 
      "Frequencies supplied:", "Rows of supplied data:", 
      "Total N in supplied data:", "Distinct covariate patterns:",
      "Empty covariate patterns:"), justify="right" )
   strB <- format( c( format(x$dataFormat),
      format( as.character(x$freqSupplied) ), 
      format(x$nrowInputData), format(x$dataFreqTot), format(x$nCovPatt),
      format(x$nCovPattEmpty) ), justify="left")
   cat( paste(strA, strB), sep="\n")
   cat("\n")
   #--------------------------------------------
   strA <- format( c("Response categories:",
      "Baseline category:"),
      justify="right")
   strB <- format( c( paste(x$responseLevels, collapse=" "),
      x$baselineLev),
      justify="left")
   cat( paste(strA, strB), sep="\n" )
   cat("\n") 
   #--------------------------------------------
   if( x$ident ) {
      cat("Assumed perturbation matrix: Identity", sep="\n")
   } else {
      cat("Assumed perturbation matrix", sep="\n")
      print( x$pertMat, digits=x$digits )
   }
   cat("\n")
   #--------------------------------------------
   strA <- format( c("Number of estimated parameters =",
      "Degrees of freedom ="), justify="right" )
   strB <- format( c(x$nParamEstimated, x$df.residual), justify = "left" )
   cat( paste(strA, strB), sep="\n")
   cat("\n")
   #--------------------------------------------
   if( x$prior == "DAP" ) {
      cat("Data-augmentation prior (DAP)", sep="\n")
      strA <- format( c( "Prior effective sample size =",
         "Prior N per pattern ="), justify="right" )
      strB <-  c( format( c( x$priorFreqTot,
            round( x$priorFreqTot / ( x$nCovPatt - x$nCovPattEmpty ),
            digits=digits) ), justify = "left" ) )
      cat( paste(strA, strB), sep="\n")
      cat("Proportions for allocating prior counts:", sep="\n")
      strA <- paste( "   ", format( x$responseLevels, justify="right"),
         sep = "" )
      strB <- format( x$priorAlloc, digits=digits, justify="left")
      cat( paste(strA, strB), sep="\n")
      cat("\n")
   }
   if( x$DAPFailed ) {
      cat( "Negative estimated probability encountered;", sep="\n")
      cat( "model fitting not attempted", sep="\n")
      return( invisible() )
   }
   #--------------------------------------------
   if( x$method %in% c("EM","NR","FS") ) {
      if( x$method == "EM" )
         cat( "Expectation-Maximization (EM) algorithm", sep="\n" )
      if( x$method == "NR" )
         cat( "Newton-Raphson procedure", sep="\n" )
      if( x$method == "FS" )
         cat( "Fisher scoring procedure", sep="\n" )
      cat( paste("Starting values:", x$startValSource), sep="\n" )
      if( x$control$startValJitter > 0 ) {
         cat("Starting values jittered with Gaussian noise", sep="\n")
	 cat( gettextf( "Jitter SD = %f", x$control$startValJitter ),
	    sep="\n" )
      }
      if( x$converged ) {
         cat( gettextf( "Converged at iteration %i", x$iter ), sep="\n")
      } else {
         cat( gettextf( "Failed to converge by iteration %i",
            x$iter ), sep="\n")
      }
      if( !is.null(x$lenGrad) ) 
          cat( gettextf( "Gradient length = %f", x$lenGrad ), sep="\n" )
      if( x$boundary ) cat("Estimate at or near boundary", sep="\n" )
      cat("\n")
      strA <- format( c("Final logP =", "Final loglik ="), justify="right")
      strB <- format( c(x$logP, x$loglik), justify="left")
      cat( paste(strA, strB), sep="\n")
      cat("\n")
      if( ! x$saturated ) {
         if( x$vhatFailed ) {
            cat("Hessian-based SEs unavailable", sep="\n")
            cat("\n")
	 }
      }
   } else {
      if( x$method == "MCMC" ) {
         if( x$saturated ) {
	    pStr <- "MCMC: Data augmentation (DA) for saturated model"
            cat( pStr, sep="\n")
	 } else {
	    if( x$control$typeMCMC == "RWM" ) {
	       pStr <- "MCMC: Random-walk Metropolis"
               cat( pStr, sep="\n")
               cat("\n")
               cat("Tuning parameters:", sep="\n")
               strA <- format( c("proposal df =", 
                  "scale factor =", "stuck limit ="), justify="right")
               strB <- c( format( x$control$tuneRWM[1L] ), 
                  format( x$control$tuneRWM[2L] ),
		  format( x$control$stuckLimit) )
               strB <- format( strB, justify="left")
               cat( paste(strA, strB), sep="\n")
               cat("\n")
               strA <- "Accept rate ="
               strB <- format(x$mhAcceptRate)
               cat( paste(strA, strB), sep="\n")
               cat("\n")
	    } else if( x$control$typeMCMC == "DA" ) {
               pStr <- "MCMC: Data augumentation (DA) with Metropolis-Hastings"
               cat( pStr, sep="\n")
               cat("\n")
               cat("Tuning parameters:", sep="\n")
               strA <- format( c("proposal df =", 
                  "step size =", "scale factor =", "stuck limit ="),
		  justify="right")
               strB <- c( format( x$control$tuneDA[1L] ), 
                  format( x$control$tuneDA[-1L] ),
		  format( x$control$stuckLimit ) )
               strB <- format( strB, justify="left")
               cat( paste(strA, strB), sep="\n")
	       cat("\n")
               strA <- "Accept rate ="
               strB <- format(x$mhAcceptRate)
               cat( paste(strA, strB), sep="\n")
               cat("\n")
	    }
	 }
         strA <- format( c("Iterations performed =", 
            "Iterations discarded as burn-in =", 
            "Iterations after burn-in =", 
            "Thinning interval for saved series =",
            "Samples in saved series =",
            "Imputation interval =",
            "Number of imputations stored ="), justify="right")
         strB <- format( c( x$iter,
            x$discarded,
            x$afterBurnIn,
            x$thin,
            x$nSampleActual,
            x$imputeEvery,
            x$nImpActual), justify="left")
         cat( paste(strA, strB), sep="\n")
         cat("\n")
      }
   }
   #--------------------------------------------
   if( ( ! x$saturated ) & x$showCoef ) {
      cat( x$headerCoef, sep="\n")
      print(x$coefficients, digits=x$digits, ...)
      cat("\n")
   }
   #--------------------------------------------
   cat( x$headerPiMarg, sep="\n" )
   print( round(x$piMarg, digits=x$digits), ...)
   #--------------------------------------------
   invisible()
}

loglik <- function(obj) {
   stopifnot( inherits(obj, "rrLogit") )
   obj$loglik
}

logP <- function(obj) {
   stopifnot( inherits(obj, "rrLogit") )
   obj$logP
}

minus2logPSeries <- function( obj, startValShift = TRUE, coda = TRUE ){
   stopifnot( inherits(obj, "rrLogit") )
   if( obj$method %in% c("EM","NR","FS") ) {
      message( gettextf(
         "No series are produced when method = '%s'", obj$method ),
         domain = NA )
      return( invisible() )
   }
   if( obj$nSampleActual < 1L ) {
      message( gettext(
         "Insufficient itertions after burn-in, no series is available" ),
         domain = NA )
      return( invisible() )
   }
   shift <- if( startValShift && obj$beganAtMode ) obj$startLogP else 0
   result <- -2 * ( obj$logPSeries[1:obj$nSampleActual] - shift )
   if(coda) coda::mcmc( result, thin=obj$control$thinMCMC ) else result
}

coefSeries <- function( obj, coda = TRUE ){
   stopifnot( inherits(obj, "rrLogit") )
   if( obj$saturated ) {
      warning( gettext(
         "Model is saturated, no coefficients are present"),
	 domain = NA )
      return( invisible() )
   }
   if( obj$method %in% c("EM","NR","FS") ) {
      message( gettextf(
         "No series are produced when method = '%s'", obj$method ),
         domain = NA )
      return( invisible() )
   }
   if( obj$nSampleActual < 1L ) {
      message( gettext(
         "Insufficient itertions after burn-in, no series is available" ),
         domain = NA )
      return( invisible() )
   }
   result <- obj$coefVecSeries[1:obj$nSampleActual,]
   if(coda) coda::mcmc( result, thin=obj$control$thinMCMC ) else result
}

margProbSeries <- function( obj, noisy = FALSE, coda = TRUE ){
   noisy <- as.logical(noisy)[1L]
   if( obj$method %in% c("EM","NR","FS") ) {
      message( gettextf(
         "No series are produced when method = '%s'", obj$method ),
         domain = NA )
      return( invisible() )
   }
   if( obj$nSampleActual < 1L ) {
      message( gettext(
         "Insufficient itertions after burn-in, no series is available" ),
         domain = NA )
      return( invisible() )
   }
   result <- obj$piMargSeries[1:obj$nSampleActual,]
   if( noisy ) result <- result %*% obj$pertMat
   if(coda) coda::mcmc( result, thin=obj$control$thinMCMC ) else result
}

impList <- function( obj ){
   stopifnot( inherits(obj, "rrLogit") )
   if( obj$nImpActual <= 0L ) {
      message( gettext("Object contains no imputations"), domain = NA )
      return( invisible( list() ) )
   }
   result <- as.list( 1:obj$nImpActual )
   if( obj$microData ) {
      for( m in 1:obj$nImpActual ) {
         tmp <- factor( obj$impVecSeries[,m,drop=TRUE],
	    levels=1:obj$dimVec["nLevels"] )
         levels(tmp) <- obj$responseLevels
         result[[m]] <- tmp
      }
   } else {
      for(m in 1:obj$nImpActual) result[[m]] <-
         obj$impMatSeries[,,m,drop=TRUE]
   }
   return(result)
}

df.residual.rrLogit <- function( object, ...){
   stopifnot( inherits(object, "rrLogit") )
   summary(object)$df.residual
}

coef.rrLogit <- function( object, covMat = FALSE, meanSeries = TRUE,
   completeData = FALSE, ...) {
   stopifnot( inherits( object, "rrLogit" ) )
   if( object$saturated ) return( NULL )
   covMat <- as.logical(covMat)[1L]
   meanSeries <- as.logical(meanSeries)[1L]
   if( covMat & completeData & ( object$method != "EM" ) ) stop( gettext(
      'completeData = TRUE may only be used when method="EM"'),
      domain = NA )
   if( object$method %in% c("EM","NR","FS") ) {
      coef <- object$coefficients
      coefVec <- object$coefVec
      if( object$vhatFailed & covMat ) {
         warning( gettext("Estimated covariance matrix not available"),
	    domain = NA )
	 vhatCoef <- NULL
      } else {
         vhatCoef <- if( completeData ) solve( -object$hessA ) else
	    object$vhatCoef
      }
   } else {
      nDraws <- max( object$nIterActual -  object$control$burnMCMC, 0L )
      if( nDraws == 0L ) {
         warning( gettext(
	    "No MCMC samples available for estimating coefficients" ),
            domain = NA )
	 return( invisible() )
      }
      if( meanSeries ) {
         coef <- object$coefMCMC
	 coefVec <- object$coefVecMCMC
         vhatCoef <- object$vhatCoefMCMC
      } else {
         coef <- object$coefficients
	 coefVec <- object$coefVec
	 if( covMat ) {
	    warning( gettext(
	       "Covariance matrix not available when 'meanSeries' is FALSE"),
	       domain = NA )
	    vhatCoef <- NULL   
         }
      }
   }
   if( covMat ) {
      return( list( coef = coefVec, covMat = vhatCoef ) )
   } else {
      return( coef )
   }
}

fitted.rrLogit <- function( object, type=c("prob", "link", "mean"),
   noisy=FALSE, covPatt=FALSE, include.predvars=FALSE, digits=NULL, ...) {
   stopifnot( inherits( object, "rrLogit" ) )
   type <- match.arg(type)
   noisy <- as.logical(noisy)[[1L]]
   covPatt <- as.logical(covPatt)[[1L]]
   if( object$method %in% c("EM","NR","FS") ) {
      fittedPi <- object$fittedPi
      fittedPistar <- object$fittedPistar
      coefficients <- if( ! object$saturated ) object$coefficients else NULL
   } else if( object$method %in% c("MCMC","approxBayes") ) {
      nDraws <- max( object$nIterActual -  object$control$burnMCMC, 0L )
      if( nDraws == 0L ) {
         warning( gettext(
	    "No MCMC samples available for getting fitted values" ),
            domain = NA )
	 return( invisible() )
      }
      fittedPi <- object$fittedPiMCMC
      fittedPistar <- object$fittedPistarMCMC
      coefficients <- if( ! object$saturated ) object$coefMCMC else NULL
   } else {
      stop( gettextf( "Method '%s' not recognized", object$method ),
         domain = NA )
   }
   if( type == "prob" ) {
      result <- if( noisy ) fittedPistar else fittedPi 
   } else if( type == "link" ) {
      if( noisy ) stop( gettext(
         'type="link" is not allowed when "noisy=TRUE"' ), domain = NA )  
      if( object$saturated ) {
         result <- fittedPi
         result <- log( result / result[, object$baselineInt] )
      } else {
         result <- object$modelMatrix %*% coefficients
      }      
   } else if( type == "mean" ) {
      result <- if( noisy ) fittedPistar else fittedPi 
      result <- result * object$freqForCovPatt
   }
   if( ! is.null(digits) ) result <- round( result, digits=digits )
   if( include.predvars ){
      tmp <- model.frame(object)
      w <- (1:NROW(tmp))[ ! duplicated( object$covPatt ) ]
      tmp <- tmp[w,]
      tmp <- tmp[ -attr(terms(object), "response") ]
      tmp$`(freq)` <- tmp$`(weight)` <- NULL
      result <- cbind(tmp, result)
   }
   rownames(result) <- NULL
   if( covPatt ) attr(result, "covPatt") <- object$covPatt
   return( result )
}

# residuals are computed with frequencies, ignoring any survey weights

residuals.rrLogit <- function( object, type=c("pearson","response"),
      dispersion=1.0, covPatt=FALSE, include.predvars=FALSE, digits=NULL, ...) {
   stopifnot( inherits( object, "rrLogit" ) )
   type <- match.arg(type)
   covPatt <- as.logical(covPatt)[[1L]]
   if( object$method %in% c("EM","NR","FS") ) {
      fittedPistar <- object$fittedPistar
      dispersion <- as.double(dispersion)[1L]
      stopifnot( dispersion > 0 )
   } else if( object$method %in% c("MCMC","approxBayes") ) {
      nDraws <- max( object$nIterActual -  object$control$burnMCMC, 0L )
      if( nDraws == 0L ) {
         warning( gettext(
	    "No MCMC samples available for getting estimated expected values" ),
            domain = NA )
	 return( invisible() )
      }
      fittedPistar <- object$fittedPistarMCMC
      dispersion <- 1.0
   } else {
      stop( gettextf( "Method '%s' not recognized", object$method ),
         domain = NA )
   }
   observed <- object$fstarMat
   expected <- fittedPistar * apply(object$fstarMat, 1, sum)
   raw <- observed - expected
   if( type == "pearson" ) {
      ident <- all( diag( object$dimVec["nLevels"] ) == object$pertMat )
      if( object$boundary & ( ! ident ) ) warning( gettext(
         "Estimate at or near boundary, Pearson residuals may be unstable" ),
	 domain = NA )
      result <- raw
      result[ raw!=0 ] <- result[ raw!=0 ] / sqrt( expected[ raw!=0 ] )
      result <- result / sqrt(dispersion)
   } else if( type == "response" ) {
      result <- raw
   }
   if( ! is.null(digits) ) result <- round( result, digits=digits )
   if( include.predvars ){
      tmp <- model.frame(object)
      w <- (1:NROW(tmp))[ ! duplicated( object$covPatt ) ]
      tmp <- tmp[w,]
      tmp <- tmp[ -attr(terms(object), "response") ]
      tmp$`(freq)` <- tmp$`(weight)` <- NULL
      result <- cbind(tmp, result)
   }
   rownames(result) <- NULL
   if( covPatt ) attr(result, "covPatt") <- object$covPatt
   return(result)
}

model.frame.rrLogit <- function( formula, ...) {
   stopifnot( inherits( formula, "rrLogit" ) )
   formula$mF
}

model.matrix.rrLogit <- function( object, ...) {
   stopifnot( inherits( object, "rrLogit" ) )
   result <- object$modelMatrix[ object$covPatt, ]
   rownames(result) <- object$rowNamesInputData
   return(result)
}

terms.rrLogit <- function(x, ...) {
   stopifnot( inherits( x, "rrLogit" ) )
   attr(x$mF, "terms")
}

formula.rrLogit <- function(x, ...) {
   stopifnot( inherits(x, "rrLogit") )
   form <- x$formula
   environment(form) <- environment(x$formula)
   form 
}

anova.rrLogit <-
   function( object, ..., method=c("lrt", "logP", "AIC", "BIC"),
      pval = FALSE, pvalDigits=4L, showRank=NULL ) {
   method <- match.arg(method)
   dotargs <- list(...)
   named <- if (is.null(names(dotargs))) 
        rep_len(FALSE, length(dotargs))
   else (names(dotargs) != "")
   if (any(named)) warning(
  "the following arguments to 'anova.rrLogit' are invalid and dropped: ", 
      paste(deparse(dotargs[named]), collapse = ", ") )
   dotargs <- dotargs[!named]
   modList <- c( list(object), dotargs )
   if( length(modList) < 2L ) stop( gettext(
      'Need at least two objects of class "rrLogit" to compare'),
      domain = NA ) 
   is.rrLogit <-
      vapply(modList, function(x) inherits(x, "rrLogit"), NA)
   if( any( !is.rrLogit ) ) stop( gettext(
      'Some supplied objects are not of class "rrLogit"'), domain = NA ) 
   summList <- lapply( modList, summary.rrLogit )
   responseLevels <- lapply( summList, `[[`, "responseLevels" )
   sameLevels <- 
      unlist( lapply( responseLevels, function(x,y) isTRUE(all.equal(x,y)),
      y=responseLevels[[1]] ) )
   if( ! all(sameLevels) ) stop( gettext(
      'Fitted models do not all have the same response levels'),
      domain = NA )
   methods <- unlist( lapply( summList, `[[`, "method" ) )
   if( any( methods == "MCMC" ) ) stop( gettext(
      'anova not defined for models fit with method="MCMC"'),
      domain = NA )
   if( any( methods == "approxBayes" ) ) stop( gettext(
      'anova not defined for models fit with method="approxBayes"'),
      domain = NA )
   priorTypes <- unlist( lapply( summList, `[[`, "prior" ) )
   if( ! all( priorTypes == priorTypes[1L] ) ) warning( gettext(
      'Fitted models do not all have the same prior distribution'),
      domain = NA )
   if( priorTypes[1L] == "DAP" ) {
      priorFreqs <- unlist( lapply( summList, `[[`, "priorFreqTot" ) )
      if( ! all( priorFreqs == priorFreqs[1L] ) ) warning( gettext(
         'Fitted models do not all have the same prior sample size'),
         domain = NA )
   }
   nTotal <- unlist( lapply( summList, `[[`, "dataFreqTot" ) )
   if( ! all( nTotal == nTotal[1L] ) ) warning( gettext(
      'Fitted models are based on different sample sizes'),
      domain = NA )
   #----------------------------------
   formulaStr <- unlist( lapply( summList, `[[`, "formulaStr" ) )
   formulaStr <- paste("Model ", format(1:length(summList)), ": ",
      formulaStr, sep="")
   saturated <- unlist( lapply( summList, `[[`, "saturated" ) )
   formulaStr[saturated] <- paste( formulaStr[saturated], "(saturated)" )
   formulaStr <- paste( formulaStr, collapse="\n" )
   nParams <- unlist( lapply( summList, `[[`, "nParamEstimated" ) )
   resid.df <- unlist( lapply( summList, `[[`, "df.residual" ) )
   if( method %in% c("lrt", "logP") ) {
      meas <- if( method == "lrt" ) 
         unlist( lapply( summList, `[[`, "loglik" ) ) else
         unlist( lapply( summList, `[[`, "logP" ) )
      meas <- -2*meas
      result <- data.frame( nParams, meas )
      rownames(result) <- NULL
      names(result)[2L] <- if( method == "lrt" ) "-2*loglik" else "-2*logP" 
      result$df <- - ( c( NA, nParams[-length(nParams)]) - nParams )
      result$change <- c( NA, meas[-length(meas)] ) - meas
      pvalDigits <- as.integer(pvalDigits)[1L]
      if( pval ) result$pval <- 
         round( 1 - pchisq( result$change, result$df ), pvalDigits )
   } else {
      meas <- -2 * unlist( lapply( summList, `[[`, "loglik" ) )
      result <- data.frame( nParams, meas )
      rownames(result) <- NULL
      names(result)[2] <- "-2*loglik" 
      IC <- if( method == "AIC" ) meas + 2*nParams else
         meas + log(nTotal) * nParams
      if( method == "AIC" ) meas <- result$AIC <- IC else
         meas <- result$BIC <- IC
   }
   showRank <- if( is.null(showRank) ) method %in% c("AIC", "BIC") else
      as.logical(showRank)[1L] 
   if( showRank ) result$rank <- rank(meas)
   structure( result,
      heading = formulaStr,
      class = c("anova", "data.frame") )
}

predict.rrLogit <- function(object, newdata=NULL, freq=NULL,
   type=c("prob", "link", "mean"), noisy=FALSE, se.fit=FALSE,
   na.action=na.pass, completeData=FALSE, ...) {
   stopifnot( inherits(object, "rrLogit") )
   #--------------------------------------------
   if( completeData & ( object$method != "EM" ) ) stop( gettext(
      'completeData = TRUE may only be used when method="EM"'),
      domain = NA )
   #--------------------------------------------
   type <- match.arg(type)
   noisy <- as.logical(noisy)[1L]
   se.fit <- as.logical(se.fit)[1L]
   if( object$saturated & se.fit ) {
      warning( gettext(
      "Standard errors not available, because model has 'saturated=TRUE'"),
         domain = NA )
      se.fit <- FALSE
   }
   if( ( ! object$saturated ) && ( object$method %in% c("EM","NR","FS") )
      && se.fit && object$vhatFailed ) {
      warning( gettext(
      "Standard errors not available, because logP is not concave"),
         domain = NA )
      se.fit <- FALSE
   }
   if( object$method %in% c("MCMC","approxBayes") ) {
      nDraws <- max( object$nIterActual -  object$control$burnMCMC, 0L )
      if( nDraws == 0L ) stop( gettext(
         "No MCMC samples available for getting estimated expected values" ),
         domain = NA )
   }
   typeInt <- match(type, c("prob","link","mean") )
   noisyInt <- as.integer(noisy)
   seFitInt <- as.integer(se.fit)
   #--------------------------------------------
   mvcode <- object$specialCodes["mvcode"]
   pertMat <- object$pertMat
   pertMatInv <- object$pertMatInv
   ident <- all( diag( object$dimVec["nLevels"] ) == pertMat )
   coefficients <- if( object$method %in% c("EM","NR","FS") )
      object$coefficients else object$coefMCMC
   vhatCoef <- if( object$method %in% c("EM","NR","FS") )
      object$vhatCoef else object$vhatCoefMCMC
   if( object$saturated ) {
      coefficients[] <- mvcode
      vhatCoef[] <- mvcode
   }
   if( ( ! object$saturated ) & ( object$method == "EM" ) & completeData ) {
      vhatCoef <- solve( - object$hessA )
   }
   #--------------------------------------------
   if( is.null( newdata ) ) {
      wideFormat <- object$wideFormat
      freqSupplied <- object$freqSupplied
      dimVec <- object$dimVec
      modelMatrix <- object$modelMatrix
      fitweightRowInputData <- object$fitweightRowInputData
      fitweightCovPatt <- object$fitweightCovPatt
      fitweightDataPatt <- object$fitweightDataPatt
      freqForDataPatt <- object$freqForDataPatt
      covPatt <- object$covPatt
      dataPatt <- if( object$wideFormat ) integer() else object$dataPatt
      covPattForDataPatt <- object$covPattForDataPatt
      responseForDataPatt <- object$responseForDataPatt
      pertMat <- object$pertMat
      pertMatInv <- object$pertMatInv
      fittedPi <- if( object$method %in% c("EM","NR","FS") )
         object$fittedPi else object$fittedPiMCMC
      fittedPi[ is.na(fittedPi) ] <- mvcode # NAs for empty cov patts
      fittedPistar <- if( object$method %in% c("EM","NR","FS") )
         object$fittedPistar else object$fittedPistarMCMC
      fittedPistar[ is.na(fittedPistar) ] <- mvcode # NAs for empty cov patts
      if( ! object$saturated ) {
         # these will be recomputed
         fittedPi[] <- mvcode
	 fittedPistar[] <- mvcode
      }
      freqRowInputData <- object$freqRowInputData
      fitted <- matrix( numeric(1L),
         object$dimVec["nrowInputData"], object$dimVec["nLevels"] )
      rownames(fitted) <- rownames(object$mF)
      colnames(fitted) <- object$responseLevels
      if(se.fit) {
         seMat <- fitted
	 vhatFittedArray <- array( numeric(1L),
	    c( object$dimVec["nrowInputData"], object$dimVec["nLevels"],
	    object$dimVec["nLevels"] ) )
	 dimnames(vhatFittedArray) <- list( rownames(object$mF),
	    object$responseLevels, object$responseLevels )
      } else {
         seMat <- matrix( numeric(), 0L, 0L)
         vhatFittedArray <- array( numeric(), c(0L,0L,0L) )
      }
   } else {
      stopifnot( is.data.frame(newdata) )
      tt <- terms(object)
      Terms <- delete.response(tt)
      mc <- match.call( expand.dots=FALSE )
      mc[[1L]] <- quote( stats::model.frame )
      m <- match( c("object", "newdata", "freq", "na.action" ),
         names(mc), nomatch=0L )
      mc <- mc[ c(1L,m) ]
      names(mc)[2L] <- "formula"
      mc[[2L]] <- Terms
      names(mc)[3L] <- "data"
      mc$drop.unused.levels <- FALSE
      mc$xlev <- object$xlevels
      if( is.null(mc$na.action) ) mc$na.action <- quote(na.pass)
      mf <- eval( mc, parent.frame() )
      if(!is.null(cl <- attr(Terms, "dataClasses"))) 
         .checkMFClasses(cl, mf, ordNotOK=TRUE)
      if( is.null( mf$`(freq)` ) ) {
         if( object$freqSupplied & (! object$surveyMode) &
            type == "mean" ) warning( gettext(
            "No 'freq' supplied with newdata, frequencies assumed to be one" ),
            domain = NA )
	 freqRowInputData <- rep(1, NROW(mf))
      } else {
         freqRowInputData <- mf$`(freq)`
	 mf$`(freq)` <- NULL
         if( any(is.na(freqRowInputData)) ) stop( gettext(
            "Missing values in 'freq' are not allowed"), domain = NA )
	 if( any( freqRowInputData < 0 ) )  stop( gettext(
            "Negative values in 'freq' are not allowed" ), domain = NA )
      }
      storage.mode(freqRowInputData) <- "double"
      # identify unique covariate patterns in newdata
      x <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)
      x[ is.na(x) ] <- mvcode
      xdf <- as.data.frame(x)
      cx <- do.call( "paste", c( structure( xdf[,,drop=FALSE], names=NULL ),
         sep="\r" ) )
      uniqueRowNumbers.x <- (1:NROW(x))[ ! duplicated(cx) ]
      cxUnique <- cx[ uniqueRowNumbers.x ]
      nCovPatt <- length(cxUnique)
      covPatt <- match(cx, cxUnique)
      names(covPatt) <- NULL
      modelMatrix <- x[uniqueRowNumbers.x,,drop=FALSE]
      # dummy data patterns
      nDataPatt <- nCovPatt
      dataPatt <- covPatt
      covPattForDataPatt <- 1:nDataPatt
      responseForDataPatt <- rep(1L,nDataPatt)
      # aggregate frequencies
      dFtmp <- data.frame( freqRowInputData, covPatt, dataPatt )
      names(dFtmp) <- c( "freq", "covPatt", "dataPatt" )
      aggF <- aggregate( freq ~ covPatt, data=dFtmp, FUN=sum )
      aggF <- aggF[ order(aggF$covPatt), ]
      freqForCovPatt <- aggF$freq
      names(freqForCovPatt) <- NULL
      storage.mode(freqForCovPatt) <- "double"
      freqForDataPatt <- freqForCovPatt
      # handle fittedPi and fittedPistar
      if( object$saturated ) {
         mmOrig <- as.data.frame( object$modelMatrix )
	 mmOrig[ is.na(mmOrig) ] <- mvcode
	 cxOrig <- do.call( "paste",
	    c( structure( mmOrig[,,drop=FALSE], names=NULL ),
            sep="\r" ) )
	 m <- match( cxUnique, cxOrig )
	 fittedPi <- if( object$method %in% c("EM","NR","FS") )
	    object$fittedPi[m,] else object$fittedPiMCMC[m,]
	 fittedPi[ is.na(fittedPi) ] <- mvcode
	 fittedPistar <- if( object$method %in% c("EM","NR","FS") )
	    object$fittedPistar[m,] else object$fittedPistarMCMC[m,]
	 fittedPistar[ is.na(fittedPistar) ] <- mvcode
	 rownames(fittedPi) <- rownames(fittedPistar) <- NULL
      } else {
         # these will be recomputed
         fittedPi <- matrix( numeric(), nCovPatt, object$dimVec["nLevels"] )
	 rownames(fittedPi) <- NULL
	 colnames(fittedPi) <- object$responseLevels
	 fittedPistar <- fittedPi
      }
      # objects to pass to Fortran
      dimVec <- c( nrowInputData = NROW(mf),
         nLevels = as.integer(object$dimVec["nLevels"]),
         nCovPatt = nCovPatt,
	 nDataPatt = nDataPatt,
         ncolMM = NCOL(modelMatrix),
         nParamThisModel = as.integer(object$dimVec["nParamThisModel"]),
	 nParamSatModel = as.integer(object$dimVec["nParamSatModel"]),
         wideFormatInt = 0L )
      fitweightRowInputData <- freqRowInputData
      fitweightCovPatt <- freqForCovPatt
      fitweightDataPatt <- freqForDataPatt
      fitted <- matrix( numeric(1L),
         dimVec["nrowInputData"], dimVec["nLevels"] )
      rownames(fitted) <- rownames(mf)
      colnames(fitted) <- object$responseLevels
      if(se.fit) {
         seMat <- fitted
	 vhatFittedArray <- array( numeric(1L),
	    c( dimVec["nrowInputData"], dimVec["nLevels"], dimVec["nLevels"] ) )
	 dimnames(vhatFittedArray) <- list( rownames(mf),
	    object$responseLevels, object$responseLevels )
      } else {
         seMat <- matrix( numeric(), 0L, 0L)
         vhatFittedArray <- array( numeric(), c(0L,0L,0L) )
      }
   }
   #--------------------------------------------
   # needed for dotCall64::.C64
   SIGNATURE <- c(
      # inputs
      dimVec = "integer",
      modelMatrix = "double",
      fitweightRowInputData = "double",
      fitweightCovPatt = "double",
      fitweightDataPatt = "double",
      freqForDataPatt = "double",
      covPatt = "integer",
      dataPatt = "integer",
      covPattForDataPatt = "integer",
      responseForDataPatt = "integer",
      surveyModeInt = "integer",
      baselineInt = "integer",
      pertMat = "double",
      pertMatInv = "double",
      priorInt = "integer",
      priorFreqTot = "double",
      priorAllocSuppliedInt = "integer",
      priorAlloc = "double",
      saturatedInt = "integer",
      methodInt = "integer",
      ctrlInt = "integer",
      ctrlReal = "double",
      specialCodes = "double",
      typeInt = "integer",
      noisyInt = "integer",
      seFitInt = "integer",
      freqRowInputData = "double",
      coefficients = "double",
      vhatCoef = "double",
      # inouts
      fittedPi = "double",
      fittedPistar = "double",
      # outputs
      fitted = "double",
      seMat = "double",
      vhatFittedArray = "double",
      # messaging
      status = "integer",
      msgLenMax = "integer",
      msgCodes = "integer",
      msgLenActual = "integer" )
   INTENT <- c(
      # inputs
      dimVec = "r",
      modelMatrix = "r",
      fitweightRowInputData = "r",
      fitweightCovPatt = "r",
      fitweightDataPatt = "r",
      freqForDataPatt = "r",
      covPatt = "r",
      dataPatt = "r",
      covPattForDataPatt = "r",
      responseForDataPatt = "r",
      surveyModeInt = "r",
      baselineInt = "r",
      pertMat = "r",
      pertMatInv = "r",
      priorInt = "r",
      priorFreqTot = "r",
      priorAllocSuppliedInt = "r",
      priorAlloc = "r",
      saturatedInt = "r",
      methodInt = "r",
      ctrlInt = "r",
      ctrlReal = "r",
      specialCodes = "r",
      typeInt = "r",
      noisyInt = "r",
      seFitInt = "r",
      freqRowInputData = "r",
      coefficients = "r",
      vhatCoef = "r",
      # inouts
      fittedPi = "rw",
      fittedPistar = "rw",
      # outputs
      fitted = "w",
      seMat = "w",
      vhatFittedArray = "w",
      # messaging
      status = "w",
      msgLenMax = "r",
      msgCodes = "w",
      msgLenActual = "w" )
   #--------------------------------------------
   # create a matrix for holding message codes
   msgLenMax <- 40L
   msgCodes <- matrix( 0L, msgLenMax, 17L )
   #--------------------------------------------
   result <- dotCall64::.C64("rrlogit_predict",
      SIGNATURE = SIGNATURE,
      # inputs
      dimVec = dimVec,
      modelMatrix = modelMatrix,
      fitweightRowInputData = fitweightRowInputData,
      fitweightCovPatt = fitweightCovPatt,
      fitweightDataPatt = fitweightDataPatt,
      freqForDataPatt = freqForDataPatt,
      covPatt = covPatt,
      dataPatt = dataPatt,
      covPattForDataPatt = covPattForDataPatt,
      responseForDataPatt = responseForDataPatt,
      surveyModeInt = 0L,
      baselineInt = object$baselineInt,
      pertMat = pertMat,
      pertMatInv = pertMatInv,
      priorInt = object$priorInt,
      priorFreqTot = object$priorFreqTot,
      priorAllocSuppliedInt = object$priorAllocSuppliedInt,
      priorAlloc = object$priorAlloc,
      saturatedInt = object$saturatedInt,
      methodInt = object$methodInt,
      ctrlInt = object$ctrlInt,
      ctrlReal = object$ctrlReal,
      specialCodes = object$specialCodes,
      typeInt = typeInt,
      noisyInt = noisyInt,
      seFitInt = seFitInt,
      freqRowInputData = freqRowInputData,
      coefficients = coefficients,
      vhatCoef = vhatCoef,
      # inouts
      fittedPi = fittedPi,
      fittedPistar = fittedPistar,
      # outputs
      fitted = fitted,
      seMat = seMat,
      vhatFittedArray = vhatFittedArray,
      # messaging
      status = integer(1L),
      msgLenMax = msgLenMax,
      msgCodes = msgCodes,
      msgLenActual = integer(1L),
      # other args to .C64
      NAOK = TRUE,
#      INTENT = INTENT,
      PACKAGE = "rrLogit")
      #--------------------------------------------
      # display message from Fortran, if present
      msg.lines <- .msg( result$msgCodes, result$msgLenActual )
      if( is.null( msg.lines ) ){
         msg <- "OK"
      } else{
         msg <- paste0( msg.lines, collapse="\n" )
      }
      msg <- paste( msg, "\n", sep="")
      if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( result$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   result$fitted[ result$fitted == mvcode ] <- NA
   result$seMat[ result$seMat == mvcode ] <- NA
   result$vhatFittedArray[ result$vhatFittedArray == mvcode ] <- NA
   #--------------------------------------------
   return( if( se.fit ) 
      list( fit=result$fitted, se.fit=result$seMat,
         cov.fit.array=result$vhatFittedArray ) else
      result$fitted )
}

impute <- function( object, ...) {
   # S3 generic function
   UseMethod("impute")
}

impute.default <- function( object, ...) {
   stop( gettext(
      'First argument must be an object of class "rrLogit"'),
      domain = NA )
}

impute.rrLogit <- function(object, newdata=NULL, freq=NULL,
   type=c("random", "condMean"), impVarName=NULL,
   meanSeries = FALSE, na.action=na.pass, ...) {
   stopifnot( inherits(object, "rrLogit") )
   #--------------------------------------------
   type <- match.arg(type)
   typeInt <- match(type, c("random", "condMean") )
   meanSeries <- as.logical(meanSeries)[1L]
   if( object$method %in% c("EM","NR","FS") ) meanSeries <- FALSE
   #--------------------------------------------
   mvcode <- object$specialCodes["mvcode"]
   pertMat <- object$pertMat
   pertMatInv <- object$pertMatInv
   ident <- all( diag( object$dimVec["nLevels"] ) == pertMat )
   coefficients <- object$coefficients
   if( object$saturated ) {
      coefficients[] <- mvcode
   } else {
      if( meanSeries ) {
         coefficients[] <- object$coefMCMC[]
      }
   }
   #--------------------------------------------
   if( is.null( newdata ) ) {
      wideFormat <- object$wideFormat
      microData <- object$microData
      if( microData & ( type=="random" ) ) {
         impVarName <- if( is.null(impVarName) ) 
            paste( object$responseVarName, ".imp", sep="" ) else
	    as.character(impVarName)[1]
	 if( impVarName == object$responseVarName ) impVarName <-
	    paste( object$responseVarName, ".imp", sep="" )
      }
      microDataInt <- as.integer(microData)
      freqSupplied <- object$freqSupplied
      dimVec <- object$dimVec
      modelMatrix <- object$modelMatrix
      fitweightRowInputData <- object$fitweightRowInputData
      fitweightCovPatt <- object$fitweightCovPatt
      fitweightDataPatt <- object$fitweightDataPatt
      freqForDataPatt <- object$freqForDataPatt
      covPatt <- object$covPatt
      dataPatt <- if( object$wideFormat ) integer() else object$dataPatt
      covPattForDataPatt <- object$covPattForDataPatt
      responseForDataPatt <- object$responseForDataPatt
      fittedPi <- if( meanSeries ) object$fittedPiMCMC else object$fittedPi
      fittedPi[ is.na(fittedPi) ] <- mvcode # NAs for empty cov patts
      fittedPistar <- if( meanSeries ) object$fittedPistarMCMC else
         object$fittedPistar
      fittedPistar[ is.na(fittedPistar) ] <- mvcode # NAs for empty cov patts
      if( ! object$saturated ) {
         # these will be recomputed
         fittedPi[] <- mvcode
	 fittedPistar[] <- mvcode
      }
      if( type == "condMean" ) {
         freqRowInputData <- object$freqRowInputData
	 storage.mode(freqRowInputData) <- "double"
	 freqRowInputDataInt <- integer(0L)
	 fMat <- object$fMat
	 storage.mode(fMat) <- "double"
	 fMatInt <- matrix( integer(), 0L, 0L )
      } else {
         freqRowInputDataInt <- as.integer(object$freqRowInputData)
	 if( any( freqRowInputDataInt != object$freqRowInputData ) )
	    warning( gettext(
	    "Some frequencies changed when integerized" ), domain = NA )
         freqRowInputData <- numeric(0L)
	 fMatInt <- object$fMat
	 storage.mode(fMatInt) <- "integer"
	 if( any( fMatInt != object$fMat ) )
	    warning( gettext(
	    "Some frequencies changed when integerized" ), domain = NA )
         fMat <- matrix( numeric(), 0L, 0L )	 
      }
      responseRowInputData <- if( object$wideFormat ) integer(0L) else
         object$responseRowInputData
      if( type == "condMean" ) {
         condMeans <- matrix( numeric(1L),
            object$dimVec["nrowInputData"], object$dimVec["nLevels"] )
         rownames(condMeans) <- rownames(object$mF)
         colnames(condMeans) <- object$responseLevels
	 impMat <- matrix( integer(), 0L, 0L )
	 impVec <- integer(0L)
      } else {
         condMeans <- matrix( numeric(), 0L, 0L )
	 impMat <- matrix( 0L,
	    object$dimVec["nrowInputData"], object$dimVec["nLevels"] )
	 impVec <- integer( object$dimVec["nrowInputData"] )
	 rownames(impMat) <- names(impVec) <- rownames(object$mF)
	 colnames(impMat) <- object$responseLevels
      }
   } else {
      stopifnot( is.data.frame(newdata) )
      Terms <- terms(object)
      mc <- match.call( expand.dots=FALSE )
      mc[[1L]] <- quote( stats::model.frame )
      m <- match( c("object", "newdata", "freq", "na.action" ),
         names(mc), nomatch=0L )
      mc <- mc[ c(1L,m) ]
      names(mc)[2L] <- "formula"
      mc[[2L]] <- Terms
      names(mc)[3L] <- "data"
      mc$drop.unused.levels <- FALSE
      mc$xlev <- object$xlevels
      if( is.null(mc$na.action) ) mc$na.action <- quote(na.pass)
      mf <- eval( mc, parent.frame() )
      if(!is.null(cl <- attr(Terms, "dataClasses"))) 
         .checkMFClasses(cl, mf, ordNotOK=TRUE)
      if( is.null( mf$`(freq)` ) ) {
         if( object$freqSupplied & (! object$wideFormat) ) message( gettext(
            "No 'freq' supplied with newdata, frequencies assumed to be one" ),
            domain = NA )
      }
      wideFormat <- object$wideFormat
      microData <- ( ! wideFormat ) & is.null( mf$`(freq)` )
      microDataInt <- as.integer(microData)
      if( microData & ( type=="random" ) ) {
         impVarName <- if( is.null(impVarName) ) 
            paste( object$responseVarName, ".imp", sep="" ) else
	    as.character(impVarName)[1]
	 if( impVarName == object$responseVarName ) impVarName <-
	    paste( object$responseVarName, ".imp", sep="" )
      }
      freqSupplied <- wideFormat | ( ! microData )
      if( wideFormat ) {
         fMat <- model.response(mf)
	 if( any(is.na(fMat)) ) stop( gettext(
	    "Missing responses not allowed in wide-format 'newdata'"),
	    domain = NA )
	 if( any( fMat < 0 ) ) stop( gettext(
	    "Negative response frequencies encountered in 'newdata'"),
	    domain = NA )
         freqRowInputData <- apply( model.response(mf), 1, sum )
	 responseRowInputData <- integer(0L)
      } else {
         freqRowInputData <- if( freqSupplied ) mf$`(freq)` else
	    rep(1, NROW(mf))
	 if( any(is.na(freqRowInputData)) ) stop( gettext(
            "Missing values in 'freq' are not allowed"), domain = NA )	    
	 if( any( freqRowInputData < 0 ) ) stop( gettext(
            "Negative values in 'freq' are not allowed"), domain = NA )
         y <- model.response(mf)
	 if( any( is.na(y) ) ) stop(gettext(
	    "Response variable in 'newdata' contains missing values"),
	    domain = NA )
	 if( ! is.factor(y) ) stop( gettext(
	    "Response variable in 'newdata' is not a factor"), domain = NA )
	 if( nlevels(y) != object$dimVec["nLevels"] ) stop( gettext(
	    "Response in 'newdata' has incorrect levels"), domain = NA )
	 if( any( levels(y) != object$responseLevels ) ) stop( gettext(
	    "Response in 'newdata' has incorrect levels"), domain = NA )
         responseRowInputData <- unclass(y)
         fMat <- matrix(0, NROW(mf), object$dimVec["nLevels"] )
	 storage.mode(fMat) <- storage.mode(freqRowInputData)
	 for( j in 1:NCOL(fMat) )  fMat[y==j,j] <- freqRowInputData[y==j]
	 rownames(fMat) <- rownames(mf)
         colnames(fMat) <- object$responseLevels
      }
      # identify unique covariate patterns in newdata
      mvcode <- object$specialCodes["mvcode"]
      x <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)
      x[ is.na(x) ] <- mvcode
      xdf <- as.data.frame(x)
      cx <- do.call( "paste", c( structure( xdf[,,drop=FALSE], names=NULL ),
         sep="\r" ) )
      uniqueRowNumbers.x <- (1:NROW(x))[ ! duplicated(cx) ]
      cxUnique <- cx[ uniqueRowNumbers.x ]
      nCovPatt <- length(cxUnique)
      covPatt <- match(cx, cxUnique)
      names(covPatt) <- NULL
      modelMatrix <- x[uniqueRowNumbers.x,,drop=FALSE]
      # dummy data patterns
      nDataPatt <- nCovPatt
      dataPatt <- covPatt
      covPattForDataPatt <- 1:nDataPatt
      responseForDataPatt <- rep(1L,nDataPatt)
      # aggregate frequencies
      dFtmp <- data.frame( freqRowInputData, covPatt, dataPatt )
      names(dFtmp) <- c( "freq", "covPatt", "dataPatt" )
      aggF <- aggregate( freq ~ covPatt, data=dFtmp, FUN=sum )
      aggF <- aggF[ order(aggF$covPatt), ]
      freqForCovPatt <- aggF$freq
      names(freqForCovPatt) <- NULL
      storage.mode(freqForCovPatt) <- "double"
      freqForDataPatt <- freqForCovPatt
      # handle fittedPi and fittedPistar
      if( object$saturated ) {
         mmOrig <- as.data.frame( object$modelMatrix )
	 mmOrig[ is.na(mmOrig) ] <- mvcode
	 cxOrig <- do.call( "paste",
	    c( structure( mmOrig[,,drop=FALSE], names=NULL ),
            sep="\r" ) )
	 m <- match( cxUnique, cxOrig )
	 fittedPi <- if( meanSeries ) object$fittedPiMCMC[m,] else
	    object$fittedPi[m,]
	 fittedPi[ is.na(fittedPi) ] <- mvcode
	 fittedPistar <- if( meanSeries ) object$fittedPistarMCMC[m,] else
	    object$fittedPistar[m,]
	 fittedPistar[ is.na(fittedPistar) ] <- mvcode
	 rownames(fittedPi) <- rownames(fittedPistar) <- NULL
      } else {
         # these will be recomputed
         fittedPi <- matrix( numeric(), nCovPatt, object$dimVec["nLevels"] )
	 rownames(fittedPi) <- NULL
	 colnames(fittedPi) <- object$responseLevels
	 fittedPi[] <- mvcode
	 fittedPistar <- fittedPi
      }
      # objects to pass to Fortran
      dimVec <- c( nrowInputData = NROW(mf),
         nLevels = as.integer(object$dimVec["nLevels"]),
         nCovPatt = nCovPatt,
	 nDataPatt = nDataPatt,
         ncolMM = NCOL(modelMatrix),
         nParamThisModel = as.integer(object$dimVec["nParamThisModel"]),
	 nParamSatModel = as.integer(object$dimVec["nParamSatModel"]),
         wideFormatInt = as.integer(wideFormat) )
      fitweightRowInputData <- as.double(freqRowInputData)
      fitweightCovPatt <- freqForCovPatt
      fitweightDataPatt <- freqForDataPatt
      if( type == "condMean" ) {
	 storage.mode(freqRowInputData) <- "double"
	 freqRowInputDataInt <- integer(0L)
	 storage.mode(fMat) <- "double"
	 fMatInt <- matrix( integer(), 0L, 0L )
      } else {
         freqRowInputDataInt <- as.integer(freqRowInputData)
	 if( any( freqRowInputDataInt != freqRowInputData ) )
	    warning( gettext(
	    "Some frequencies changed when integerized" ), domain = NA )
         freqRowInputData <- numeric(0L)
	 fMatInt <- fMat
	 storage.mode(fMatInt) <- "integer"
	 if( any( fMatInt != fMat ) )
	    warning( gettext(
	    "Some frequencies changed when integerized" ), domain = NA )
         fMat <- matrix( numeric(), 0L, 0L )	 
      }
      if( type == "condMean" ) {
         condMeans <- matrix( numeric(1L),
            dimVec["nrowInputData"], dimVec["nLevels"] )
         rownames(condMeans) <- rownames(mf)
         colnames(condMeans) <- object$responseLevels
	 impMat <- matrix( integer(), 0L, 0L )
	 impVec <- integer(0L)
      } else {
         condMeans <- matrix( numeric(), 0L, 0L )
	 impMat <- matrix( 0L,
	    dimVec["nrowInputData"], dimVec["nLevels"] )
	 impVec <- integer( dimVec["nrowInputData"] )
	 rownames(impMat) <- names(impVec) <- rownames(mf)
	 colnames(impMat) <- object$responseLevels
      }
   }
   #--------------------------------------------
   # needed for dotCall64::.C64
   SIGNATURE <- c(
      # inputs
      dimVec = "integer",
      modelMatrix = "double",
      fitweightRowInputData = "double",
      fitweightCovPatt = "double",
      fitweightDataPatt = "double",
      freqForDataPatt = "double",
      covPatt = "integer",
      dataPatt = "integer",
      covPattForDataPatt = "integer",
      responseForDataPatt = "integer",
      surveyModeInt = "integer",
      baselineInt = "integer",
      pertMat = "double",
      pertMatInv = "double",
      priorInt = "integer",
      priorFreqTot = "double",
      priorAllocSuppliedInt = "integer",
      priorAlloc = "double",
      saturatedInt = "integer",
      methodInt = "integer",
      ctrlInt = "integer",
      ctrlReal = "double",
      specialCodes = "double",
      typeInt = "integer",
      condMeanInt = "integer",
      microDataInt = "integer",
      freqRowInputData = "double",
      freqRowInputDataInt = "integer",
      fMat = "double",
      fMatInt = "integer",
      responseRowInputData = "integer",
      coefficients = "double",
      # inouts
      fittedPi = "double",
      fittedPistar = "double",
      # outputs
      condMeans = "double",
      impMat = "integer",
      impVec = "integer",
      # messaging
      status = "integer",
      msgLenMax = "integer",
      msgCodes = "integer",
      msgLenActual = "integer" )
   INTENT <- c(
      # inputs
      dimVec = "r",
      modelMatrix = "r",
      fitweightRowInputData = "r",
      fitweightCovPatt = "r",
      fitweightDataPatt = "r",
      freqForDataPatt = "r",
      covPatt = "r",
      dataPatt = "r",
      covPattForDataPatt = "r",
      responseForDataPatt = "r",
      surveyModeInt = "r",
      baselineInt = "r",
      pertMat = "r",
      pertMatInv = "r",
      priorInt = "r",
      priorFreqTot = "r",
      priorAllocSuppliedInt = "r",
      priorAlloc = "r",
      saturatedInt = "r",
      methodInt = "r",
      ctrlInt = "r",
      ctrlReal = "r",
      specialCodes = "r",
      typeInt = "r",
      condMeanInt = "r",
      microDataInt = "r",
      freqRowInputData = "r",
      freqRowInputDataInt = "r",
      fMat = "r",
      fMatInt = "r",
      responseRowInputData = "r",
      coefficients = "r",
      # inouts
      fittedPi = "rw",
      fittedPistar = "rw",
      # outputs
      condMeans = "w",
      impMat = "w",
      impVec = "w",
      # messaging
      status = "w",
      msgLenMax = "r",
      msgCodes = "w",
      msgLenActual = "w" )
   #--------------------------------------------
   # create a matrix for holding message codes
   msgLenMax <- 40L
   msgCodes <- matrix( 0L, msgLenMax, 17L )
   #--------------------------------------------
   result <- dotCall64::.C64("rrlogit_impute",
      SIGNATURE = SIGNATURE,
      # inputs
      dimVec = dimVec,
      modelMatrix = modelMatrix,
      fitweightRowInputData = fitweightRowInputData,
      fitweightCovPatt = fitweightCovPatt,
      fitweightDataPatt = fitweightDataPatt,
      freqForDataPatt = freqForDataPatt,
      covPatt = covPatt,
      dataPatt = dataPatt,
      covPattForDataPatt = covPattForDataPatt,
      responseForDataPatt = responseForDataPatt,
      surveyModeInt = 0L,
      baselineInt = object$baselineInt,
      pertMat = pertMat,
      pertMatInv = pertMatInv,
      priorInt = object$priorInt,
      priorFreqTot = object$priorFreqTot,
      priorAllocSuppliedInt = object$priorAllocSuppliedInt,
      priorAlloc = object$priorAlloc,
      saturatedInt = object$saturatedInt,
      methodInt = object$methodInt,
      ctrlInt = object$ctrlInt,
      ctrlReal = object$ctrlReal,
      specialCodes = object$specialCodes,
      typeInt = typeInt,
      condMeanInt = if( type == "condMean" ) 1L else 0L,
      microDataInt = microDataInt,
      freqRowInputData = freqRowInputData,
      freqRowInputDataInt = freqRowInputDataInt,
      fMat = fMat,
      fMatInt = fMatInt,
      responseRowInputData = responseRowInputData,
      coefficients = coefficients,
      # inouts
      fittedPi = fittedPi,
      fittedPistar = fittedPistar,
      # outputs
      condMeans = condMeans,
      impMat = impMat,
      impVec = impVec,
      # messaging
      status = integer(1L),
      msgLenMax = msgLenMax,
      msgCodes = msgCodes,
      msgLenActual = integer(1L),
      # other args to .C64
      NAOK = TRUE,
#      INTENT = INTENT,
      PACKAGE = "rrLogit")
      #--------------------------------------------
      # display message from Fortran, if present
      msg.lines <- .msg( result$msgCodes, result$msgLenActual )
      if( is.null( msg.lines ) ){
         msg <- "OK"
      } else{
         msg <- paste0( msg.lines, collapse="\n" )
      }
      msg <- paste( msg, "\n", sep="")
      if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #--------------------------------------------
   if( result$status != 0 ) stop( gettext( 
      "Procedure aborted" ), domain = NA )
   #--------------------------------------------
   if( type == "condMean" ) {
      condMeans <- result$condMeans
      condMeans[ condMeans == mvcode ] <- NA
      return(condMeans)
   } else if( type == "random" ) {
      if( microData ) {
         impVec <- result$impVec
	 impVec[ impVec == -1L ] <- NA
         impVec <- factor( impVec, levels = 1:dimVec["nLevels"] )
	 levels( impVec ) <- object$responseLevels
	 impVec <- as.data.frame( impVec )
	 names( impVec ) <- impVarName
	 if( is.null(newdata) ) {
	    rownames(impVec) <- rownames(object$mF)
	    return(impVec)
	 } else {
	    rownames(impVec) <- rownames(newdata)
	    return( cbind( newdata, impVec ) )
	 }
      } else {
         impMat <- result$impMat
	 impMat[ impMat == -1L ] <- NA
	 impMat <- as.data.frame(impMat)
	 names(impMat) <- paste( names(impMat), ".imp", sep="" )
	 if( is.null(newdata) ) {
	    rownames(impMat) <- rownames(object$mF)
	    return(impMat)
	 } else {
	    rownames(impMat) <- rownames(newdata)
	    return( cbind( newdata, impMat ) )
	 }
      }
   }
}
