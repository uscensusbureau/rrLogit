# tests involving a response with no noise added, using alligator data

library(rrLogit)

# example in wide format
fitA <- rrLogit( cbind(Fish, Inv, Rept, Bird, Other) ~ Lake + Sex + Size,
   data=alligatorWide)
# same model, narrow format
fitB <- rrLogit( Food ~ Lake + Sex + Size,
   data=alligatorNarrow, freq=Freq)
# same model, microdata
fitC <- rrLogit( Food ~ Lake + Sex + Size,
   data=alligatorMicro)

# verify that the resulting coefficients are the same
coef(fitA)
all.equal(coef(fitA), coef(fitB))
all.equal(coef(fitA), coef(fitC))

# verify that the resulting fitted values are the same
head(fitted.values(fitA))
all.equal( fitted.values(fitA), fitted.values(fitB) )
all.equal( fitted.values(fitA), fitted.values(fitC) )

# verify that the resulting residuals are the same
head(residuals(fitA))
all.equal( residuals(fitA), residuals(fitB) )
all.equal( residuals(fitA), residuals(fitC) )

# display predictions
head(predict(fitA))
head(predict(fitB))
head(predict(fitC))
all.equal( predict(fitA)[1,], predict(fitB)[1,] )
all.equal( predict(fitA)[1,], predict(fitC)[1,] )

# fit a smaller model without Sex and compare
fitD <- rrLogit( cbind(Fish, Inv, Rept, Bird, Other) ~ Lake + Size,
   data=alligatorWide)
anova(fitD, fitA)

# apply DAP prior with default settings
fitA <- rrLogit( cbind(Fish, Inv, Rept, Bird, Other) ~ Lake + Sex + Size,
   data=alligatorWide, prior="DAP")
fitB <- rrLogit( Food ~ Lake + Sex + Size,
   data=alligatorNarrow, freq=Freq, prior="DAP")
fitC <- rrLogit( Food ~ Lake + Sex + Size,
   data=alligatorMicro, prior="DAP")
coef(fitA)
all.equal(coef(fitA), coef(fitB))
all.equal(coef(fitA), coef(fitC))

# apply DAP prior with custom settings
fitA <- rrLogit( cbind(Fish, Inv, Rept, Bird, Other) ~ Lake + Sex + Size,
   data=alligatorWide,
   prior="DAP", priorFreqTot=3.5, priorAlloc=rep(1/5,5) )
fitB <- rrLogit( Food ~ Lake + Sex + Size,
   data=alligatorNarrow, freq=Freq,
   prior="DAP", priorFreqTot=3.5, priorAlloc=rep(1/5,5) )
fitC <- rrLogit( Food ~ Lake + Sex + Size,
   data=alligatorMicro,
   prior="DAP", priorFreqTot=3.5, priorAlloc=rep(1/5,5) )
coef(fitA)
all.equal(coef(fitA), coef(fitB))
all.equal(coef(fitA), coef(fitC))

# fit intercept-only models
fitA <- rrLogit( cbind(Fish, Inv, Rept, Bird, Other) ~ 1, data=alligatorWide)
fitB <- rrLogit( Food ~ 1, data=alligatorNarrow, freq=Freq)
fitC <- rrLogit( Food ~ 1, data=alligatorMicro)
coef(fitA)
all.equal(coef(fitA), coef(fitB))
all.equal(coef(fitA), coef(fitC))

# fit intercept-only models
fitA <- rrLogit( cbind(Fish, Inv, Rept, Bird, Other) ~ 1, data=alligatorWide)
fitB <- rrLogit( Food ~ 1, data=alligatorNarrow, freq=Freq)
fitC <- rrLogit( Food ~ 1, data=alligatorMicro)
coef(fitA)
all.equal(coef(fitA), coef(fitB))
all.equal(coef(fitA), coef(fitC))
margProportions <- table(alligatorMicro$Food) / NROW(alligatorMicro)
all.equal(as.vector(margProportions), as.vector(fitted(fitA)))
all.equal(as.vector(margProportions), as.vector(fitted(fitB)))
all.equal(as.vector(margProportions), as.vector(fitted(fitC)))


# change the baseline
fitD <- rrLogit( cbind(Fish, Inv, Rept, Bird, Other) ~ 1, data=alligatorWide,
   baseline="Rept")
all(coef(fitA) == coef(fitD))      # should be FALSE
all.equal(fitted(fitA), fitted(fitD))  # should be TRUE

# a few model summaries using print.summary.rrLogit
fitA <- rrLogit( cbind(Fish, Inv, Rept, Bird, Other) ~ Lake + Sex + Size,
   data=alligatorWide)
summary(fitA)
fitB <- rrLogit( Food ~ Lake + Sex + Size,
   data=alligatorNarrow, freq=Freq)
summary(fitB)
fitC <- rrLogit( Food ~ Lake + Sex + Size,
   data=alligatorMicro)
summary(fitC)
