# WIND devel

Y <- as.matrix(assays(scEx_log)[[1]])
trueclass <- projections$dbCluster

ctStruct = createRef(Y, trueclass)
plot(ctStruct$hc, xlab="", axes=FALSE, ylab="", ann=FALSE)


wNMI(ctStruct, trueclass, projections$sampleNames)

weights = createWeights(Y, trueclass)

ri = wRI(trueclass, projections$sampleNames, 
    weights$W0, weights$W1)[1:3]

barplot(ri, beside=TRUE, ylim=c(0.6,1.05), 
        legend.text=TRUE, xpd=FALSE)
