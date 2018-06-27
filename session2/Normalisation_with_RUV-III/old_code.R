#
old code

```{r}
boxplot(ExpressionData - rowMedians(ExpressionData), 
        ylim = c(-2.5 , 2.5), col = 'black', xaxt = 'n', yaxt = 'n' ,
        ylab = '', xlab = '', cex.lab = 3, cex.main = 2, cex.axis = 2.5 ,
        las = 1, outline = FALSE, frame = FALSE, whisklty = 0, staplelty = 0, names=FALSE,
        main = 'Unnormalized data')
axis(2, mgp = c(3.5,.9,0), lwd.ticks = 3, las = 1, cex.axis = 1.5)
mtext('RLE', 2 , line = 2 , cex = 1.5)
med_RLE <- apply(ExpressionData - rowMedians(ExpressionData), 2, median)
points(c(1:959), med_RLE, bg = ColorBatch[factor(SampleInformation$Batch)], col = 'black', pch = 21, cex = .6)
abline(h = 0, col = 'cyan', lwd = 3, lty = 2)
box(lwd = 5, bty = 'l')
```


boxplot(t(RUVcorrected) - rowMedians(t(RUVcorrected)), 
        ylim = c(-2.5 , 2.5), col = 'black', xaxt = 'n', yaxt = 'n' ,
        ylab = '', xlab = '', cex.lab = 3, cex.main = 2, cex.axis = 2.5 ,
        las = 1, outline = FALSE, frame = FALSE, whisklty = 0, staplelty = 0, names=FALSE,
        main = 'RUV-III normalized')
axis(2, mgp = c(3.5,.9,0), lwd.ticks = 3, las = 1, cex.axis = 1.5)
mtext('RLE', 2 , line = 2 , cex = 1.5)
med_RLE <- apply(t(RUVcorrected) - rowMedians(t(RUVcorrected)), 2, median)
points(c(1:959), med_RLE, bg = ColorBatch[factor(SampleInformation$Batch)], col = 'black', pch = 21, cex = .6)
abline(h = 0, col = 'cyan', lwd = 3, lty = 2)
box(lwd = 5, bty = 'l')