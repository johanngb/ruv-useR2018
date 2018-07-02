---
title: "Removing unwanted variation - III (RUV-III)"
subtitle: "useR! 2018 - Brisbane"
output:
  html_document:
    df_print: paged
    toc: yes
    toc_depth: 2
    toc_float: yes
  github_document:
    toc: yes
    toc_depth: 3
---


```{r}
library(ruv)
library(matrixStats)
library(ggplot2)
library(knitr)
library(corrplot)
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE,eval=TRUE,cache=TRUE)
library(here)
```

# Load expression data and samples information

```{r}
ExpressionData <- read.csv('ExpressionData_RawData_Main.csv', stringsAsFactors = FALSE, as.is = TRUE, row.names = 1)
dim(ExpressionData)
ExpressionData <- as.matrix(log2(ExpressionData))
```

The gene expression data matrix contains `r dim(ExpressionData)[1]` genes (rows) and `r dim(ExpressionData)[2]` samples (columns). This data has been generated with the [Nanostring technology](https://www.nanostring.com/) which only targets up to 800 genes of interest.

```{r}
SampleInformation <- read.csv('SampleInformation_Main.csv', stringsAsFactors = FALSE, as.is = TRUE)
dim(SampleInformation)
kable(head(SampleInformation))
```

The variable `Tissue` is the biology that we are interested in and we would expect that gene expression separates the different tissues. The variable `Batch` distinguishes samples processed with reagents bought at different times. The component of the regeants itself are the same but they were provided in different batches. 

# Visualise your data

## Where `ruv` meets `ggplot`: `ruv::ruv_svdplot`

```{r}
ColorBatch <- c("#c51b8a", "blue", "darkgreen") 
ruv_svdplot(t(ExpressionData), k = c(1,2)) + list(aes(color=SampleInformation$Batch), 
                                                  scale_color_manual(values=ColorBatch), labs(color="Batches") )
```

```{r}
ruv_svdplot(t(ExpressionData), k = c(1,2)) + list(aes(color=SampleInformation$Tissue), labs(color="Tissues"))

# Also works
ruv_svdplot(t(ExpressionData), k = c(1,2),info = SampleInformation) + geom_point(aes(color=Tissue)) + labs(color="Tissues")

```

```{r}
ruv_svdplot(t(ExpressionData), k = c(1,3)) + 
  list(aes(color=SampleInformation$Batch), 
       scale_color_manual(values=ColorBatch), labs(color="Batches") )
```

```{r}
ruv_svdplot(t(ExpressionData), k = c(1,3)) + list(aes(color=SampleInformation$Tissue), labs(color="Tissues"))
```

## Relative Log Expression (RLE) plots

**Add sentence about RLE** Ramyar.

It is important to notice that a *bad looking* RLE plot means that there is unwanted variantion in your data but a *good looking* RLE plot does not necesserily mean that your data is correctly normalised. For more information about RLE plot see reference paper [
RLE Plots: Visualising Unwanted Variation in High Dimensional Data](https://arxiv.org/abs/1704.03590), Luke C. Gandolfo, Terence P. Speed 2017.

In the plot below, medians have been highlited for better visualisation.

```{r}
ColorBatch <- c("#c51b8a", "blue", "darkgreen") 
ruv_rle(Y = t(ExpressionData),ylim = c(-3 , 3)) + geom_point(aes(x = rle.x.factor,y= middle,colour=SampleInformation$Batch)) + theme(legend.position = "bottom") + labs(colour="Batches") +
       scale_color_manual(values=ColorBatch) + geom_hline(yintercept = 0,linetype="dotted",colour="cyan")+ggtitle('Unnormalized data')
```


## Any questions?

# Normalise your data 

## Build the matrix of replicates 

RUV-III performs normalisation by exploiting technical replicates which are replicates of the same biological sample. The function `ruv::replicate.matrix` will create a matrix with as many rows as the initial `SampleInformation` data frame and as many columns as the unique biological samples in the study. The order of the rows of `ReplicateMatrix` should be exactly the same as the order of the samples names in the `ExpressionData` matrix. The entries of the `replicate.matrix` are either 1 or 0 to indicate if two samples are or aren't replicates of one another. 

```{r}
ReplicateMatrix <- replicate.matrix(SampleInformation$SamleIds)
dim(ReplicateMatrix)
dim(SampleInformation)
```

The number of technical replicates in the study can be found as the difference between the number of rows and columns of `ReplicateMatrix`. In this study there are 46 techical replicates.


The plot below is one way to understand how the matrix of replicates is structured. In the example below (a subset of `ReplicateMatrix`) one can see that rows 1 and 5 are replicates of the same sample `Sample_1` and therefore they have entry equal to 1. 

```{r}
corrplot(ReplicateMatrix[1:10,1:10],is.corr = FALSE)
```

Every row has always one entry equal to 1. 

```{r}
barplot(rowSums(ReplicateMatrix))
```

Every columns has at least one entry equal to 1. Samples with technical replicates will have more than one entry equal to 1. 

```{r}
barplot(colSums(ReplicateMatrix))
```

## RUV-III 

Here, RUV-III is run with the maximum possible value of `k = 46` which is equal to the number of technical replicates within the study.

```{r}
RUVcorrected <- RUVIII(Y = t(ExpressionData), M = ReplicateMatrix, ctl = c(1:nrow(ExpressionData)), k = 46)
```

Now we can visualise the data after correction and understand if we improved the data.

```{r}
ColorBatch <- c("magenta", "blue", "darkgreen") 
ruv_svdplot(RUVcorrected, k = c(1,2)) + list(aes(color=SampleInformation$Batch), 
                                                  scale_color_manual(values=ColorBatch), labs(color="Batches") )
```

```{r}
ruv_svdplot(RUVcorrected, k = c(1,2)) + list(aes(color=SampleInformation$Tissue), labs(color="Tissues"))
```

```{r}
ColorBatch <- c("#c51b8a", "blue", "darkgreen") 
ruv_rle(Y = RUVcorrected,ylim = c(-3 , 3)) + geom_point(aes(x = rle.x.factor,y= middle,colour=SampleInformation$Batch)) + theme(legend.position = "bottom") + labs(colour="Batches") +
       scale_color_manual(values=ColorBatch) + geom_hline(yintercept = 0,linetype="dotted",colour="cyan")+ggtitle('Normalized data with RUV-III')
```

# RUV-III exercise

Now try to run `RUV-III` on the unnormalised data using different values for `k`. Try `k=1`, `k=10` and `k=20` and comment on the result by visualising the normalised data. Which one do you prefer? And why?


# Session Info

```{r}
sessionInfo()
```

