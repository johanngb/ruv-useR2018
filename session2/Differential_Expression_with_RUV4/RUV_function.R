


################################################################################################################################
############################################# Moemeneh Foroutan -- Dec 2016 ####################################################
################################################################################################################################

##----------------------- The RUVs function takes these inputs:

## matrix: gene expression matrix that has genes in cols and samples in rows:
## info: a dta frame of annotations, where the biological factor of interest is in a column
## biological_group: the column name for the biological factor of interest in the info data frame
## hk: housekeeping genes with the same ID format as col names in expression matrix (both Entrez IDs for example)
## method = one of the RUV4, RUVinv and RUVrinv
## out_path: a path to a location where the output files and figures are saved


##----------------------- The RUVs function does:

## RUV4/inv/rinv to obtain DEGs in a data with unwanted variations 
## uses housekeeping genes in the first step to obtain DEGs, and 
## then uses the least significant genes as empirical control genes, 
## then re-do the normalisation this time using these emprical DEGs,
## It produces statistics for DE analysis along with visualisation plots


##----------------------- The RUVs function generates these outputs:

## A txt file containing al the genes with their p-values, adj-pvalues, logFC and t statistics
## A txt file containing al only DEGs  with their p-values, adj-pvalues, logFC and t statistics

## png files for p-value histograms obtained using housekepping genes as negative control
## png files for p-value histograms obtained using emprical negative control genes
## png file for three types of p-values generated in RUV methods
## png files for volcano plots when using housekeeping genes as negative control genes
## png files for volcano plots when using emprical negative control genes

################################################################################################################################


RUVs<- function(matrix=mat, info=info, biological_group="treatment", hk=hk, method, out_path){
  
  library(ruv)
  
  ctrl<- colnames(mat) %in% hk
  
  group<- factor(info[,biological_group])  ## take the biological factor of interest
  g<- matrix(group, ncol=1)
  g<- as.numeric(as.factor(g))  ## biological classes(e.g. 1=control   2=TGFb)
  g<- matrix(g, ncol=1)
  
  ##======================================================= RUV4
  
  if (method=="RUV4"){
    estimateK<- getK(mat, g, ctrl, Z = 1, eta = NULL, fullW0 = NULL, cutoff = NULL,
                     method="select", l=1, inputcheck = TRUE)  
    
    k<- estimateK$k
    
    ruv<- RUV4(mat, g, ctrl, k=k, Z = 1, eta = NULL, fullW0 = NULL,
               inputcheck = TRUE)
    
    ## extract p values
    pvals<- t(ruv$p)
    pvals<- data.frame(pvals)
    
    ## adjust p values
    ruv.adj<- variance_adjust(ruv)
    
    ## extract adjusted p values
    adjpvals<- t(ruv.adj$p.BH)
    adjpvals<- data.frame(adjpvals)
    head(adjpvals)
    adjpvals$genes<- row.names(pvals)
    
    ## we want to take the least significant genes as emprical negative control genes
    ## order data based on adjusted p-values
    adjpvals.ordered <- data.frame(adjpvals[order(adjpvals$adjpvals),])
    
    ## take 30% of genes with highest adjusted p-values
    nGenes<- 0.3 * length(ruv.adj$p.BH)  
    empCtrl<- tail(adjpvals.ordered$genes, nGenes)
    
    ## re-run RUV using empical control genes
    empCtrl<- colnames(mat) %in% empCtrl
    
    ruvEmp<- RUV4(mat, g, empCtrl, k=k, Z = 1, eta = NULL, fullW0 = NULL,
                  inputcheck = TRUE)
  }
  
  ##====================================================== RUVinv
  
  if(method== "RUVinv"){
    
    ruv<- RUVinv(mat, g, ctrl, Z = 1, fullW0 = NULL, lambda = NULL, iterN = 100000)
    pvals<- t(ruv$p)
    pvals<- data.frame(pvals)
    
    ruv.adj <- variance_adjust(ruv, ebayes = TRUE, evar = TRUE, rsvar = TRUE, bin = 10, rescaleconst = NULL)
    
    adjpvals<- t(ruv.adj$p.BH)
    adjpvals<- data.frame(adjpvals)
    adjpvals$genes<- row.names(pvals)
    
    adjpvals.ordered <- data.frame(adjpvals[order(adjpvals$adjpvals),])
    
    nGenes<- 0.3 * length(ruv.adj$p.BH) 
    empCtrl<- tail(adjpvals.ordered$genes, nGenes)
    empCtrl<- colnames(mat) %in% empCtrl
    
    ruvEmp<- RUVinv(mat, g, empCtrl, Z = 1, fullW0 = NULL, lambda = NULL, iterN = 100000)
    
  }
  
  ##======================================================= RUVrinv
  
  if(method=="RUVrinv"){
    
    ruv<-RUVrinv(mat, g, ctrl, Z = 1, fullW0 = NULL, lambda = NULL, k=NULL, iterN = 100000)
    pvals<- t(ruv$p)
    pvals<- data.frame(pvals)
    
    ruv.adj <- variance_adjust(ruv, ebayes = TRUE, evar = TRUE, rsvar = TRUE, bin = 10, rescaleconst = NULL)
    
    adjpvals<- t(ruv.adj$p.BH)
    adjpvals<- data.frame(adjpvals)
    head(adjpvals)
    adjpvals$genes<- row.names(pvals)
    
    adjpvals.ordered <- data.frame(adjpvals[order(adjpvals$adjpvals),])
    
    nGenes<- 0.3 * length(ruv.adj$p.BH)  
    empCtrl <- tail( adjpvals.ordered$genes, nGenes)
    empCtrl <- colnames(mat) %in% empCtrl
    
    ruvEmp<- RUVrinv(mat, g, empCtrl, Z = 1, fullW0 = NULL, lambda = NULL, iterN = 100000)
  }
  
  ##---- adjust p values:
 
  ruv.adjEmp <- variance_adjust ( ruvEmp, ebayes = TRUE, evar = TRUE, rsvar = TRUE, bin = 10, rescaleconst = NULL)
  
  ruv.adjpvalsEmp <- t(ruv.adjEmp$p.BH)
  ruv.adjpvalsEmp<- data.frame(ruv.adjpvalsEmp)
  head(ruv.adjpvalsEmp)
  ruv.adjpvalsEmp$genes<- row.names(pvals)
  
  ##---- add betahats (logFC) to this:
  ruv.adjpvalsEmp$logFC<- t(ruv.adjEmp$betahat)
  
  ##---- add t statistics (t) to this:
  ruv.adjpvalsEmp$t<- t(ruv.adjEmp$t)
  
  ## oreder them based on the logFC:
  ruv.adjpvalsEmp.ordered <- data.frame(ruv.adjpvalsEmp[order(ruv.adjpvalsEmp$logFC, decreasing=T),])
  
  ## export all the p vals and logFCs:
     write.table(ruv.adjpvalsEmp.ordered, paste0(out_path, method, "_all_adjPvals_logFC_EmpricalCtrl.txt"), sep="\t", row.names=F)
  
  ##---------------- obtain DEGs:
  DEGs_ruvEmp <- ruv.adjpvalsEmp[ruv.adjpvalsEmp$ruv.adjpvalsEmp < 0.05 &
                                  abs(ruv.adjpvalsEmp$logFC) >  1 ,]  ##  
  
  DEGs_ruvEmp <- data.frame(DEGs_ruvEmp[order(DEGs_ruvEmp$logFC, decreasing=T),])
  
  ## export the p vals, logFCs and t of only DEGs:
     write.table(DEGs_ruvEmp, paste0(out_path, method, "_DEGs_adjPvals_logFC_EmpricalCtrl.txt"), sep="\t", row.names=F)
  
  
  ##-------------- distribution of p values using HK genes:
  
  png(paste0(out_path, "hist_pvals_", method, "_HK_ylim.png"), height=600, width=800)
  par(mar=c(6,6,4,1))
  hist(ruv$p, breaks=100, col="light blue", main= method, ylim=c(1,500), 
       xlab="p-values", cex.lab=1.5, cex.axis=1.5, cex.main=2)
  dev.off()
  
  ##-------------- distribution of p values using emprical control genes:
  
  png(paste0(out_path, "hist_pvals_", method, "_EmpricalCtrl_ylim.png"), height=600, width=800)
  par(mar=c(6,6,4,1))
    hist(ruvEmp$p, breaks=100, col="light blue", main= method, ylim=c(1,500), 
       xlab="p-values", cex.lab=1.5, cex.axis=1.5, cex.main=2)
  dev.off()
  
  png(paste0(out_path, "hist_pvals_", method, "_EmpricalCtrl.png"), height=600, width=800)
  par(mar=c(6,6,4,1))
    hist(ruvEmp$p, breaks=100, col="light blue", main= method,
       xlab="p-values", cex.lab=1.5, cex.axis=1.5, cex.main=2)
  dev.off()

  ##--------------- distributions of the three other p values:
  
  png(paste0(out_path, "hist_threePvals_", method, "_EmprCtrl.png"), height=600, width=2000)
    par(mfrow=c(1,3),mar=c(7, 3, 3, 1))
    hist( ruv.adjEmp $p.evar[1,], breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
        main="p.evar", xlab="",  ylim=c(0,5), cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
    hist( ruv.adjEmp $p.ebayes[1,], breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
        main="p.ebayes",  xlab="", ylim=c(0,5), cex.main=2, cex.lab=1.5, cex.axis=1.5)
  
    hist( ruv.adjEmp $p.rsvar[1,], breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
        main="p.rsvar",  xlab="", ylim=c(0,5), cex.main=2, cex.lab=1.5, cex.axis=1.5)
  dev.off()
  
  ##--------------- volcano plot colored by HK genes:
  
  png(paste0(out_path, "volcano_", method, "_HK.png"), height=600, width=800)
  par(mar=c(6,6,1,1))
    plot(ruv.adj$betahat[1,],-log10(ruv.adj$p.BH[1,]),main="" , col="blue", 
         xlab="logFC", ylab="-log10(adjusted p-values)", cex.axis=1.5, cex.lab=2)
    
    abline(h=1.3, v=c(-1,1),col='red')  ## 1.3 for p-value < 0.05
    points(ruv.adj$betahat[1,ctrl],-log10(ruv.adj$p.BH [1,ctrl]),pch=16, col="green")
  dev.off()
  
  ##--------------- volcano plot colored by emprical control genes:
  
  png(paste0(out_path, "volcano_", method, "_Emprical_controls.png"), height=600, width=800)
  par(mar=c(6,6,1,1))
    plot(ruv.adjEmp$betahat[1,],-log10(ruv.adjEmp$p.BH[1,]),main="" ,col="blue",
         xlab="logFC", ylab="-log10(adjusted p-values)", cex.axis=1.5, cex.lab=2)
    
    abline(h=1.3, v=c(-1,1),col='red')  ## 1.3 for p-value < 0.05
    points(ruv.adjEmp$betahat[1,empCtrl],-log10(ruv.adjEmp$p.BH [1,empCtrl]),pch=16,col="green")
  dev.off()
  
}