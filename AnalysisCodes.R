comparison.f <- function(te.mat,te.samp,g1,g2,exp=FALSE){    
    te.mat.TPM <- te.mat[,grep("TPM",colnames(te.mat))]
    colnames(te.mat.TPM) <- gsub("_TPM","",colnames(te.mat.TPM))
    rownames(te.mat.TPM) <- te.mat[,"Gene_Symbol"]
    if (exp)    return (te.mat.TPM)
    g1.samp <- te.samp[te.samp[,"Sample.Group"] == g1,1]
    g2.samp <- te.samp[te.samp[,"Sample.Group"] == g2,1]
    te.mat.TPM <- te.mat.TPM[,c(g1.samp,g2.samp)]
    com.re <- apply(te.mat.TPM,1,function(x){
        g1.exp <- as.double(as.matrix(x[g1.samp]))
        g2.exp <- as.double(as.matrix(x[g2.samp]))
        p.v <- t.test(g1.exp,g2.exp)$p.value
        fc <- mean(g2.exp) / mean(g1.exp)
        if (fc < 1)    fc <- -1 * (1/fc)
        cbind(p.v,fc)
    })
    if (nrow(com.re) == 2)    com.re <- t(com.re)
    colnames(com.re) <- c("P.value","FoldChange")
    f.re <- cbind(com.re,te.mat.TPM)
    #rownames(f.re) <- te.mat[,"Gene_Symbol"]
    return (f.re)
}

deseq.f <- function(te.mat,te.samp,g1=NULL,g2=NULL,exp=FALSE){
    te.mat.ct <- te.mat[,grep("_Read_Count",colnames(te.mat))]
    rownames(te.mat.ct) <- te.mat[,"Gene_Symbol"]
    colnames(te.mat.ct) <- gsub("_Read_Count","",colnames(te.mat.ct))
    te.mat.ct <- te.mat.ct[,te.samp[,"Sample.ID"]]
    if (length(g1) & length(g2)){
        te.samp <- te.samp[is.element(te.samp[,"Sample.Group"],c(g1,g2)),]
        te.mat.ct <- te.mat.ct[,te.samp[,"Sample.ID"]]
    }
    dds <- DESeqDataSetFromMatrix(countData=te.mat.ct,colData=te.samp,design=~Sample.Group)
    dds <- DESeq(dds)
    res <- results(dds)
    sig.mat <- as.matrix(res[,c("padj","log2FoldChange")])
    nms <- sig.mat[,"log2FoldChange"] < 0
    sig.mat[,"log2FoldChange"] <- 2^sig.mat[,"log2FoldChange"]
    sig.mat[nms,"log2FoldChange"] <- -1 * (1/sig.mat[nms,"log2FoldChange"])    
    dds.exp <- counts(dds,normalized=TRUE)
    exp.mat <- cbind(sig.mat,dds.exp[rownames(sig.mat),])
    colnames(exp.mat)[1:2] <- c("P.value","FoldChange")
    if (exp){
        exp.mat <- exp.mat[,-c(1:2)]
    }
    if (exp == "result")    return (res)
    return (exp.mat)
}

DEG.f <- function(dds.re,p=0.05,fc=1.5,met=NULL){
    p.v <- as.double(as.matrix(dds.re[,1]))
    if (length(met)){
        p.v <- p.adjust(p.v,met)
    }
    sig.p <- which(p.v < p)
    sig.fc <- which(abs(as.double(as.matrix(dds.re[,2]))) > fc)
    sig <- intersect(sig.p,sig.fc)
    sig.dds <- dds.re[sig,]
    return (sig.dds)
}

venn.f <- function(venn.li,out=NULL){
    f.li <- lapply(venn.li,function(x){
        rownames(x)
    })
    names(f.li) <- names(venn.li)
    venn.re <- venn.diagram(x=f.li,
        category.names = names(f.li),
        filename = out)
    if (!length(out)){
        grid.newpage()
        grid.draw(venn.re)
    }
}

subset.union.f <- function(total.li,total.exp){
    u.genes <- NULL
    for (i in 1:length(total.li)){
        u.genes <- c(u.genes,rownames(total.li[[i]]))
    }
    u.exp <- total.exp[unique(u.genes),]
    return (u.exp[,-c(1:2)])
}


vol.f <- function(te.exp,cut.pv,cut.fc,xl,yl,out=NULL){
    #cut.pv <- -log10(0.05)
    #cut.fc <- 1

    logp <- -log10(te.exp[,1])

    fc <- te.exp[,2]
    neg <- fc < 0
    fc[neg] <- 1/abs(fc[neg])
    fc <- log2(fc)

    thrs <- rep("Stable",length(logp))
    thrs[which(fc > cut.fc & logp > cut.pv)] <- "Up"
    thrs[which(fc < -1 * cut.fc & logp > cut.pv)] <- "Down"
    thrs[which(fc > -1 & fc < 1 * cut.fc & logp > cut.pv)] <- "SignificantNOtDEG"
    DEG <- factor(thrs,levels=c("Up","Down","SignificantNOtDEG","Stable"))

    gplot.df <- data.frame(x=fc,y=logp)
    gp <- ggplot(data=gplot.df, aes(x=x, y=y)) +
        geom_point(size=2,aes(color=DEG)) +
        theme(legend.position = "none") + 
        xlim(xl) + ylim(yl) + 
        scale_color_manual(values=c("#F8766D","#7CAE00","#00BFC4","grey")) + 
        geom_vline(xintercept=c(-1*cut.fc,cut.fc),lty=4,col="black",lwd=0.8) + 
        geom_hline(yintercept=cut.pv,lty=4,col="black",lwd=0.8) +
        xlab("log2 fold change") + ylab("-log10 p-value") + theme_bw()
    if (length(out)){
        png(out,width=2000,height=2000,res=300)
        print (gp)
        dev.off()
    }
    else{
        gp
    }
}

heatmap.f <- function(te.exp,Groups,rown=FALSE,out=NULL){
    te.exp <- te.exp[,is.element(colnames(te.exp),Groups[,1])]
    Groups <- cbind(Groups,Groups=gsub("[0-9]$","",Groups[,2]))
    colnames(Groups)[2] <- "Subgroups"
    rownames(Groups) <- Groups[,1]
    Groups <- Groups[colnames(te.exp),]
    #if (length(genes)){
    #    te.exp[genes,]
    #}
    gp <- pheatmap(te.exp,annotation_col=Groups[,2:3],scale="row",show_rownames=rown)
    if (length(out)){
        png(out,width=2000,height=2000,res=300)
        print (gp)
        dev.off()
    }
    else{
        gp
    }
}

cluster.f <- function(te.exp,Groups,g,met=NULL,out=NULL){
    te.exp <- te.exp[,is.element(colnames(te.exp),Groups[,1])]
    Groups <- cbind(Groups,Groups=gsub("[0-9]$","",Groups[,2]))
    colnames(Groups)[2] <- "Subgroups"
    rownames(Groups) <- Groups[,1]
    Groups <- Groups[colnames(te.exp),]
    #if (length(genes)){
    #    te.exp[genes,]
    #}
    if (met == "PCA"){
        pca.cal.mat <- prcomp(t(te.exp),scale. = TRUE)
        gp <- autoplot(pca.cal.mat, data=Groups, colour = g) + theme_bw()
    }
    if (met == "tsne"){
        tsne.mat <- tsne(t(te.exp))
        colnames(tsne.mat) <- c("X1","X2")
        gp <- ggplot(data=tsne.mat,aes(x=X1,y=X2)) + geom_point(aes(color=Groups[,g])) + theme_bw() +
            scale_colour_discrete(name=g)
    }
    if (met == "MDS"){
        d <- dist(t(te.exp))
        fit <- cmdscale(d, eig=TRUE, k=length(unique(Groups[,g])))
        mds.mat <- fit$points
        colnames(mds.mat) <- paste0("X",1:length(unique(Groups[,g])))
        gp <- ggplot(data=mds.mat,aes(x=X1,y=X2)) + geom_point(aes(color=Groups[,g])) + theme_bw() +
            scale_colour_discrete(name=g)
    }
    if (length(out)){
        png(out,width=2000,height=2000,res=300)
        print (gp)
        dev.off()
    }
    else{
        gp
    }
}
