# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0, DESeq2 1.38.3
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(fgsea)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GSA)
library(DOSE)
library(enrichR)
library(ggplot2)
library(enrichplot)
library(pathview)




# load series and platform data from GEO

gset <- getGEO("GSE30310", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL9392", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pData(gset)
fData(gset)
colnames(gset)
rownames(gset)

attributes(pData(gset))
attributes(fData(gset))
# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11111011111001011011101010111111111111111111111110",
               "11100001100001100110011100110111001001101110101101",
               "11111111111111111011101111111111011100110111011111",
               "1111101101111111")
sml <- strsplit(gsms, split="")[[1]]

exprs(gset)


range(exprs(gset))


# log2 transformation
ex <- exprs(gset)
hist(exprs(gset), main="Histogram of not log-transformed data")
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

hist(exprs(gset), main="Histogram of log-transformed data")

# log2 transformation
ex <- exprs(gset)
ex

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("HIVnegative","HIVpositive"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)


t_test_results <- apply(gset, 1, function(x) t.test(x[gs=="HIVpositive"], x[gs=="HIVnegative"]))
t_test_results
getOption("max.print")
# extract t-statistic and p-value from t-test results
t_statistic <- sapply(t_test_results, function(x) x$statistic)
t_statistic
p_value <- sapply(t_test_results, function(x) x$p.value)
write.csv(ex, file = "expression_data.csv", row.names = TRUE, col.names = TRUE)
df = read.csv("expression_data.csv", header=TRUE)
df$log_fold_change <- apply(gset, 1, function(x) log2(mean(x[gs=="HIVpositive"])/mean(x[gs=="HIVnegative"])))

# Correct the p-values for multiple comparisons using Holm correction
p_values_corrected <- p.adjust(p_value, method = "holm")

# Create a volcano plot to visualize the results
plot(log_fold_change, -log10(p_values_corrected), pch=20, main="Volcano Plot", xlab="Log2 Fold Change", ylab="-log10(P-value)", xlim=c(-6,6), ylim=c(0,10), col=ifelse(abs(log_fold_change) > 1 & p_values_corrected < 0.05, "red", "black"))
abline(h=-log10(0.05), col="blue", lty=2)
abline(v=c(-1,1), col="blue", lty=2)


fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="holm", sort.by="B", p.value=0.05)
View(tT)



# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="holm", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="holm", p.value=0.05)

dT
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))


# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))



################################################################
# General expression data analysis
ex <- exprs(gset)


                   
                   
                   DEGs <- tT$Gene.symbol
              
                   
                   geneList <- bitr(DEGs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
                   
                   # Perform enrichment analysis using Gene Ontology biological process
                   res <- enrichGO(gene     = geneList$ENTREZID,
                                   keyType  = "ENTREZID",
                                   OrgDb    = org.Hs.eg.db,
                                   ont      = "ALL",
                                   pAdjustMethod = "holm",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   minGSSize     = 15,
                                   maxGSSize     = 500)
                   res
                   # Plot the results
                   
                   barplot(res, showCategory = 20, title = "GSEA Results", xlab = "Enrichment Score")
                   
                   dotplot(res, showCategory = 20)
                   
                   
                   enr_Results<-enrichr(gene=DEGs,database="KEGG_2019_Human")
                   enr_Results
                   
                   heatplot(res, showCategory = 30)
                   cnetplot(res, showCategory = 30, layout = "kk")
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   
                   