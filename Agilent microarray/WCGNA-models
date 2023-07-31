#WGCNA
MA.reomoveCTprobe$genes <- MA.reomoveCTprobe$genes %>% 
  mutate(ENTREZID = bitr(MA.reomoveCTprobe$genes$GeneName, fromType = 'SYMBOL',OrgDb=org.Mm.eg.db, toType = 'ENTREZID',drop = F)$ENTREZID)
MA.reomoveCTprobe <- MA.reomoveCTprobe[!is.na(MA.reomoveCTprobe$genes$ENTREZID),]
eset = exprs.MA(MA.reomoveCTprobe) #get the expression from MA
colnames(eset) <- rownames(targets2)
rownames(eset) <- MA.reomoveCTprobe$genes$GeneName 
eset <- removeBatchEffect(eset, as.factor(targets2$Batch))
sampleinfo <- targets2$Target %>% as.data.frame()%>%
  mutate(rownames(targets2)) %>%
  `colnames<-`(c('group','filenames'))
rownames(sampleinfo) <- sampleinfo$filenames
table(sampleinfo$group)
#3个大组做WCGNA
sampleinfo <- sampleinfo %>% mutate(threegroups= case_when(sampleinfo$group %in% c("C7plus6", "C12") ~ "CERmodel",
                                            sampleinfo$group %in% c("FAEE150", "FAEE50") ~ "FAEEmodel",
                                            sampleinfo$group %in% c("TLCS3", "TLCS5") ~ "TLCSmodel",
                                            TRUE ~ sampleinfo$group))
x=c("CERmodel",'FAEEmodel','TLCSmodel')
WGCNAinfo <- sampleinfo[sampleinfo$group%in% x,]
WGCNAinfo <- sampleinfo[sampleinfo$threegroups%in% x,]
eset_WGCNA <- eset[,WGCNAinfo$filenames]
exprMat <- eset_WGCNA
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
type = "unsigned"
# corType: pearson or bicor
corType = "bicor"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- exprMat
dim(dataExpr)
head(dataExpr)[,1:8]
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                                  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
        printFlush(paste("Removing genes:", 
                         paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", 
                         paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
# sample network based on squared Euclidean distance note that we
# transpose the data
A = adjacency(t(dataExpr), type = "distance")
# this calculates the whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized connectivity
Z.k = scale(k)
# Designate samples as outlying if their Z.k value is below the threshold
thresholdZ.k = -5  # often -2.5
# the color vector indicates outlyingness (red)
outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")
# calculate the cluster tree using flahsClust or hclust
sampleTree = hclust(as.dist(1 - A), method = "average")
# Convert traits to a color representation: where red indicates high
# values
# traitColors = data.frame(numbers2colors(datTraits, signed = FALSE))
# dimnames(traitColors)[[2]] = paste(names(datTraits), "C", sep = "")
# datColors = data.frame(outlierC = outlierColor, traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, groupLabels = names(outlierColor), 
                    colors = outlierColor, 
                    main = "Sample dendrogram and trait heatmap")
plotDendroAndColors(sampleTree, groupLabels = names(outlierColor), 
                    colors = outlierColor, 
                    main = "Sample dendrogram and trait heatmap")
dev.off()
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# R-square=0.85
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off()
power = 20#get power from plot
power
if (is.na(power)){
    print("Using experience power since no suitable power found.")
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))       
                   )
    )
}

print(paste("Finally chooosed power is :", power))
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 25,
                       networkType = type,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOM=TRUE,
                       TOMDenom = "min",  deepSplit = 1,
                       stabilityCriterion = "Individual fraction", 
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3, randomSeed=1117)
table(net$colors)
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
pdf(file="plotDendroAndColors.pdf", onefile=F, paper="special", 
    bg="white", pointsize=6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.5,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
dynamicColors <- labels2colors(net$unmergedColors)
plotDendroAndColors(net$dendrograms[[1]], cbind(dynamicColors,moduleColors),
                    c("Dynamic Tree Cut", "Module colors"),
                    dendroLabels = FALSE, hang = 0.5,
                    addGuide = TRUE, guideHang = 0.05)
gene_module <- data.frame(ID=colnames(dataExpr), module=moduleColors)
gene_module = gene_module[order(gene_module$module),]
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
MEs_colt = as.data.frame(t(MEs_col))
colnames(MEs_colt) = rownames(dataExpr)
write.table(MEs_colt,file="module_eipgengene.xls",
            sep="\t",quote=F)
pdf(file="plotEigengeneNetworks.pdf", onefile=F, paper="special", 
    bg="white", pointsize=6)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
group <- factor(WGCNAinfo$group, levels=unique(WGCNAinfo$group))
group <- factor(WGCNAinfo$threegroups, levels=unique(WGCNAinfo$threegroups))
traitData<- model.matrix(~0+group) #make the design with each color of total sample
colnames(traitData) <- unique(WGCNAinfo$group)
colnames(traitData) <- unique(WGCNAinfo$threegroups)
rownames(traitData) <- rownames(WGCNAinfo)
traitData <- traitData[,c("Control","C7plus6","FAEE150","TLCS3")]
MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                     marDendro = c(4,4,4,4),
                      marHeatmap = c(5,5,5,5), plotDendrograms = T, 
                      xLabelsAngle = 90)
hubs = chooseTopHubInEachModule(dataExpr, colorh=moduleColors, power=power, type=type)
hubs
con <- nearestNeighborConnectivity(dataExpr, nNeighbors=50, power=power,
                                   type=type, corFnc = corType)

#TOM plot
TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
TOM <- as.matrix(TOM)
dissTOM = 1-TOMsimilarityFromExpr(dataExpr, power = power)
nSelect = 400 
#设置种子序列，保证结果具有可重复性
set.seed(10)
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
#get 400 genes for TOM plot
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
#heatmap for TOM plot
pdf(file="TOM_plot.pdf")
plotDiss = selectTOM^7
diag(plotDiss) = NA
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotDiss, selectTree,col=myheatcol, selectColors, main = "Network heatmap plot, selected genes")
dev.off()
if (corType=="pearsoon") {
        modTraitCor = cor(MEs_col, traitData, use = "p")
        modTraitP = corPvalueStudent(modTraitCor, nSamples)
    } else {
        modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
        modTraitCor = modTraitCorP$bicor
        modTraitP   = modTraitCorP$p
    }
# signif表示保留几位小数 
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")  
dim(textMatrix) = dim(modTraitCor)
pdf(file="labeledHeatmap.pdf", onefile=F, paper="special", 
        pointsize=6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
                   yLabels = colnames(MEs_col), 
                   cex.lab = 1, 
                   ySymbols = colnames(MEs_col), colorLabels = FALSE, 
                   colors = blueWhiteRed(50), 
                   textMatrix = textMatrix, setStdMargins = T, 
                   cex.text = 0.5, zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
dev.off()
    
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
                   yLabels = colnames(MEs_col), 
                   cex.lab = 0.75, 
                   ySymbols = colnames(MEs_col), colorLabels = FALSE, 
                   colors = blueWhiteRed(50), 
                   textMatrix = textMatrix, setStdMargins = T, 
                   cex.text = 0.5, zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
modTraitCorMelt = as.data.frame(modTraitCor)
write.table(modTraitCorMelt,file=("module_trait_correlation.xls"),
                sep="\t",quote=F)
modTraitCorMelt$ID = rownames(modTraitCor)
modTraitCorMelt = melt(modTraitCorMelt)
colnames(modTraitCorMelt) <- c("Module","Trait","PersonCorrelationValue")
modTraitPMelt = as.data.frame(modTraitP)
write.table(modTraitPMelt,file=("module_trait_correlationPvalue.xls"),
                sep="\t",quote=F)
modTraitPMelt$ID = rownames(modTraitP)
modTraitPMelt = melt(modTraitPMelt)
colnames(modTraitPMelt) <- c("Module","Trait","Pvalue")
modTraitCorP = merge(modTraitCorMelt, modTraitPMelt, by=c("Module","Trait"))
write.table(modTraitCorP,file=("module_trait_correlationPvalueMelt.xls"),
                sep="\t",quote=F,row.names=F)

#GS and MM relationship
if (corType=="pearsoon") {
        geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
        MMPvalue = as.data.frame(corPvalueStudent(
            as.matrix(geneModuleMembership), nSamples))
    } else {
        geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
        geneModuleMembership = geneModuleMembershipA$bicor
        MMPvalue   = geneModuleMembershipA$p
    }
    
    if (corType=="pearsoon") {
        geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
        geneTraitP = as.data.frame(corPvalueStudent(
            as.matrix(geneTraitCor), nSamples))
    } else {
        geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
        geneTraitCor = as.data.frame(geneTraitCorA$bicor)
        geneTraitP   = as.data.frame(geneTraitCorA$p)
    }
geneTraitCorMelt = as.data.frame(geneTraitCor)
write.table(geneTraitCorMelt,file=("gene_trait_correlation.xls"),
                sep="\t",quote=F)
geneTraitCorMelt$ID = rownames(geneTraitCor)
geneTraitCorMelt = melt(geneTraitCorMelt)
colnames(geneTraitCorMelt) <- c("Gene","Trait","PersonCorrelationValue")
geneTraitPMelt = as.data.frame(geneTraitP)
write.table(geneTraitPMelt,file=("gene_trait_correlationPvalue.xls"),
                sep="\t",quote=F)
geneTraitPMelt$ID = rownames(geneTraitP)
geneTraitPMelt = melt(geneTraitPMelt)
colnames(geneTraitPMelt) <- c("Gene","Trait","Pvalue")
#geneTraitCorP = cbind(geneTraitCorMelt, Pvalue=geneTraitPMelt$Pvalue)
geneTraitCorP = merge(geneTraitCorMelt, geneTraitPMelt, by=c("Gene","Trait"))
write.table(geneTraitCorP,
                file=("gene_trait_correlationPvalueMelt.xls"),
                sep="\t",quote=F,row.names=F)
    
plot_me_trat <- cbind(dynamicColors,moduleColors,geneTraitCor)
geneTraitCorColor <- numbers2colors(geneTraitCor)
plotDendroAndColors(net$dendrograms[[1]],
                        cbind(dynamicColors,moduleColors,geneTraitCorColor),
                        c("Dynamic Tree Cut", "Module colors", colnames(geneTraitCor)),
                        dendroLabels = FALSE, hang = 0.5,
                        addGuide = TRUE, guideHang = 0.05)
    
pdf(file = "plotDendroAndColors2.pdf")
plotDendroAndColors(net$dendrograms[[1]],
                        cbind(dynamicColors,moduleColors,geneTraitCorColor),
                        c("Dynamic Tree Cut", "Module colors", colnames(geneTraitCor)),
                        dendroLabels = FALSE, hang = 0.5,
                        addGuide = TRUE, guideHang = 0.05)
dev.off()
#define the model
module = "turquoise" #for FAEE-AP
#module = "brown"#for CER-&TLCS-AP
pheno = "FAEEmodel"
#pheno = "CERmodel"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module,modNames)
pheno_column = match(pheno,colnames(traitData))
moduleGenes = moduleColors == module
rownames()
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
if (corType=="pearsoon") {
        geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
        MMPvalue = as.data.frame(corPvalueStudent(
            as.matrix(geneModuleMembership), nSamples))
    } else {
        geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
        geneModuleMembership = geneModuleMembershipA$bicor
        MMPvalue   = geneModuleMembershipA$p
    }
if (corType=="pearsoon") {
        geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
        geneTraitP = as.data.frame(corPvalueStudent(
            as.matrix(geneTraitCor), nSamples))
    } else {
        geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
        geneTraitCor = as.data.frame(geneTraitCorA$bicor)
        geneTraitP   = as.data.frame(geneTraitCorA$p)
    }
#verboseScatterplot
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),ylim = c(0, 1), xlim = c(0, 1),
                       abs(geneTraitCor[moduleGenes, pheno_column]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", pheno),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
pdf(file = "verboseScatterplot.pdf")
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                       abs(geneTraitCor[moduleGenes, pheno_column]),ylim = c(0, 1), xlim = c(0, 1),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", pheno),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,showPValue = F)
dev.off()
#get the most related genes based from the verboseScatterplot
geneselect1 <- as.data.frame(geneModuleMembership)
#geneselect1 <- geneselect1[geneselect1$MEbrown>0.7,]
geneselect1 <- geneselect1[geneselect1$MEturquoise>0.8,]
geneselect2 <- as.data.frame(geneTraitCor)
#geneselect2 <- geneselect2[geneselect2$CERmodel>0.35&geneselect2$TLCSmodel>0.35,]
geneselect2 <- geneselect2[geneselect2$FAEEmodel>0.85,]
geneselect <- intersect(rownames(geneselect1),rownames(geneselect2))
geneselect
FAEEmodel <- geneselect
#CERmodel<- geneselect
#annotation for genes, and then directly to do the ORA enrichment analyze
gene =bitr(FAEEmodel, fromType = "SYMBOL", toType = c("ENTREZID"), 
        OrgDb = org.Mm.eg.db ) #%>%  # from 'symbol' to 'ENTREZID'
