###Read in and normalise data
target <- readTargets("finalFAEETLCS.txt" ) 
target <- readTargets("finaltargetCER.txt" ) 
target <- readTargets("badarray.txt" ) # the faee and tlcs model
target <- readTargets("CERdose.txt" ) #
rownames(target) <- removeExt(target$FileName)
target
RG <- read.maimages(target, source ="agilent") #read files, as two-color array

###Quality Assessment
summary(RG$R)
plotMD(RG)
boxplot(data.frame(log2(RG$Gb)),main="Green background") #quality characteristics of each array
boxplot(data.frame(log2(RG$Rb)),main="Red background")
imageplot(log2(RG$Gb[,1]),RG$printer)
table(RG$genes$ControlType)
imageplot(log2(RG$Rb[,1]), RG$printer, low="white", high="red")
positive <- RG$genes$ControlType %in% '1'
RG <- RG[!positive,] # All positive control probes were filtered before background correction and normalization.

### Perform background correction on the fluorescent intensities
RG.bgcorrect <- backgroundCorrect(RG, method = 'normexp', offset = 25) #offset may depend on the number of sample (25%-50%)

### Normalize the data in array with the 'loess' method
RG.bgcorrect.norm <- normalizeWithinArrays(RG.bgcorrect, method = 'loess')
plotMD(RG.bgcorrect.norm)

### Normalize the data between arrays 
MA <- normalizeBetweenArrays(RG.bgcorrect.norm, method = "Aquantile")
plotDensities(RG.bgcorrect.norm)#check the signal density map
plotDensities(MA)
MA$genes$ControlType %>% unique() #include what category of probe

### remove the control probe
negative <- MA$genes$ControlType %in% '-1'
a <- MA[negative,]
quantile(a$A)
keepProbe <- MA$genes$ControlType == 0
MA.reomoveCTprobe <- MA[keepProbe,] #remove the control probe
MA.reomoveCTprobe$genes$ProbeName %>% table() %>% sort(decreasing = TRUE) %>% head(n = 3) 
MA.reomoveCTprobe.average <- avereps(MA.reomoveCTprobe, ID = MA.reomoveCTprobe$genes$ProbeName) # take an average of probe
IsExpr <- rowSums(MA.reomoveCTprobe$A > 0) >= 5.913736;table(IsExpr) # there is no genes less than 75% Negetive genes

### For replicate gene in each sample, replace values with max one
Amean <- rowMeans(MA.reomoveCTprobe$A)
o <- order(Amean, decreasing=TRUE)
MA.reomoveCTprobe <- MA.reomoveCTprobe[o,]
dup <- duplicated(MA.reomoveCTprobe$genes$GeneName)
MA.reomoveCTprobe <- MA.reomoveCTprobe[!dup,] 

### annotation gene names and filter no significance probe
Annotationname = MA.reomoveCTprobe$genes %>% 
  filter( !grepl( "Unknown", Description ) ) %>% #there are many genes with Description 'unknown' 
  filter( !grepl( "^chr", GeneName ) ) # there are many gene names begin with 'chr'
rownames(Annotationname) <- Annotationname$GeneName

### Array Quality Weights
arrayw <- arrayWeights(MA.reomoveCTprobe)
barplot(arrayw, xlab="Array", ylab="Weight", col="white", las=2)
abline(h=1, lwd=1, lty=2)
w <- matvec(MA.reomoveCTprobe$weights,aw)

### analyze from MA list #in this way, i only can do it up to DEG via topTable function 
targets2 <- targetsA2C(target)

targets2 <- targets2 %>% 
  mutate(modelgroup=ifelse(grepl('Control',targets2$Target),'Control', ifelse(grepl('SalineFAEE',targets2$Target), 'SalineFAEE',
                                                                              ifelse(grepl('SalineTLCS',targets2$Target),'SalineTLCS',
                                                                                     ifelse(grepl('SC',targets2$Target),'SalineCER',
                                                                                            ifelse(grepl('Control',targets2$Target),'Control', 
                                                                                                   ifelse(grepl('TLCS',targets2$Target),'TLCS',
                                                                                                          ifelse(grepl('POAFAEE',targets2$Target),'POAFAEE','CER'))))))))%>% 
  mutate(knockout=ifelse(grepl('ppif',targets2$Target), 'ppif', 'wt'))%>% 
  mutate(rehybridised=ifelse(grepl('8',targets2$Batch), 'badarray', 'normal'))%>%
  mutate(rowname=rownames(targets2))%>%
  mutate(SCcombin = case_when(str_detect(targets2$Target,'wt.SC4plus9')~ 'SC',
                              str_detect(targets2$Target,'wt.SC7plus6')~ 'SC',
                              str_detect(targets2$Target,'wt.C7plus6')~ 'C',
                              str_detect(targets2$Target,'wt.C4plus9')~ 'C',
                              TRUE ~ targets2$Target))

u <- unique(targets2$Target)
f <- factor(targets2$Target, levels=u)
design <- model.matrix(~0+f) #make the design with each color of total sample
colnames(design) <- u
rownames(design) <- rownames(targets2)
corfit <- intraspotCorrelation(MA.reomoveCTprobe, design)
fit <- lmscFit(MA.reomoveCTprobe, design,correlation=corfit$consensus)

cont.matrix <- makeContrasts("wt.SC7-wt.C7","wt.SC7-wt.Control2","wt.C7-wt.Control2", levels=design)

### for 3 model
cont.matrix <- makeContrasts("wtTLCS-wtControl","wtSalineTLCS-wtControl", 'wtTLCS-wtSalineTLCS',levels=design)
cont.matrix <- makeContrasts("wtPOAFAEE-wtControl","wtSalineFAEE-wtControl", 'wtPOAFAEE-wtSalineFAEE',levels=design)
cont.matrix <- makeContrasts("wt.C7rep-wt.Controlrep","wt.SC7rep-wt.Controlrep", 'wt.C7rep-wt.SC7rep',levels=design)
### for different dose of CER
cont.matrix <- makeContrasts('wt.C7-wt.Control', 'wt.SC7-wt.Control','wt.C7-wt.SC7',levels=design)
cont.matrix <- makeContrasts('wt.C7plus6-wt.Control', 'wt.SC7plus6-wt.Control','wt.C7plus6-wt.SC7plus6',levels=design)
cont.matrix <- makeContrasts('wt.C4plus9-wt.Control', 'wt.SC4plus9-wt.Control','wt.C4plus9-wt.SC4plus9',levels=design)
cont.matrix <- makeContrasts('wt.SC7plus6-wt.SC7', 'wt.SC4plus9-wt.SC7','wt.SC7plus6-wt.SC4plus9',levels=design)
cont.matrix <- makeContrasts('wt.C7plus6-wt.C7', 'wt.C4plus9-wt.C7','wt.C7plus6-wt.C4plus9',levels=design)
cont.matrix <- makeContrasts('wt.C7rep-wt.SC7rep',levels=design)
### for bad array (FAEE+TLCS)
cont.matrix <- makeContrasts('wtSalineInj-wtControl','wtEtOH-wtSalineInj','wtEtOH-wtControl',levels=design)
cont.matrix <- makeContrasts('wtFAEE150-wtControl', 'wtFAEE150-wtSalineInj','wtFAEE150-wtEtOH',levels=design)
cont.matrix <- makeContrasts('wtFAEE50-wtControl','wtFAEE50-wtSalineInj','wtFAEE50-wtEtOH',levels=design)
cont.matrix <- makeContrasts('wtFAEE150-wtFAEE50',levels=design)
cont.matrix <- makeContrasts('wtTLCS5-wtControl', 'wtTLCS3-wtControl','wtSalinePerf-wtControl','wtTLCS5-wtSalinePerf', 'wtTLCS3-wtSalinePerf',levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
vennDiagram(results)
topTable(fit2, coef=1, adjust="BH")

degs_temp <- topTable(fit2,coef=1, adjust="BH", n = Inf)%>%na.omit
row.names(degs_temp)=degs_temp$GeneName
logFC_cut = with(degs_temp, mean(abs(logFC))+2*sd(abs(logFC))) # take 95% CI for logFC
degs_temp = degs_temp %>%
  mutate( DEG = factor( ifelse( abs(degs_temp$logFC) > logFC_cut & degs_temp$adj.P.Val < 0.05,
                                ifelse(degs_temp$logFC > 0, 'up','down'),
                                'ns' ), 
                        levels = c("up", "ns", "down" ), ordered = F ) )%>%
  mutate( row = rownames(degs_temp))



###heatmap for DEG validation
chs_x = degs_temp %>%
  slice_max( abs(logFC), n = 5 ) 
#slice_min( P.value, n = 500 )
x <- c("wt.C7plus6","wt.SC7plus6")
heat_matrix = eset[ chs_x$GeneName,targets2$Target%in%x ]

legend_col = data.frame( row.names = sampleinfo$filenames,
                         Group = sampleinfo$group )  # set up the group
bk = 2 # set up the color
brk <- c(seq(-bk,-0.01,by=0.01),seq(0,bk,by=0.01))
heat1 = pheatmap(heat_matrix,
                 scale = "row",
                 annotation_col = legend_col,
                 color = c(colorRampPalette(colors = c("dodgerblue4","white"))(length(brk)/2),
                           colorRampPalette(colors = c("white","brown"))(length(brk)/2)),
                 legend_breaks=seq(-bk,bk,1),
                 breaks=brk,
                 treeheight_row = 40,
                 treeheight_col = 40,
                 border_color = NA,
                 cluster_cols = F,
                 cluster_rows = F,
                 display_numbers = round(heat_matrix,2), number_format = "%.2f", number_color = "grey30",
                 show_rownames = F,
                 show_colnames = F
)
p6 = as.ggplot(heat1)+
  ggtitle('Top 5 DEGs of wt.C7plus6 and wt.SC7plus6')+  
  xlab('Sample') + ylab('Gene')+
  theme(
    plot.title = element_text(hjust = 0.4),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
p6

### annotation for DEG
gpl_anno=read.table('GPL10787.txt', sep = "\t", quote = "\"",header = TRUE, fill = TRUE)  # read the annotation file for the GPL10787 platform (downloaded)
head(gpl_anno)
gpl_anno = gpl_anno %>% distinct(GENE_SYMBOL, .keep_all = TRUE ) # delete the repeat probe
rownames(gpl_anno) = gpl_anno$GENE_SYMBOL
degs_temp = degs_temp[rownames(degs_temp) %in% gpl_anno$GENE_SYMBOL,] #  delete the unmatched in eset
dim(degs_temp)
gpl_anno = gpl_anno[gpl_anno$GENE_SYMBOL %in% rownames(degs_temp),] #  delete the unmatched in annotation
identical( gpl_anno$GENE_SYMBOL, rownames(degs_temp)) #check the arrange of eset and annotation
gpl_anno = gpl_anno[rownames(degs_temp),] #  arrange 
identical( gpl_anno$GENE_SYMBOL, rownames(degs_temp)) #recheck
degs_temp = cbind(gpl_anno,degs_temp)

### volcano plot for DEG
degs_temp$label <- c(rownames(degs_temp)[1:10],rep(NA,(nrow(degs_temp)-10))) 
p1 <- ggplot(degs_temp,aes(logFC, -log10(adj.P.Val)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  geom_point(aes(size=-log10(adj.P.Val), color= -log10(adj.P.Val)))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+    # color change
  scale_size_continuous(range = c(1,4))+    # plot change
  #scale_y_continuous(limits = c(0, 30))+
  ggtitle("wt.Control-ppif.Control") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.5),
        legend.justification = c(0,1)    # adject the location of theme and label
  )+
  guides(col = guide_colourbar(title = "-Log10_q-value"),
         size = "none")+    # label not show
  geom_text(aes(label=label, color = -log10(adj.P.Val)), size = 3, vjust = 1.5, hjust=1)+    # add label
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")    # anno of axis
p1
save.image(file = "DEG.RData")












