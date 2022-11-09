###Read in and normalise data
target <- readTargets("CER.txt" ) #import all the CER group and control group, including 4, 7, 12 injections. I just change the contrast matrix for the different comparation.
rownames(target) <- removeExt(target$FileName)
target
RG <- read.maimages(target, source ="agilent") #read files, as two-color array

###Quality Assessment
summary(RG$R)
plotMD(RG)
boxplot(data.frame(log2(RG$Gb)),main="Green background") #quality characteristics of each array
boxplot(data.frame(log2(RG$Rb)),main="Red background")
imageplot(log2(RG$Gb[,1]),RG$printer)

### Perform background correction on the fluorescent intensities
RG.bgcorrect <- backgroundCorrect(RG, method = 'normexp', offset = 16) #offset may depend on the number of sample (25%-50%)

### Normalize the data in array with the 'loess' method
RG.bgcorrect.norm <- normalizeWithinArrays(RG.bgcorrect, method = 'loess')
plotMD(RG.bgcorrect.norm)

### Normalize the data between arrays 
MA <- normalizeBetweenArrays(RG.bgcorrect.norm, method = "Aquantile")
plotDensities(RG.bgcorrect.norm)#check the signal density map
plotDensities(MA)
MA$genes$ControlType %>% unique() #include what category of probe

### remove the control probe
keepProbe <- MA$genes$ControlType == 0
MA.reomoveCTprobe <- MA[keepProbe,] #remove the control probe
MA.reomoveCTprobe$genes$ProbeName %>% table() %>% sort(decreasing = TRUE) %>% head(n = 3) 
MA.reomoveCTprobe.average <- avereps(MA.reomoveCTprobe, ID = MA.reomoveCTprobe$genes$ProbeName) # take an average of probe

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
u <- unique(targets2$Target)
f <- factor(targets2$Target, levels=u)
design <- model.matrix(~0+f) #make the design with each color of total sample
colnames(design) <- u
rownames(design) <- rownames(targets2)
corfit <- intraspotCorrelation(MA.reomoveCTprobe, design)
fit <- lmscFit(MA.reomoveCTprobe, design,correlation=corfit$consensus)

cont.matrix <- makeContrasts("wt.SC7-wt.C7","wt.SC7-wt.Control2","wt.C7-wt.Control2", levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
volcanoplot(fit2)
topTable(fit2)
glimmaVolcano(fit2)
glimmaVolcano(fit2,
  main = 'wt.SC7-wt.C7',
  xlab = "logFC",
  ylab = "negLog10PValue",
  html = NULL,
  width = 920,
  height = 920)
degs_temp <- topTable(fit2,coef=1, adjust="BH", n = Inf)%>%na.omit
logFC_cut = with(degs_temp, mean(abs(logFC))+2*sd(abs(logFC))) # take 95% CI for logFC
degs_temp = degs_temp %>%
  mutate( DEG = factor( ifelse( abs(degs_temp$logFC) > logFC_cut & degs_temp$adj.P.Val < 0.05,
                                ifelse(degs_temp$logFC > 0, 'up','down'),
                                'ns' ), 
                        levels = c("up", "ns", "down" ), ordered = F ) )%>%
  mutate( row = rownames(degs_temp))
results <- decideTests(fit2)
vennDiagram(results)














