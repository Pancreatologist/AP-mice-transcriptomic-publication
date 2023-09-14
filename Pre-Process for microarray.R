###Read in and normalise data
target <- readTargets("designwithbatch.txt" ) 
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
RG.bgcorrect <- backgroundCorrect(RG, method = 'normexp', offset = 50) #offset may depend on the number of sample (25%-50%)

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

### check the batch effect
negativeprobe <- exprs.MA(a) 
targets2 <- targetsA2C(target)
colnames(negativeprobe) <- rownames(targets2)
rownames(negativeprobe) <- a$genes$GeneName
my_data <- data.frame(
  Batch = as.factor(targets2$Batch),
  value = colMeans(negativeprobe))
ggplot(my_data, aes(x = Batch, y = value)) +
  geom_boxplot()

### For replicate gene in each sample, replace values with max one
Amean <- rowMeans(MA.reomoveCTprobe$A)
o <- order(Amean, decreasing=TRUE)
MA.reomoveCTprobe <- MA.reomoveCTprobe[o,]
dup <- duplicated(MA.reomoveCTprobe$genes$GeneName)
MA.reomoveCTprobe <- MA.reomoveCTprobe[!dup,] 

### filter the mRNA
MA.reomoveCTprobe$genes <- MA.reomoveCTprobe$genes %>% 
  mutate(ENTREZID = bitr(MA.reomoveCTprobe$genes$GeneName, fromType = 'SYMBOL',OrgDb=org.Mm.eg.db, toType = 'ENTREZID',drop = F)$ENTREZID)
MA.reomoveCTprobe <- MA.reomoveCTprobe[!is.na(MA.reomoveCTprobe$genes$ENTREZID),]
