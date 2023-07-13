### DE analyse, include the batch effect into designmatrix
targets2 <- targetsA2C(target)
table(targets2$Target)
u <- unique(targets2$Target)
f <- factor(targets2$Target, levels=u)
b <- factor(targets2$Batch, levels=unique(targets2$Batch))
design <- model.matrix(~b+f) #make the design with each color of total sample
is.fullrank(design)
nonEstimable(design)
col_names <- colnames(design) %>% gsub("^f", "", .)
colnames(design) <- col_names
colnames(design)[1] <- "Intercept"
rownames(design) <- rownames(targets2)
corfit <- intraspotCorrelation(MA.reomoveCTprobe, design)
fit <- lmscFit(MA.reomoveCTprobe, design,correlation=corfit$consensus) 
table(targets2$Target)
#following code is designed for different DE contrast
contrastgroup <- 'C7plus6-Control'
contrastgroup <- 'SC7plus6-Control'
contrastgroup <- 'C7-Control'
contrastgroup <- 'C12-Control'
contrastgroup <- 'C4plus9-Control'
contrastgroup <- 'C4-Control'
contrastgroup <- 'EtOH-Control'
contrastgroup <- 'FAEE50-Control'
contrastgroup <- 'FAEE150-Control'
contrastgroup <- 'TLCS3-Control'
contrastgroup <- 'SalinePerf-Control'
contrastgroup <- 'TLCS5-Control'
contrastgroup <- 'ppifC7-ppifControl'
contrastgroup <- 'ppifC7plus6-ppifControl'
contrastgroup <- 'ppifTLCS-ppifControl'
contrastgroup <- 'ppifPOAFAEE-ppifControl'
contrastgroup <- 'ppifControl-Control'
contrastgroup <- '(ppifC7plus6-C7plus6-ppifControl+Control)/2'
contrastgroup <- '(ppifC7−C7−ppifControl+Control)/2'
contrastgroup <- '(ppifPOAFAEE-POAFAEE-ppifControl+Control)/2'
contrastgroup <- '(ppifTLCS-TLCS-ppifControl+Control)/2'
contrastgroup <- 'C12-C7plus6'
contrastgroup <- 'FAEE150-FAEE50'
contrastgroup <- 'TLCS5-TLCS3'
contrastgroup <- 'ppifTLCS-TLCS'
contrastgroup <- 'ppifC7plus6-C7plus6'
contrastgroup <- 'ppifPOAFAEE-POAFAEE'
contrastgroup <- 'FAEE150-C7plus6'
contrastgroup <- 'FAEE150-TLCS3'
contrastgroup <- 'C7plus6-TLCS3'
contrastgroup <- 'C7-C7plus6'
cont.matrix <- makeContrasts(contrastgroup,levels=design)
#compare different models via venn plot
cont.matrix <- makeContrasts(C7plus6-Control,FAEE150-Control,TLCS3-Control,levels=design)
cont.matrix <- makeContrasts(C12-C7plus6,TLCS5-TLCS3,FAEE150-FAEE50,levels=design)
# then get the DE results
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2,adjust.method='BH',lfc = 1) #logFC more than 1
vennDiagram(results,include=c("up","down"),counts.col=c("#CB9B8F", "#5B7493"))
degs_temp <- topTable(fit2,coef=1, adjust="BH", n = Inf,sort.by = 'logFC')%>%na.omit
hist(degs_temp$P.Value, col = brewer.pal(3, name = "Set2")[2], xlab = "p-values")# check the distribution of P-value
row.names(degs_temp)=degs_temp$GeneName
degs_temp = degs_temp %>%
  mutate( Condition = factor( ifelse( abs(degs_temp$logFC) > 1 & degs_temp$adj.P.Val < 0.05,
                                ifelse(degs_temp$logFC > 0, 'Up','Down'),
                                'NS' ), 
                        levels = c("Up", "NS", "Down" ), ordered = F ))%>%
  mutate( row = rownames(degs_temp))
summary(degs_temp)
degs_temp$logP <- -log10(degs_temp$adj.P.Val)
write.csv(degs_temp[c('GeneName','Description','ENTREZID','logFC','AveExpr','t','P.Value','B','Condition')],'DEG.csv')#export each DE results for the Table S1


