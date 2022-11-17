###return the expression data of each color from two color array, i think which may be correct but need to valid.
#Gordon Smyth said that we do not need to disassemble the limma data objects.https://www.biostars.org/p/492110/
targets2 <- targetsA2C(target)
eset = exprs.MA(MA.reomoveCTprobe.average) #get the expression from MA. In the guide, I think it will from aCY3, aCy5, bCy3, bCy5...... So it may be arranged as targets2 list? 
colnames(eset) <- rownames(targets2) # If it arranged correctly, this code may be right. But if it doesnt arranged as targets, this is the thing i dont sure.
rownames(eset) <- MA.reomoveCTprobe.average$genes$GeneName #make the rowname to the genenames from the rawdata
plot(MA.reomoveCTprobe$A)
plot(eset)
eset <- eset[row.names(eset) %in% row.names(Annotationname), ] #filter eset to only mRNA based on filter genes which i do it in the pre-process(excluding some begining with 'chr')


### Gene ID change (from "SYMBOL" to "ENTREZID") to further using 'clusterProfiler'
gpl_anno=read.table('GPL10787.txt', sep = "\t", quote = "\"",header = TRUE, fill = TRUE)  # read the annotation file for the GPL10787 platform (downloaded)
head(gpl_anno)
gpl_anno = gpl_anno %>% distinct(GENE_SYMBOL, .keep_all = TRUE ) # delete the repeat probe
rownames(gpl_anno) = gpl_anno$GENE_SYMBOL

### annotation for eset
eset = eset[rownames(eset) %in% gpl_anno$GENE_SYMBOL,] #  delete the unmatched in eset
dim(eset)
gpl_anno = gpl_anno[gpl_anno$GENE_SYMBOL %in% rownames(eset),] #  delete the unmatched in annotation
identical( gpl_anno$GENE_SYMBOL, rownames(eset)) #check the arrange of eset and annotation
gpl_anno = gpl_anno[rownames(eset),] #  arrange 
identical( gpl_anno$GENE_SYMBOL, rownames(eset)) #recheck
eset2 = cbind(gpl_anno,eset)

### remove batch effect, dont use
#batch <- targets2$Batch
#eset <- removeBatchEffect(eset, batch)

### make the sampleinfo
sampleinfo <- targets2$Target %>% as.data.frame()%>%
  mutate(rownames(targets2)) %>%
  `colnames<-`(c('group','filenames'))
rownames(sampleinfo) <- sampleinfo$filenames

### check the PCA of sample
PCA <- prcomp(t(eset), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Model =
                       targets2$Target)
ggplot(dataGG, aes(PC1, PC2)) +
  #geom_point(aes(shape = Model)) +
  geom_point(aes(colour = Model)) +
  ggtitle("PCA (cluster experimental groups)") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  #scale_shape_manual(values = c(1:5))
  scale_colour_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000','#CCFFFF', '#CCFFCC',	'#CCFF99','#CCFF66','#CCFF33','#CCFF00')) 

### DEG
### firstly check the PCA of each comtrast, because total PCA may difficult to see the detail of subgroups 
x <- c("TLCS","SC")# choose the target col which filename cooresponding to group in target
eset_new = eset[, targets2$Target %in% x ] 
sampleinfonew <- sampleinfo[sampleinfo$group%in%x,]
PCA <- prcomp(t(eset_new), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Model =
                       sampleinfonew$group)
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = Model)) +
  ggtitle("PCA (cluster experimental groups)") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(1:6))

### correlation between samples
r <- cor(eset_new,method = "pearson") # methods have "pearson", "spearman", "kendall"
pheatmap(r, 
         show_colnames = TRUE,   
         show_rownames=TRUE,     
         fontsize=5,             
         color = colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(50), 
         annotation_legend=TRUE, 
         border_color=NA,        
         scale="none",           
         cluster_rows = TRUE,    
         cluster_cols = TRUE)

### make contrast for DEG
u <- unique(targets2$Target)
f <- factor(targets2$Target, levels=u)
Batch <- targets2$Batch
design <- model.matrix(~0+f) #make the design with each color of total sample
colnames(design) <- u
rownames(design) <- rownames(targets2)
cont.matrix <- makeContrasts(wt.Control-wt.ControlREHY,wt.Control-wtControl,wtControl-wt.ControlREHY,levels=design) #it will be based on what need to contrast 
fit2 <- lmFit(eset, design) %>% # Fit linear model
  contrasts.fit( cont.matrix ) %>% # compute estimated coefficients and standard errors for a given set of contrasts
  eBayes # compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression
results <- decideTests(fit2)
vennDiagram(results)

degs_temp <- topTable(fit2,coef=1, adjust="BH", n = Inf)%>%na.omit # get the detail of different DE
logFC_cut = with(degs_temp, mean(abs(logFC))+2*sd(abs(logFC))) # get a 95%CI for logFC
degs_temp = degs_temp %>%
  mutate( DEG = factor( ifelse( abs(degs_temp$logFC) > 1 & degs_temp$adj.P.Val < 0.05,
                                ifelse(degs_temp$logFC > 0, 'up','down'),
                                'ns' ), 
                        levels = c("up", "ns", "down" ), ordered = F ) )%>%
  mutate( row = rownames(degs_temp)) # for the volcano plot
summary(degs_temp)
head(degs_temp)
table(degs_temp$DEG) #to see the number of up and down genes


### volcano plot 
#first choice for volcano plot
p1 <- gradual_volcano(degs_temp, x = "logFC", y = "P.Value", # because adjP plot may unsatisfied, so i plot it with P-value.
                      fills = brewer.pal(5, "RdYlBu"),
                      colors = brewer.pal(8, "RdYlBu"),
                      x_lab = 'logFC',
                      y_lab = '-log10(P.Value)',
                      legend_title = '-log10',
                      pointSizeRange = c(0.5, 4),
                      label = "row", label_number = 10, output = FALSE) # with 10 genes with names
p1 
save(degs_temp, file = paste('wt.C7plus6-wt.SC7plus6','.Rdata',sep = ""))
#second choice for volcano plot
degs_temp$label <- c(rownames(degs_temp)[1:10],rep(NA,(nrow(degs_temp)-10))) #get the label need to clear
p1 <- ggplot(degs_temp,aes(logFC, -log10(adj.P.Val)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  geom_point(aes(size=-log10(adj.P.Val), color= -log10(adj.P.Val)))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+    # 指定颜色渐变模式：
  scale_size_continuous(range = c(1,3))+    # 指定散点大小渐变模式：
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.7),
        legend.justification = c(0,1)    # 调整主题和图例位置：
  )+
  guides(col = guide_colourbar(title = "-Log10_q-value"),
         size = "none")+    # 设置部分图例不显示：
  geom_text(aes(label=label, color = -log10(adj.P.Val)), size = 3, vjust = 1.5, hjust=1)+    # 添加标签：
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")    # 修改坐标轴：

###heatmap for DEG validation
chs_x = degs_temp %>%
  slice_max( abs(logFC), n = 50 ) 
heat_matrix <- eset_new [row.names(eset_new) %in% row.names(chs_x), ] 
legend_col = data.frame(row.names = colnames(eset_new),
  Group= sampleinfonew$group)  # set the label of group
rownames(legend_col) <- colnames(heat_matrix)
heat1 <- pheatmap(heat_matrix, scale="row", name = 'Expression level',
                  cellwidth = 10,cellheight = 5, # set the high and wide of box
                  border="white",
                  annotation_col = legend_col,
                  annotation_colors = list(Group=(c(wt.C7="#1B9E77",wt.Control="#D95F02"))),
                  cluster_rows = F,# del the cluster of genes
                  cluster_cols = TRUE,
                  show_rownames = F, #delete the row and col id
                  show_colnames = F, 
                  clustering_distance_rows = "minkowski", # set the class of cluster
                  color = c(colorRampPalette(colors = c("dodgerblue4","white"))(length(brk)/2),
                            colorRampPalette(colors = c("white","brown"))(length(brk)/2)),
                  #cutree_cols = 3, #划分列
                  #cutree_rows =5, #划分行
                  #cluster_cols = TRUE,treeheight_col = 20, #按照列（样本名）划分区域
                  #cluster_rows = T,treeheight_row = 20,按照行（基因名）划分区域
                  fontsize_row = 12, # set the size of font
                  fontsize_col = 16)
heat1
p2 = as.ggplot(heat1)+
  ggtitle('Top 500 DEGs')+  
  theme(
    plot.title = element_text(hjust = 0.4),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
p2
