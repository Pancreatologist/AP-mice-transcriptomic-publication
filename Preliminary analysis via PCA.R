###return the expression data of each color from two color array
targets2 <- targetsA2C(target)
eset = exprs.MA(MA.reomoveCTprobe) #get the expression from MA
colnames(eset) <- rownames(targets2)
rownames(eset) <- MA.reomoveCTprobe$genes$GeneName 
eset <- removeBatchEffect(eset, as.factor(targets2$Batch))

PCA_raw <- prcomp(t(eset), scale. = T) 
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
x <- c('C4','C4plus9','C7','C7plus6','C12','FAEE150','FAEE50','TLCS3','TLCS5')
targets3 <- targets2[targets2$Target%in% x,]
dataGG <- data.frame(PC1 = PCA_raw$x[targets2$Target%in% x,1], PC2 = PCA_raw$x[targets2$Target%in% x,2],
                     Batch =as.factor(targets3$Batch),
                     Group = as.factor(targets3$Target))
dataGG <- dataGG %>% mutate(Model= ifelse(grepl('TLCS',dataGG$Group),'TLCS',
                                          ifelse(grepl('FAEE',dataGG$Group),'FAEE','CER')))
#dataGG <- dataGG %>% mutate(Model= ifelse(grepl('TLCS',dataGG$Group),'TLCS',
#                                          ifelse(grepl('FAEE',dataGG$Group),'FAEE',
#                                                 ifelse(grepl('Control',dataGG$Group),'Control','CER'))))
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Batch,color=Group)) +
  ggtitle("PCA plot") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(1:5)) +
  scale_color_manual(values = c('#e6194B','#fabed4','#f58231','#ffd8b1','#bfef45',"#BFBFC1",'#42d4f4','#4363d8','#469990','#808000'))+
  stat_ellipse(aes(fill=Model),type="t",geom="polygon",alpha=0.3)+
  theme(panel.grid.major = element_blank(),  # 隐藏主要网格线
        panel.grid.minor = element_blank(),  # 隐藏次要网格线
        panel.background = element_blank(),  # 设置空白背景
        axis.line = element_line(color = "black"))
