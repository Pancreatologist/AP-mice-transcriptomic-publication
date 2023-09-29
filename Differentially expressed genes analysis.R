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
#following code is designed for different DE contrast, you have to change each group for different contrasts
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
contrastgroup <- 'C12-C7plus6'
contrastgroup <- 'FAEE150-FAEE50'
contrastgroup <- 'TLCS5-TLCS3'
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

### visualization for DE

volcanoplot<-ggscatter(degs_temp,x = "logFC", y = "logP",color = "Condition",palette = c("#CB9B8F", "gray" ,"#5B7493"  ),
                size=2,
                xlab = "Log2FC",
                ylab = "-Log10(Adj.P-value)",)+ 
  ggtitle(paste("Volcano plot for \n", contrastgroup))+
  ylim(0,15)+ 
  xlim(-3,3)+ 
  geom_hline(yintercept = 1.3 , linetype ="dashed")+ 
  geom_vline(xintercept = c(-1,1), linetype= "dashed")+
  labs(subtitle = paste(sprintf("up: %d", table(degs_temp$Condition)["Up"]),  sprintf("down: %d", table(degs_temp$Condition)["Down"]), sep = " | "))+
  theme(legend.position = "right")+
  theme_cowplot(20)

### upset plot for DE
# remove genes not signif. in any contrast
results <- results[rowSums(abs(results)) > 0, ]
# calculate the intersection between the differentially expressed gene sets
m <- make_comb_mat(abs(results), mode = "distinct")
# exclude self-intersects (total # of diff. genes will be displayed separately)
m <- m[comb_degree(m) > 1]
(ss <- set_size(m))
(cs <- comb_size(m))
ht <- UpSet(m, 
            set_order = colnames(m), comb_col = "#f9d099",
            comb_order = order(comb_degree(m)),
            top_annotation = HeatmapAnnotation(
              "Distinct diff. genes" = anno_barplot(
                cs, 
                ylim = c(0, max(cs)*1.1),
                border = FALSE, 
                gp = gpar(fill = "black"), 
                height = unit(4, "cm")
              ), 
              annotation_name_side = "left", 
              annotation_name_rot = 90),
            right_annotation = HeatmapAnnotation(
              which = "row",
              "Total" = anno_barplot(
                ss, 
                ylim = c(0, max(ss)*1.1),
                border = FALSE, 
                gp = gpar(fill = "black"), 
                width = unit(4, "cm")
              )
            ),
            column_title = "Intersection between contrasts"
)
ht = draw(ht)
od = column_order(ht)
rod = row_order(ht)
decorate_annotation("Distinct diff. genes", {
  grid.text(cs[od], 
            x = seq_along(cs), 
            y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("center", "bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"))
})
decorate_annotation("Total", {
  grid.text(ss[rod], 
            x = unit(ss[rod], "native") + unit(20, "pt"), 
            y = rev(seq_along(ss)), 
            default.units = "native", just = c("right", "bottom"), 
            gp = gpar(fontsize = 8, col = "#404040"))
})
