# First off, it's a good idea to have all your library() calls in one place at the start.
# This makes any errors from not having a package installed appear quickly, when you run the script.
library(BiocManager)
library(affy)
library(here)
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(ggsci)
library(limma)
library(ggfortify)
library(ggplot2)
library(ggsci)
library(pheatmap)# 热图
library(ggplotify) # 转换 pheatmap 对象为 ggplot2 对象
library(dplyr)
library(AnnoProbe)
library(hgu133plus2.db)
library(clusterProfiler)
library(enrichplot) 
library(DOSE)
library(randomForest)
library(ROCR)
library(STRINGdb) 
library(igraph)
library(magrittr)
###  set work dictionary, read files, do RMA and take the expression
library(here)
dir_cels='C:\\Users\\wudi\\Desktop\\Di Wu project\\msap+sap vs map\\Human AP Affymetrix Arrays\\cel_files'
eset = ReadAffy(celfile.path=dir_cels)
raw.names <- sampleNames(eset) # change sample name
eset_rma = rma(eset) # RMA
# Make sure you're aware of what microarray type your data come from, and what level rma() is summarising at. For example, some Affymetrix arrays have an option to summarise at an 'exon' or a 'probeset' level with this step.
exprset1 = exprs(eset_rma) # take the expression data from CEL files

### look over the expr
summary( as.numeric( unlist(exprset1) ) ) # whether is after log2 (used to be is log2)
table( is.na(exprset1) ) # chech if there is NA
table( exprset1 < 0 ) # check if there is data less than 0
hist(unlist( exprset1 ))

### 1.1 get the sampleinfo----
pData(eset)
sampleinfo = pData(eset) # Where did this pData come from? This is the first time it appears, should this line come after another step? I wanna to use the pData to check the information of group and get the clear group information
term_x = sampleinfo$sample 
table(term_x)
sort(c( "MSAPandSAP","MAP") )
level_x = c( "MSAPandSAP","MAP")
contrast = c( "MSAPandSAP-MAP" )
b <- read.csv(file = 'APnumber.csv')#read the grouplist
term_x <- b$RAC
sampleinfo = sampleinfo %>% 
  mutate( Group = b$RAC )
identical( rownames(sampleinfo), colnames(exprset1))

### 1.2 look over the total data
boxplot( exprset1 )
library(ggsci)
n_select = 5000  
sample( 1:nrow(exprset1), n_select )
expr_l = exprset1[ sample(1:nrow(exprset1),n_select) , ] %>% 
  as.data.frame(  ) %>% 
  gather(  )  %>%  
  setNames( c( 'Sample','Expression') ) %>% 
  mutate( Group = sampleinfo$Group[ match( Sample, sampleinfo$Group) ] ) %>%  
  arrange( Group ) %>% 
  mutate( Sample = factor(Sample, levels = unique(Sample), ordered = T ) ) 
unique(expr_l$Sample)
p1 = ggplot(expr_l,aes(x= Sample, y= Expression, fill= Group, color = Group))+ 
  geom_boxplot(outlier.shape = 0.05,outlier.size = 0.05) +
  theme_light() +
  ggtitle("Expression levels of samples") +
  theme(
    plot.title = element_text( hjust = 0.5 ),
    axis.text.x = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_nejm(alpha = 1)+
  scale_fill_nejm(alpha = 0.8)

# normalize based on the boxplot if necessary, based on the  rma's quantile normalisation
#exprset1 = normalizeBetweenArrays( exprset1, method = "quantile" ) # normalizeBetweenArrays 可用于芯片数据标准化

### 1.3 PCA
pca1 = prcomp( as.data.frame( t( exprset1) ) ) 
xxx = pca1[["x"]]
p3 = autoplot( pca1,
               data = sampleinfo, colour = "Group",
               size=2.5, frame = FALSE, frame.type = 'norm')+
  theme_light() +
  ggtitle("PCA of samples") +
  theme(
    plot.title = element_text( hjust = 0.5 ),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_nejm(alpha = 0.8) +
  scale_fill_nejm(alpha = 0.8) 
p3
ggsave(p3, filename =paste("PCA for XXX",".pdf",sep=""),width = 5.2, height = 4)
rm( pca1 )

### 1.4.1 probe annotation----
ids=toTable(hgu133plus2SYMBOL)#take the probe ID and corespending Gene from hgu133plus2.db
gpl_anno = ids %>% 
  as.data.frame() %>%
  distinct( probe_id, .keep_all = T ) # delete the dulpucated probe 
rownames(gpl_anno) = gpl_anno$probe_id 
head( gpl_anno )
# clean the meaningless probe
exprset2 = exprset1[rownames(exprset1) %in% gpl_anno$probe_id,] # delete the unmatched probe expre
dim(exprset2)
gpl_anno = gpl_anno[gpl_anno$probe_id %in% rownames(exprset2),] # delete the unmatched probe

### 1.4.2 expression annotation----
identical( gpl_anno$probe_id, rownames(exprset2))# Check the rownames matched or not
gpl_anno = gpl_anno[rownames(exprset2),] #rearrangement

identical( gpl_anno$probe_id, rownames(exprset2))# Check again, if true can cbind
exprset2 = cbind(gpl_anno,exprset2)
table(duplicated(exprset2$symbol)) #check the duplicated gene symbol
annoted_exprsetwithGeneid <- aggregate(x = exprset2[,3:ncol(exprset2)],
                             by = list(exprset2$symbol),
                             FUN = median)  # give the duplicated genes as mean expression
table(duplicated(annoted_exprsetwithGeneid$symbol)) #check  if all the false
rownames(annoted_exprsetwithGeneid) <- annoted_exprsetwithGeneid[,1] # make the first col to the rownames
annoted_exprsetwithGeneid <- as.matrix(annoted_exprsetwithGeneid[,-1]) # delete the symbol col and change dataframe to the matrix


### 1.5 see the total difference of genes
# I would also recommend filtering transcripts by signal intensity prior to any further analysis or annotation. Answer:I do the filter as above, please check if right
# histograms can help here. Hist what?
top500_sd = annoted_exprsetwithGeneid %>%
  as.data.frame() %>%
  mutate( SD = apply(., 1, sd ) ) %>%
  # arrange( desc(SD) ) %>%
  slice_max( SD, n= 400 )
legend_col = data.frame( row.names = rownames(sampleinfo),
                         Group = sampleinfo$Group )  # set the grouplabel for the heatmap
bk = 2 
brk <- c(seq(-bk,-0.01,by=0.01),seq(0,bk,by=0.01))
heat1 = pheatmap(sd,
                 scale = "row", # 
                 annotation_col = legend_col,
                 color = c(colorRampPalette(colors = c("dodgerblue4","white"))(length(brk)/2),
                           colorRampPalette(colors = c("white","brown"))(length(brk)/2)),
                 legend_breaks=seq(-bk,bk,1),
                 breaks=brk,
                 border_color = NA,
                 show_rownames = F,
                 show_colnames = F)
p4 = as.ggplot(heat1)+
  ggtitle('Top 500 gene with large SD')+  
  xlab('Sample') + ylab('Gene')+
  theme(
    plot.title = element_text(hjust = 0.4),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) 
p4

rm(sd, legend_col, bk, brk, heat1 )

### 2 DEG
### 2.1 make the contrast----
# group
table(sampleinfo$Group)
contrast #check the righ contrast

# contrast matrix
design <- model.matrix(~0+factor(sampleinfo$Group) ) %>%
  as.data.frame() %>% 
  setNames( unique(sampleinfo$Group) )
rownames(design) = rownames(sampleinfo)
design

c("MSAPandSAP","MAP-MSAP+SAP")#contrast

contrast.matrix <-makeContrasts( contrasts = c("MSAPandSAP-MAP","MAP-MSAPandSAP"), levels = design ) 
contrast.matrix

### 2.2 count the DEG----
degs_temp <- lmFit(annoted_exprsetwithGeneid, design ) %>% # Fit linear model
  contrasts.fit( contrast.matrix ) %>% # compute estimated coefficients and standard errors for a given set of contrasts
  eBayes %>% # compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression
  topTable( coef = 1, n = Inf ) %>% # this mean the first contrast
  na.omit %>% # 
  setNames( c("log2FC", "Mean.Expr", "t", "P.value", "adj.P", "B" ) ) %>%
  mutate( probe_id = rownames(.) ) 

logFC_cut = with(degs_temp, mean(abs(log2FC))+2*sd(abs(log2FC))) # take a 95%CI for log2FC

# There are problems here — the P.value from topTable is not adjusted for multiple comparisons,
# meaning each transcript coefficient is being tested entirely independently of the 1000s of other tests.
# I would recommend using adj.P.Val instead, which are adjusted to control the proportion of Type I errors (false discovery rate, FDR). Answer:I change the P to adj.P in the degs_temp follow this comment. And there are all ns. I don't know if doing correctly. 
# In this case, adj.P.Vals of <0.05 would reflect a FDR of <0.05.
degs_temp = degs_temp %>%
  mutate( DEG = factor( ifelse( abs(degs_temp$log2FC) > logFC_cut & degs_temp$adj.P < 0.05,
                                ifelse(degs_temp$log2FC > 0, 'up','down'),
                                'ns' ), 
                        levels = c("up", "ns", "down" ), ordered = T ) )
table( degs_temp$DEG )


rm( design, contrast.matrix )

### 2.3 volcano plot for DEG
p5 = ggplot(data = degs_temp,
            aes(x = log2FC, y = -log10(adj.P), color = DEG)) +
  # ylim(0,25)+ 
  #  xlim(-0.005,0.005)+ 
  geom_point(size = 1, shape = 16, alpha = 0.9) +
  ggtitle('Overall distribution of DEGs') +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_manual(values = c('brown','grey','dodgerblue3' ))+
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_vline( xintercept = c(-logFC_cut, logFC_cut), linetype = "dashed" ) +
  geom_hline( yintercept = -log10(0.05), linetype = "dashed" )
p5
ggsave( p5, filename = paste("2.3 volcano plot",".pdf",sep= ""),width = 5, height = 4)

rm( logFC_cut )

### 2.4 heatmap----
chs_x = degs_temp %>%
  slice_max( abs(log2FC), n = 500 ) 

heat_matrix <- as.matrix(row.names(degs_temp))
legend_col = data.frame( row.names = rownames(sampleinfo),
                         Group = sampleinfo$Group )  
bk = 2 
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
                 show_rownames = F,
                 show_colnames = F
)
p6 = as.ggplot(heat1)+
  ggtitle('Top 500 DEGs')+  
  xlab('Sample') + ylab('Gene')+
  theme(
    plot.title = element_text(hjust = 0.4),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
p6
rm(chs_x, heat_matrix, legend_col, bk, brk, heat1 )


### Gene ID change (from "SYMBOL" to "ENTREZID") to further using 'clusterProfiler'
degswithENTRE = degs_temp %>%
  mutate(symbol = row.names(degs_temp))%>%
  pull(symbol) %>% 
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), # from  'symbol' to 'ENTREZID'
        OrgDb = hgu133plus2.db ) %>%  
# if DEG too many, that filter to 300, make a threshold
#qu1 = quantile(degs_final$Mean.Expr, probs = 0.25 ) # expr : IQR 25%
#fc1 = mean(abs(degs_final$log2FC))  # Log2FC :  mean
#qu1;fc1
# degswithENTRE = filter( Mean.Expr > qu1 & abs(log2FC) > fc1 & adj.P < 0.05 ) %>%  # expr > IQR 25%，Log2FC > mean，adj.P < 0.05
#  arrange( adj.P ) %>% # arrange by adj.P
#  arrange( desc( abs(log2FC) ) ) %>%
#  mutate(symbol = row.names(degs_temp))%>%
# pull( symbol ) %>% # 
#  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), 
#        OrgDb = org.Mm.eg.db ) %>%  # from  'symbol' to 'ENTREZID'
#  slice_head( n = 300 ) # get the top 300 or other


### 3 enrich based on DEG

### 3.1 GOenrich
go_data <- enrichGO( gene = degswithENTRE$ENTREZID, 
                     OrgDb="hgu133plus2.db", 
                     ont ="ALL", 
                     pvalueCutoff = 0.05, 
                     qvalueCutoff = 1,
                     readable= TRUE 
) %>% 
  clusterProfiler::simplify( . ) %>% 
  pairwise_termsim() 

xxx = go_data@result
write.csv(xxx, file = 'GO-richment.csv')
### 3.2 KEGG enrich----
kegg_data <- enrichKEGG(gene = degswithENTRE$ENTREZID, 
                        keyType = "kegg",
                        organism   = 'hsa', # mmu is mouse，human is  hsa
                        pvalueCutoff = 1, # 1 means with all the result without filter
                        qvalueCutoff = 1,
                        use_internal_data = F 
) %>%
  setReadable( ., OrgDb = 'hgu133plus2.db', keyType = 'ENTREZID' )  %>% # from ENTRZID to SYMBOL
  pairwise_termsim() # to do the further emap 
write.csv(kegg_data, file = 'KEGG-richment.csv')

### 3.3 visualization for enrich----
enrichDO()
enrichWP()
dotplot(kegg_data, showCategory= 20, title="Top enriched terms") + 
  theme_light()

emapplot(go_data, showCategory = 20) +
  ggtitle("Top enriched terms")+
  theme_light()+
  theme(axis.text = element_blank(), axis.title = element_blank())+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

cnetplot(go_data,showCategory = 10, 
         circular=F, 
         colorEdge= T)+
  ggtitle("Terms and DEGs") +
  theme_light()+
  theme(axis.text = element_blank(), axis.title = element_blank())

### 3.4.1 make a list for GSEA
# go_gmt = read.gmt("/usr/local/lib/R/Rlib/c5.go.v7.2.symbols.gmt")
# kegg_gmt = read.gmt("/usr/local/lib/R/Rlib/c2.cp.kegg.v7.2.symbols.gmt")
btr_gsea_order = bitr( degs_temp$symbol, 
                       fromType = "SYMBOL", toType = c("ENTREZID"),
                       OrgDb = hgu133plus2.db) %>% # change the GENE ID
  mutate( log2FC = degs_temp$log2FC[ match(SYMBOL, degs_temp$symbol ) ] ) %>% # mutate the col of log2FC
  arrange( desc( log2FC ) ) # from high to low

table( is.na(btr_gsea_order$ENTREZID) )  # check na in ENTREZID
table( is.na(btr_gsea_order$log2FC) )  # check na in log2FC
head(btr_gsea_order) ; tail(btr_gsea_order) # check na in the order


gsea_list = setNames( btr_gsea_order$log2FC, btr_gsea_order$ENTREZID ) # use log2FC to make the list

head(gsea_list);tail(gsea_list)
### 3.4.2 GSEA GO enrich
gsea_go <- gseGO(gsea_list, 
                 ont = "ALL",
                 OrgDb = hgu133plus2.db, #human
                 eps = 0 ) %>%
  clusterProfiler::simplify( )  



### 3.4.3 GSEA KEGG enrich
gsea_kegg <- gseKEGG(gsea_list, 
                     keyType = "kegg",
                     organism = 'hsa', # human
                     pvalueCutoff = 1,
                     eps = 0,
                     use_internal_data = F ) 
# gsea_wp <- gseWP(gsea_list, organism = 'Homo sapiens', eps = 0)
# gsea_do <- gseDO(gsea_list, by = "fgsea", eps = 0)


### 3.4.4 visualization for GSEA GO and GSEA KEGG
ridgeplot(gsea_go,showCategory = 30)+
  theme_light()+
  ggtitle("Top GSEA-GO terms")+
  theme(
    plot.title = element_text(hjust = 0.4),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_fill_gradient( low = "brown", high = "blue" )

ridgeplot(gsea_kegg,showCategory = 30)+
  theme_light()+
  ggtitle("Top GSEA-KEGG terms")+
  theme(
    plot.title = element_text(hjust = 0.4),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_fill_gradient( low = "brown", high = "blue"  )

enrichplot::gseaplot2(gsea_go,1:5)+
  theme_light()+
  ggtitle("Top 5 GSEA-GO terms")+
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )


