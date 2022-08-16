# First off, it's a good idea to have all your library() calls in one place at the start.
# This makes any errors from not having a package installed appear quickly, when you run the script.

#  set work dictionary ----

rm( list = ls() ) # I would recommend NOT saving the R environment when closing the session. This will make this line unnecessary.
getwd # Instead of using setwd(), I would recommend using the approach given here: https://www.tidyverse.org/blog/2017/12/workflow-vs-script/. I use the 'here' package a lot.
setwd("C:\\Users\\wudi\\Desktop\\Di Wu project\\msap+sap vs map") 
library(affy) # If this works for your microarray type then carry on, but if you have problems the 'oligo' package is the more modern version of 'affy', and may be more compatible with newer arrays.
dir_cels='C:\\Users\\wudi\\Desktop\\Di Wu project\\msap+sap vs map\\Human AP Affymetrix Arrays\\cel_files'
eset = ReadAffy(celfile.path=dir_cels)

# Are lines 14-18 necessary? Or just experimenting? They seem to lead to a redundant result.
raw.names<-sampleNames(eset)  
raw.names <- a$x
write.csv(sampleNames(eset), file = 'a.txt')
a <- read.csv(file = 'a.txt')
sampleNames( eset ) = raw.names # change sample name

### RMA and take the expression
eset_rma = rma(eset) # RMA
# Make sure you're aware of what microarray type your data come from, and what level rma() is summarising at. For example, some Affymetrix arrays have an option to summarise at an 'exon' or a 'probeset' level with this step.
exprset1 = exprs(eset_rma) # take the expression data from CEL files

### look over the expr
summary( as.numeric( unlist(exprset1) ) ) # whether is after log2 (used to be is log2)
table( is.na(exprset1) ) # chech if there is NA
table( exprset1 < 0 ) # check if there is data less than 0
hist(unlist( exprset1 ))

# 1.1 get the sampleinfo----
sampleinfo = pData(eset) # Where did this pData come from? This is the first time it appears, should this line come after another step?
# get the clear group information
library(tidyverse)
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

### 1.3 look over the total data
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


# normalize based on the boxplot
library(limma)
# normaliseBetweenArrays should not be necessary, quantile normalisation is performed as part of rma().
# Could you send me an image of the boxplot from line 48 please? This should show whether rma's quantile normalisation worked.
exprset1 = normalizeBetweenArrays( exprset1, method = "quantile" ) # normalizeBetweenArrays 可用于芯片数据标准化

### 1.5 PCA
library(devtools)
install_github('sinhrks/ggfortify') # install commands should be left out of scripts, and only used in console. Stick to just 'library' calls, and the end user can install the package themselves if necessary.
library(ggfortify); library(ggplot2)
library(ggplot2) # This shou;d already be loaded from library(tidyverse) earlier
library(ggsci) # You have already loaded this package
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

### see the total difference of genes
# I assume you mean 'standard deviation' rather than "difference" here? So these are the top 500 most variable transcripts?
# Another point is that these aren't necessarily "genes" yet. It depends on what this microarray targets,
# and what level rma() summarised at for these data. E.g. these rows could represent 'probesets', 'exons', 'transcript clusters' etc.
# Not all of these would match the definition of a gene. 
# Your arrays here appear to be Human Genome U133 Plus 2.0 (https://www.affymetrix.com/products_services/arrays/specific/hgu133plus.affx#1_2)
# which cover 47000 transcripts, including genes and variants. This is obviously a lot more than we would expect if these were protein-coding genes alone.

# If you would like to describe the output of the following analyses as "genes",
# you should do your probe annotation from section 4 HERE.
# Otherwise, what you are finding and clustering are differentially expressed TRANSCRIPTS, not genes.
# You will also have multiple transcripts from the same gene in the analysis, which might not be what you want.

# I would also recommend filtering transcripts by signal intensity prior to any further analysis or annotation,
# histograms can help here.

sd = exprset1 %>% # As 'sd()' is a function in R, this call will overwrite that function. Perhaps name it something like top500_sd
  as.data.frame() %>%
  mutate( SD = apply(., 1, sd ) ) %>%
  # arrange( desc(SD) ) %>%
  slice_max( SD, n= 500 ) %>% 
  dplyr::select( 1:ncol(exprset1) ) # Is this step necessary?


library(pheatmap) 
library(ggplotify) 
legend_col = data.frame( row.names = rownames(sampleinfo),
                         Group = sampleinfo$Group )  
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
                 show_colnames = F
)
p4 = as.ggplot(heat1)+
  ggtitle('Top 500 differentiated gene ')+  
  xlab('Sample') + ylab('Gene')+
  theme(
    plot.title = element_text(hjust = 0.4),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) 
p4

rm(sd, legend_col, bk, brk, heat1 )

  ## 3、DEG----
library( limma )
### 3.1 make the contrast----
# group
table(sampleinfo$Group)
contrast


# contrast matrix
design <- model.matrix(~0+factor(sampleinfo$Group) ) %>%
  as.data.frame() %>% 
  setNames( unique(sampleinfo$Group) )
rownames(design) = rownames(sampleinfo)
design

c("MSAPandSAP","MAP-MSAP+SAP")#contrast

contrast.matrix <-makeContrasts( contrasts = c("MSAPandSAP-MAP","MAP-MSAPandSAP"), levels = design ) 
contrast.matrix

### 3.2 count the deg-probes----
degs_temp <- lmFit( exprset1, design ) %>% # Fit linear model
  contrasts.fit( contrast.matrix ) %>% # compute estimated coefficients and standard errors for a given set of contrasts
  eBayes %>% # compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression
  topTable( coef = 1, n = Inf ) %>% # this mean the first contrast
  na.omit %>% # 
  setNames( c("log2FC", "Mean.Expr", "t", "P.value", "adj.P", "B" ) ) %>%
  mutate( probe_id = rownames(.) ) 

logFC_cut = with(degs_temp, mean(abs(log2FC))+2*sd(abs(log2FC))) # 取一个 log2FC 的 95% 置信区间阈值

# There are problems here — the P.value from topTable is not adjusted for multiple comparisons,
# meaning each transcript coefficient is being tested entirely independently of the 1000s of other tests.
# I would recommend using adj.P.Val instead, which are adjusted to control the proportion of Type I errors (false discovery rate, FDR).
# In this case, adj.P.Vals of <0.05 would reflect a FDR of <0.05.

degs_temp = degs_temp %>%
  mutate( DEG = factor( ifelse( abs(degs_temp$log2FC) > logFC_cut & degs_temp$P.value < 0.05,
                                ifelse(degs_temp$log2FC > 0, 'up','down'),
                                'ns' ), 
                        levels = c("up", "ns", "down" ), ordered = T ) )
table( degs_temp$DEG )


rm( design, contrast.matrix )

### 3.3 volcano plot----
p5 = ggplot(data = degs_temp,
            aes(x = log2FC, y = -log10(P.value), color = DEG)) +
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
ggsave( p5, filename = paste("3.3_DEG火山图_",".pdf",sep= ""),width = 5, height = 4)

rm( logFC_cut )

# just take 500 genes, base on log2FC
# Make sure all 500 genes meet adj.P-value < 0.05 after making the suggested changes above.
chs_x = degs_temp %>%
  slice_max( abs(log2FC), n = 500 ) 

heat_matrix = exprset1[ chs_x$probe_id, ]

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

### 3.4 heatmap----

chs_x = degs_temp %>%
  slice_max( abs(log2FC), n = 500 ) 

heat_matrix = exprset1[ chs_x$probe_id, ]

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

## 4、probe annotation----

# This should be done before the analyses above.
# I'm not familiar with AnnoProbe so perhaps the code below accounts for this already,
# but I would suggest annotating with Entrez (NCBI) gene IDs as well as gene symbols.
# This would help with compatibility with packages like clusterProfiler.

library(AnnoProbe)
if (!requireNamespace("hgu133plus2.db", quietly = TRUE))
  BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL)
gpl_anno = ids %>% 
  as.data.frame() %>%
  distinct( probe_id, .keep_all = T ) 
rownames(gpl_anno) = gpl_anno$probe_id
head( gpl_anno )

exprset2 = exprset1[rownames(exprset1) %in% gpl_anno$probe_id,] 
dim(exprset2)
gpl_anno = gpl_anno[gpl_anno$probe_id %in% rownames(exprset2),] 

### 4.2 expression annotation----
identical( gpl_anno$probe_id, rownames(exprset2))
gpl_anno = gpl_anno[rownames(exprset2),] 

identical( gpl_anno$probe_id, rownames(exprset2))
exprset2 = cbind(gpl_anno,exprset2)



###  DEG annotation without deduplication----
degs_anno = degs_temp %>% 
  filter( rownames(.) %in% gpl_anno$probe_id )
gpl_anno2 = gpl_anno[ degs_anno$probe_id, ] 

identical( gpl_anno2$probe_id,rownames(degs_anno))
degs_anno = cbind( gpl_anno2,degs_anno) %>%
  dplyr::select(-1)


### 4.4 deduplication for DEG----
summary( degs_anno$Mean.Expr)
plot(table(sort(table(degs_anno$symbol))))
table( is.na(degs_anno$P.value) )
table( degs_anno$P.value == 0 )

# deduplication for genes (the expression larger will be preserved)
degs_final = degs_anno %>% 
  group_by( symbol ) %>%
  mutate( mean = mean(Mean.Expr) ) %>% 
  group_by( ) %>% # 
  filter( Mean.Expr >= mean | Mean.Expr >= mean(Mean.Expr) ) %>%
  arrange( P.value ) %>%
  distinct( symbol, .keep_all = T ) 

# deduplication for DEG matrix 
exprset3 = exprset1[ degs_final$probe_id, ] 
identical( degs_final$probe_id, rownames(exprset3))
#
rownames(exprset3) = degs_final$symbol

rm(pos, exprset1, exprset2, degs_temp, degs_anno, gpl_anno, gpl_anno2, xxx )


## 4.5 filter via Mean.Expr and log2FC----
summary(degs_final$Mean.Expr)
summary(abs(degs_final$log2FC))

# set the threshold value
qu1 = quantile(degs_final$Mean.Expr, probs = 0.25 ) 
fc1 = mean(abs(degs_final$log2FC))  
qu1;fc1

### DEG filter +ID change----
library(hgu133plus2.db)
library(clusterProfiler) 
library(dplyr)
#
degs_top = degs_final %>%
  filter( Mean.Expr > qu1 & abs(log2FC) > fc1 & P.value < 0.05 ) %>%  
  arrange( P.value ) %>% 
  arrange( desc( abs(log2FC) ) ) %>%
  pull( symbol ) %>% 
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), 
        OrgDb = hgu133plus2.db ) %>%  
  slice_head( n = 300 ) 

## e1.3 enrich based on TOP300----
library(clusterProfiler) 
library(enrichplot) 

load( paste('4.5_Top300DEG_','.Rdata',sep = "") )

### e1.3 GOenrich----
go_data <- enrichGO( gene = degs_top$ENTREZID, 
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
### e1.4 KEGG enrich----
kegg_data <- enrichKEGG(gene = degs_top$ENTREZID, 
                        keyType = "kegg",
                        organism   = 'hsa', # mmu is mouse，hsa is human
                        pvalueCutoff = 1, # 1 means no selected
                        qvalueCutoff = 1,
                        use_internal_data = F 
) %>%
  setReadable( ., OrgDb = 'hgu133plus2.db', keyType = 'ENTREZID' )  %>% # ENTRZID change to SYMBOL
  pairwise_termsim() 
write.csv(kegg_data, file = 'KEGG-richment.csv')


### e1.5 visualization for enrich----
library(DOSE)
enrichDO()
enrichWP()
library(ggplot2)
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

