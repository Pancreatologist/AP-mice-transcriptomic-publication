### Gene ID change (from "SYMBOL" to "ENTREZID") to further using 'clusterProfiler'
gpl_anno=read.table('GPL10787.txt', sep = "\t", quote = "\"",header = TRUE, fill = TRUE)  # read the annotation file for the GPL10787 platform (downloaded)
head(gpl_anno)
#gpl_anno <- gpl_anno[,c(1,6,7,9)] if it need to filter the col
gpl_anno = gpl_anno %>% distinct(GENE_SYMBOL, .keep_all = TRUE ) # delete the repeat probe
rownames(gpl_anno) = gpl_anno$GENE_SYMBOL

### annotation for eset
eset_anno = eset[rownames(eset) %in% gpl_anno$GENE_SYMBOL,] #  delete the unmatched in eset
dim(eset_anno)
gpl_anno = gpl_anno[gpl_anno$GENE_SYMBOL %in% rownames(eset_anno),] #  delete the unmatched in annotation
identical( gpl_anno$GENE_SYMBOL, rownames(eset_anno)) #check the arrange of eset and annotation
gpl_anno = gpl_anno[rownames(eset_anno),] #  arrange 
identical( gpl_anno$GENE_SYMBOL, rownames(eset_anno)) #recheck
eset2 = cbind(gpl_anno,eset_anno)

### annotation for DEG and change into ENTREZID for the clusterProfiler
degs_anno = degs_temp %>% 
  filter( rownames(.) %in% gpl_anno$GENE_SYMBOL )%>%
  mutate(GENE_SYMBOL= rownames(.))
gpl_anno2 = gpl_anno[ degs_anno$GENE_SYMBOL, ] # arrange
identical( gpl_anno2$GENE_SYMBOL,rownames(degs_anno))
gpl_anno2 <- gpl_anno2[,c(1,2,4)] 
degs_anno = cbind( gpl_anno2,degs_anno) 
summary( degs_anno$AveExpr)
table( is.na(degs_anno$adj.P.Val) )
table( degs_anno$adj.P.Val == 0 )

degs_top  = degs_anno %>%
  filter(adj.P.Val < 0.05 ) %>%
  #filter( AveExpr > qu1 & abs(log2FC) > fc1) %>%  # filter the expr more than 25% and log2FC more than mean with adj.P.Val less than 0.05 if needed
  arrange( adj.P.Val ) %>% # arrange via adj.P value
  arrange( desc( abs(logFC) ) ) %>%
  pull( GENE_SYMBOL ) %>% # choose the col of 'SYMBOL' and make it to a vector
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), 
        OrgDb = org.Mm.eg.db ) #%>%  # from 'symbol' to 'ENTREZID'
  #slice_head( n = 300 ) # slice for the top 300 if need
  
### Enrichment
### Enrichment Based the TOP DEGs
### GO pathway
go_data <- enrichGO( gene = degs_top$ENTREZID, 
                     OrgDb="org.Mm.eg.db", # Mm is mouse，human is Hs
                     ont ="ALL", # ALL includes BP、CC、MF
                     pvalueCutoff = 0.05, # make P-value to 0.05
                     qvalueCutoff = 1,
                     readable= TRUE # change ENTREZID to SYMBOL in the result 
) %>% 
  clusterProfiler::simplify( . ) %>% # delete the repeat or highly similar col
  pairwise_termsim() # for the further emap 
### GO pathway visualization
emapplot(go_data, showCategory = 20) +
  ggtitle("Top enriched terms")+
  theme_light()+
  theme(axis.text = element_blank(), axis.title = element_blank())+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

cnetplot(go_data,showCategory = 10, 
         circular=F, 
         colorEdge= F)+
  ggtitle("Terms and DEGs") +
  theme_light()+
  theme(axis.text = element_blank(), axis.title = element_blank())
### KEGG pathway
kegg_data <- enrichKEGG(gene = degs_top$ENTREZID, 
                        keyType = "kegg",
                        organism   = 'mmu', # mmu is mouse，human is hsa
                        pvalueCutoff = 0.05, # filter via pvalue
                        qvalueCutoff = 1,
                        use_internal_data = F 
) %>%
  setReadable( ., OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID' )  %>% # from ENTRZID to SYMBOL in result
  pairwise_termsim() # emap need this transform
### KEGG pathway visualization
dotplot(kegg_data, showCategory= 20, title="Top enriched terms") + 
  theme_light()



