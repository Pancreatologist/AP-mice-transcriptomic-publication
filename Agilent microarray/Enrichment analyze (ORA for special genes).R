#ORA for GO-BP
gene = degs_temp[degs_temp$Condition=='Up',] %>% 
  pull( GeneName) %>% # choose the col of 'SYMBOL' and make it to a vector
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), 
        OrgDb = org.Mm.eg.db ) #%>%  # from 'symbol' to 'ENTREZID'
gene = degs_temp[degs_temp$Condition=='Down',] %>% 
  pull( GeneName) %>% # choose the col of 'SYMBOL' and make it to a vector
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), 
        OrgDb = org.Mm.eg.db ) #%>%  # from 'symbol' to 'ENTREZID'
geneList = degs_temp$GeneName %>% 
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"), 
        OrgDb = org.Mm.eg.db ) #%>%  # from 'symbol' to 'ENTREZID'

ego <- enrichGO(gene          = gene$ENTREZID,
                universe      = geneList$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)%>% 
  clusterProfiler::simplify( . ) %>% # delete the repeat or highly similar col
  pairwise_termsim() # for the further emap 



#Chord plot
enrichgo <- DOSE::setReadable(ego, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
GO <- enrichgo@result[1:10,c('Description','qvalue','geneID')] #choose the top 10 pathways
tmp=do.call(rbind,
            apply(GO, 1,function(x){
              data.frame(go=x[1],
                         gene=strsplit(x[3],'/')[[1]])
            })
)
tmp2=dcast(tmp,go~gene)
tmp2[is.na(tmp2)]=0
rownames(tmp2)=tmp2[,1]
tmp2=tmp2[,-1]
tmp2=t(tmp2)
tmp2[tmp2!=0]=1
tmp2=as.data.frame(tmp2)
table(rownames(degs_temp[rownames(tmp2),]) == rownames(tmp2))
cg=rownames(tmp2)
tmp2=apply(tmp2,2,as.numeric)
rownames(tmp2)=cg
cor_mat <- cor(tmp2) 
col_mat <- rand_color(nrow(cor_mat), transparency = 0.3) 
#col_fun <- colorRamp2(c(-1, 0, 1),c("#42d4f4", "white", "#e6194B"))
chordDiagram(cor_mat, link.sort = TRUE, link.decreasing = TRUE, grid.col = col_mat,#directional = 1,
             transparency = 0.3,
             annotationTrack = c("grid"),
             symmetric = TRUE)
circos.clear()



