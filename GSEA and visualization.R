###gseaGO
alldiff <- degs_temp[order(degs_temp$logFC,decreasing = T),]
genelist <- alldiff$logFC
names(genelist) <- alldiff$ENTREZID
ego2 <- gseGO(geneList = genelist,
              ont = 'BP',
              OrgDb = org.Mm.eg.db,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              minGSSize = 15,
              maxGSSize = 500,
              seed = FALSE,
              by = "fgsea",
              keyType = "ENTREZID") %>% 
  DOSE::setReadable(OrgDb='org.Mm.eg.db',keyType='ENTREZID')
ego2 <- simplify(ego2,  cutoff = 0.7,
                 by = "p.adjust",
                 select_fun = min,
                 measure = "Wang",
                 semData = NULL)
ego2 <- ego2 %>% as.data.frame()
ego2$Description


###chord plot
enrichgo <- DOSE::setReadable(ego2, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
#GO <- enrichgo@result[1:10,c('Description','qvalue','geneID')] 
GO <- enrichgo@result[1:10,c('Description','qvalue','core_enrichment')] 
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
col_mat <- rand_color(10, transparency = 0.4) 
#col_fun <- colorRamp2(c(-1, 0, 1),c("#42d4f4", "white", "#e6194B"))
chordDiagram(cor_mat, link.sort = TRUE, link.decreasing = TRUE, grid.col = col_mat,#directional = 1,
             transparency = 0.3,
             annotationTrack = c("name","grid"),
             symmetric = TRUE)
circos.clear()

### GSEA plot
ego2@result[grepl('B cell|MHC|tumor necrosis factor production|neutrophil migration|myeloid leukocyte activation',ego2@result$Description),c('NES','ID','Description','p.adjust')]
terms <- c('GO:0019724',
           'GO:1990266',
           'GO:0032680',
           'GO:0002274',
           'GO:0019886')
ego2@result[ego2@result$ID %in% terms,c('NES','ID','Description','p.adjust')]
gseaplot2(ego2, geneSetID = terms, pvalue_table = FALSE,rel_heights = c(1.5, 0.5, 0.5),base_size = 15,subplots=1:2,
          color = c("#4A475C", "#ACD4D6", "#A6BAAF",'#B98A82','#E4DBD2'), ES_geom = "line")

### summarized results via heatmap
innate <- c('neutrophil migration','neutrophil activation','myeloid leukocyte activation','leukocyte migration','leukocyte cell-cell adhesion',
'leukocyte aggregation','myeloid leukocyte cytokine production','regulation of leukocyte proliferation','neutrophil chemotaxis','granulocyte activation',
'granulocyte chemotaxis','granulocyte migration','tumor necrosis factor superfamily cytokine production','tumor necrosis factor production')
adaptive <- c('antigen processing and presentation of exogenous peptide antigen via MHC class II','antigen processing and presentation of peptide antigen via MHC class II',
              'B cell mediated immunity','T cell receptor signaling pathway','T cell mediated immunity','regulation of B cell activation',
              'B cell receptor signaling pathway','B cell activation involved in immune response','regulation of T cell activation',
              'B cell activation','B cell differentiation')
all <- c(innate,adaptive)
filtered_rows <- ego2[ego2$Description %in% all,c('Description','NES')]
write.csv(filtered_rows, file = "filtered_data.csv", row.names = FALSE)
#make the summarized result from the each filtered_data.csv, and put them into the GSEAheatmap.txt
heatmatrix_long <- read.delim2('GSEAheatmap.txt')#read from the summarized result
heatmatrix_long <- as_tibble(heatmatrix_long)
heatmatrix_long$NES <- as.numeric(heatmatrix_long$NES)
heatmatrix_long
innate <- c('neutrophil migration','neutrophil activation','myeloid leukocyte activation','leukocyte migration','leukocyte cell-cell adhesion',
            'leukocyte aggregation','myeloid leukocyte cytokine production','regulation of leukocyte proliferation','neutrophil chemotaxis','granulocyte activation',
            'granulocyte chemotaxis','granulocyte migration','tumor necrosis factor superfamily cytokine production','tumor necrosis factor production')
adaptive <- c('antigen processing and presentation of exogenous peptide antigen via MHC class II','antigen processing and presentation of peptide antigen via MHC class II',
              'B cell mediated immunity','T cell receptor signaling pathway','T cell mediated immunity','regulation of B cell activation',
              'B cell receptor signaling pathway','B cell activation involved in immune response','regulation of T cell activation',
              'B cell activation','B cell differentiation')
heatmatrix_long_groupings = 
  heatmatrix_long |>
  mutate(Variable_group = if_else(Group %in% c("C7plus6", "TLCS3",'FAEE150'), "M",'S'))|>
  mutate(immune = if_else(Description %in% innate, "innate",'adaptive'))
heatmatrix_long_groupings$Variable_group <- factor(heatmatrix_long_groupings$Variable_group, levels = c( "M",'S'))
heatmatrix_long_groupings$Description <- factor(heatmatrix_long_groupings$Description,levels = c(innate,adaptive))
heatmatrix_long_groupings$Group <- factor(heatmatrix_long_groupings$Group, levels = c('C7plus6','TLCS3','FAEE150', 'CER','TLCS','FAEE'
))
heatmatrix_long_groupings
heatmatrix_long_groupings|> 
  group_by(Variable_group,immune) |>
  heatmap(Description,Group, NES,na_col = "#F5F4F0",rect_gp = gpar(col="black"),
          palette_value = circlize::colorRamp2(c(-2,0,2), c("#5B7493", "#BFBFC1","#E7ADAC")),
          cluster_col = FALSE,cluster_rows = F,
          palette_grouping = list(
            # For first grouping (vs)
            c("#66C2A5", "#FC8D62"), 
            # For second grouping (property_group)
            c("#b58b4c","#74a6aa")
          )
  ) |> 
  save_pdf("GSEA_heatmap.pdf")



