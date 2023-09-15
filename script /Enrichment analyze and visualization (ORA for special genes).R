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


### visualization
ego <- ego %>% as.data.frame() %>% 
  mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
ego$Description
ggplot(ego[1:5,], #c(4,7,10,9,12)
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
    scale_fill_manual(values = as.vector(palette.colors()))+  
  #scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )+
  scale_x_continuous(expand = expansion(c(0, 0.1)))+
  background_grid(major = "y", minor = "y")+
  theme_cowplot(20) + 
  xlab("rich factor") + ylab("Terms") +
  ggtitle("ORA for GO")


