### pre-process of expression matrix for ImmuCC AI
eset = exprs.MA(MA.reomoveCTprobe)#get the expression from MA
colnames(eset) <- rownames(targets2)
rownames(eset) <- MA.reomoveCTprobe$genes$GeneName 
eset <- removeBatchEffect(eset, as.factor(targets2$Batch))
esetwithoutlog <- exp(eset)

### import the results and visualization
cell_prop <-read.csv("ImmuCC.txt",sep = '') # ImmuCC
cell_prop <- cell_prop %>%  dplyr::select(-Infiltration_score)
Group <- rownames(cell_prop)
for (i in seq_along(Group)) {
  row_name_without_suffix <- strsplit(Group[i], "_")[[1]][1]
  Group[i] <- row_name_without_suffix
}
cell_prop <- cbind(cell_prop,Group)
cell_prop <- cell_prop[,c("B_cell","Dendritic_cells","Granulocytes","Macrophage","Monocytes","NK","T_cell",'Group')]

heatmatrix_long <- cell_prop
heatmatrix_long <- heatmatrix_long %>%
  group_by(Group) %>%
  summarize(across(everything(), mean))
heatmatrix_long <- heatmatrix_long %>%
  gather(key = "variable", value = "value", -Group)
heatmatrix_long <- as_tibble(heatmatrix_long)
heatmatrix_long$value <- as.numeric(heatmatrix_long$value)
heatmatrix_long
heatmatrix_long <- heatmatrix_long[heatmatrix_long$Group%in% c('Control','C7plus6',"C12","FAEE50",'FAEE150','TLCS3','TLCS5'),]
heatmatrix_long_groupings = 
  heatmatrix_long |>
  mutate(Variable_group = if_else(Group %in% c("C12", 'C7plus6'), "CER",
                                  if_else(Group %in% c("FAEE50",'FAEE150'), "FAEE", #
                                         if_else(Group %in% c('TLCS3','TLCS5'),'TLCS','Control')
                                  )))

heatmatrix_long_groupings$Variable_group <- factor(heatmatrix_long_groupings$Variable_group, levels = c('Control',"CER", "FAEE",'TLCS'))
heatmatrix_long_groupings$Group <- factor(heatmatrix_long_groupings$Group, levels = c('Control','C7plus6',"C12","FAEE50",'FAEE150','TLCS3','TLCS5'))
heatmatrix_long_groupings$variable <- factor(heatmatrix_long_groupings$variable,levels = c("Granulocytes","Macrophage","Monocytes","B_cell","Dendritic_cells","NK","T_cell",'Group'))

heatmatrix_long_groupings
heatmatrix_long_groupings|> 
  group_by(Variable_group) |>
  heatmap(variable,Group, value,na_col = "white",rect_gp = gpar(col="black"),
          palette_value = circlize::colorRamp2(c(-2,0,2), c('#8D91AA','white','#B98A82' )),
          scale = 'row',
          cluster_col = F,cluster_rows = F,
          palette_grouping = list(
            # For first grouping (vs)
            c("#BFCAC2", "#013E41"), 
            # For second grouping (property_group)
            c("#C6DEE0","#F7EDEB")
          )
  ) |> save_pdf("ImmuCell.pdf",width =7,height =4,units = c("in"))
