###### Figure 1 #######

merged_tpm_matrix<- read.table("HUMAN_MODELS_TPM.txt")
merged_annotation_df<- read.table("HUMAN_MODELS_METADATA.txt")

merged_tpm_matrix<- unlog_matrix(merged_tpm_matrix)

nm_merged_annotation_df<- merged_annotation_df[merged_annotation_df$CNA!="Malignant",]
cell_prop <- data.frame(table(nm_merged_annotation_df[,c("Organoid_line","Cell_type")])/rowSums(table(nm_merged_annotation_df[,c("Organoid_line","Cell_type")])))

cols<-pal_d3("category10")(7)
names(cols)<- c("aRG"  , "oRG.Astroglia" ,"Cycling.oRG.Astroglia",  "Immature.interneurons", "IP", "CPN" ,"CFuPN")  

pdf("Barplot_cell_types.pdf",width = 4, height = 6,onefile=F)
ggbarplot(cell_prop, x = "Organoid_line",y="Freq",fill = "Cell_type",palette = cols)
dev.off()

# CNA_plots

all_ell_groups<- split(rownames(merged_annotation_df),merged_annotation_df$Cell_type)

ref_cell_groups<- all_ell_groups[names(all_ell_groups)!="Malignant"]

log_expression_matrix<- log2(1+merged_tpm_matrix/10)

top.genes<- get_top_genes(merged_tpm_matrix)


cna = infercna::infercna(m = log_expression_matrix[top.genes,unlist(all_ell_groups)],refCells = ref_cell_groups, n = 10000, noise = 0.1, isLog = TRUE, verbose = FALSE)


gfp_order_list<-list("GFP+_nonmalignant"= rownames(merged_annotation_df)[merged_annotation_df$CNA=="Non-malignant"& merged_annotation_df$Type=="GFP_POS"],"GFP-_nonmalignant"= rownames(merged_annotation_df)[merged_annotation_df$CNA=="Non-malignant"& merged_annotation_df$Type=="GFP_NEG"], "Malignant"= rownames(merged_annotation_df)[merged_annotation_df$CNA=="Malignant"])


gfp_order_list$`GFP+_nonmalignant`<- sample(gfp_order_list$`GFP+_nonmalignant`, 500)
gfp_order_list$`GFP-_nonmalignant`<- sample(gfp_order_list$`GFP-_nonmalignant`, 500)

lapply(gfp_order_list, length)


cna_plot<- infercna::cnaPlot(cna = cna[,unlist(gfp_order_list)], order.cells = gfp_order_list)

png("Final_paper_plots/CNA_merged.png",width = 9, height = 6,units = 'in', res = 300,type="cairo")
print(cna_plot$p)
dev.off()


###### Figure 2 #########

### malignant cell states

NMF_programs<- read.table("TableS3_DeNovoPrograms.txt")
NMF_programs<-as.list(NMF_programs)[-1]

##inter tumor

malignant_expression_matrix<- merged_tpm_matrix[,rownames(merged_annotation_df)[merged_annotation_df$CNA=="Malignant"]]
log_expression_matrix<- data.frame(log2(1+malignant_expression_matrix/10))
colnames(log_expression_matrix)<- str_replace(colnames(log_expression_matrix),"[.]","-")

##intra-tumor

Signatures_GBM<- scalop::Signatures_GBM

intra_tumor_log_expression_matrix<-data.frame(row.names = rownames(log_expression_matrix))
score_df<-sigScores(log_expression_matrix, Signatures_GBM)
score_df_NMF<-sigScores(log_expression_matrix, NMF_programs)


for(tumor in unique(merged_annotation_df$Sample))
{
  temp_matrix<-  merged_tpm_matrix[,rownames(merged_annotation_df)[merged_annotation_df$CNA=="Malignant" & merged_annotation_df$Sample==tumor]]
  
  
  temp_matrix<- data.frame(temp_matrix)
  colnames(temp_matrix)<- str_replace(colnames(temp_matrix),"[.]","-")
  temp_score_df<- sigScores(temp_matrix, Signatures_GBM)
  score_df[rownames(temp_score_df),]<- temp_score_df
  temp_score_df<- sigScores(temp_matrix, NMF_programs)
  score_df_NMF[rownames(temp_score_df),]<- temp_score_df
  
  temp_matrix<- temp_matrix-rowMeans(temp_matrix)
  intra_tumor_log_expression_matrix<- cbind(intra_tumor_log_expression_matrix, temp_matrix)
}



NMF_heatmap_all<- NMF_heatmap(NMF_programs = NMF_programs,log_expression_matrix = intra_tumor_log_expression_matrix,showAllGenes1 = F, center = T,num_genes = 20 )

NMF_heatmap_all$heatmap


p1<-basic_heatmap_wrapper(intra_tumor_log_expression_matrix[unlist(NMF_programs),NMF_heatmap_all$cells_order[merged_annotation_df[NMF_heatmap_all$cells_order,]$Sample=="MGG123"]])+ggtitle("MGG123")+NoLegend()

p2<-basic_heatmap_wrapper(intra_tumor_log_expression_matrix[unlist(NMF_programs),NMF_heatmap_all$cells_order[merged_annotation_df[NMF_heatmap_all$cells_order,]$Sample=="MGH143"]])+ggtitle("MGH143")+NoLegend()


p3<-basic_heatmap_wrapper(intra_tumor_log_expression_matrix[unlist(NMF_programs),NMF_heatmap_all$cells_order[merged_annotation_df[NMF_heatmap_all$cells_order,]$Sample=="MGG23"]])+ggtitle("MGG23")+NoLegend()


p4<-basic_heatmap_wrapper(intra_tumor_log_expression_matrix[unlist(NMF_programs),NMF_heatmap_all$cells_order[merged_annotation_df[NMF_heatmap_all$cells_order,]$Sample=="MGG101"]])+ggtitle("MGG101")+NoLegend()


p5<-basic_heatmap_wrapper(intra_tumor_log_expression_matrix[unlist(NMF_programs),NMF_heatmap_all$cells_order[merged_annotation_df[NMF_heatmap_all$cells_order,]$Sample=="MGG125"]])+ggtitle("MGG125")+NoLegend()


p6<- basic_heatmap_wrapper(intra_tumor_log_expression_matrix[unlist(NMF_programs),NMF_heatmap_all$cells_order[merged_annotation_df[NMF_heatmap_all$cells_order,]$Sample=="MGG65"]])+ggtitle("MGG65")+NoLegend()


p7<-basic_heatmap_wrapper(intra_tumor_log_expression_matrix[unlist(NMF_programs),NMF_heatmap_all$cells_order[merged_annotation_df[NMF_heatmap_all$cells_order,]$Sample=="MGG70"]])+ggtitle("MGG70")+NoLegend()


p8<-basic_heatmap_wrapper(intra_tumor_log_expression_matrix[unlist(NMF_programs),NMF_heatmap_all$cells_order[merged_annotation_df[NMF_heatmap_all$cells_order,]$Sample=="MGG75"]])+ggtitle("MGG75")+NoLegend()


p9<-basic_heatmap_wrapper(intra_tumor_log_expression_matrix[unlist(NMF_programs),NMF_heatmap_all$cells_order[merged_annotation_df[NMF_heatmap_all$cells_order,]$Sample=="MGG87"]])+ggtitle("MGG87")+NoLegend()

merged_annotation_df$cells<-rownames(merged_annotation_df)

annotation_sample <- ggplot(merged_annotation_df[NMF_heatmap_all$cells_order,], aes(x=factor(cells,levels = cells), y = 1)) + geom_raster(aes(fill=Sample))+scale_fill_manual(values = c("#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",   "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF"), drop = FALSE) + theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text=element_blank(), axis.line = element_blank(), axis.text.x = element_blank(),legend.position = "top")+ labs(x="",y="",fill="") 

annotation_Organoid <- ggplot(merged_annotation_df[NMF_heatmap_all$cells_order,], aes(x=factor(cells,levels = cells), y = 1)) + geom_raster(aes(fill=Organoid_line)) +scale_fill_manual(values = c("#3B4992FF","#EE0000FF","#008B45FF"), drop = FALSE)+ theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text=element_blank(), axis.line = element_blank(), axis.text.x = element_blank())+ labs(x="",y="",fill="")


score_plot<- basic_heatmap_wrapper(t(score_df[NMF_heatmap_all$cells_order,c("G1S","G2M","OPC","NPC1","NPC2","AC","MES1","MES2")]),showAllGenes = T, limits = 1.5) + theme(axis.text.x = element_blank(), axis.title=element_blank())

nmf_heatmap_plot<-NMF_heatmap_all$heatmap


pdf("malignant_cells_heatmap.pdf",width = 10, height = 10,onefile=F)
egg::ggarrange(annotation_sample,annotation_Organoid,nmf_heatmap_plot,score_plot, nrow=4, heights = c(1,1,50,8))
dev.off()

pdf("malignant_cells_heatmap_pertumor.pdf",width = 10, height = 10,onefile=F)
egg::ggarrange(p3,p2,p5,p6,p7,p8,p9,p4,p1)
dev.off()


vector<- maxcol_strict(as_four_state_gbm(score_df[,c("MES1","MES2","NPC1","NPC2","OPC","AC")]),min=0.2)
state_list<-split(vector,factor(names(vector)))

state_list$unclassified<- rownames(score_df)[!rownames(score_df)%in% vector]

group_cols<- c("#0bb5fd","#fce232","#ff80a0","#63bf2a","grey")
names(group_cols)<- c("NPC","AC","MES","OPC","unclassified")

pie_chart_list<-list()

for(cell_line in unique(merged_annotation_df$Sample))
{
  temp_state_list<- lapply(state_list,function(x){x[x %in% rownames(merged_annotation_df)[merged_annotation_df$Sample==cell_line]]})
  
  pie_chart_list[[cell_line]]<-create_state_piechart(temp_state_list)+ scale_fill_manual(values=group_cols)+ggtitle(cell_line)
  
}


pdf("canonical_pie_pertumor.pdf",width = 10, height = 15,onefile=F)
pie_chart_list$MGH143+NoLegend()+pie_chart_list$MGG23+NoLegend()+pie_chart_list$MGG65+NoLegend()+pie_chart_list$MGG75+NoLegend()+pie_chart_list$MGG87+NoLegend()+pie_chart_list$MGG70+NoLegend()+pie_chart_list$MGG123+NoLegend()+pie_chart_list$MGG125+NoLegend()+pie_chart_list$MGG101
dev.off()



vector<- maxcol_strict((score_df_NMF),min=0.2)
state_list<-split(vector,factor(names(vector)))

state_list$unclassified<- rownames(score_df_NMF)[!rownames(score_df_NMF)%in% vector]


group_cols<- c("#0bb5fd","#fce232","#ff80a0","#63bf2a")
names(group_cols)<- c("NPC","AC","MES","OPC")


group_cols<- c("#0bb5fd","#fce232","#ff80a0","pink","#63bf2a","#E3F9A6","grey")
names(group_cols)<- c("NPC","AC","MES1","MES2","OPC","CC","unclassified")

pie_chart_list<-list()

for(cell_line in unique(merged_annotation_df$Sample))
{
  temp_state_list<- lapply(state_list,function(x){x[x %in% rownames(merged_annotation_df)[merged_annotation_df$Sample==cell_line]]})
  
  pie_chart_list[[cell_line]]<-create_state_piechart(temp_state_list)+ scale_fill_manual(values=group_cols)+ggtitle(cell_line)
  
}


pdf("NMFl_pie_pertumor.pdf",width = 10, height = 15,onefile=F)
pie_chart_list$MGH143+NoLegend()+pie_chart_list$MGG23+NoLegend()+pie_chart_list$MGG65+NoLegend()+pie_chart_list$MGG75+NoLegend()+pie_chart_list$MGG87+NoLegend()+pie_chart_list$MGG70+NoLegend()+pie_chart_list$MGG123+NoLegend()+pie_chart_list$MGG125+NoLegend()+pie_chart_list$MGG101
dev.off()


## comparing in vitro vs in organoid

Signatures_GBM<-Signatures_GBM[c("AC","MES1","MES2","OPC","NPC1","NPC2","G1S","G2M")] 



####### Figure 3 ########

merged_annotation_df$CNA<- factor(merged_annotation_df$CNA,levels=c("Malignant","Non-malignant"))

merged_annotation_df$Type<- factor(merged_annotation_df$Type,levels=c("GFP_POS","GFP_NEG"))

merged_annotation_df$Cell_type<-factor(merged_annotation_df$Cell_type, levels=c("aRG"  , "oRG.Astroglia" ,"Cycling.oRG.Astroglia",  "Immature.interneurons", "IP", "CPN","CFuPN"))

pdf("Final_paper_plots/GFPboxplot.pdf",width = 6.5, height = 5,onefile=F)
ggboxplot(merged_annotation_df, x = "CNA", y = "log_counts", color = "black",facet.by = "Type",
          add = "jitter", add.params = list(size=0.1))+ ylab("log (# detected counts of GFP transcript)")+xlab("Cells")+ stat_compare_means(aes(group = CNA), label = "p.signif")
dev.off()


pdf("Final_paper_plots/GFPboxplot_nomalignant.pdf",width = 7, height = 5,onefile=F)
ggboxplot(merged_annotation_df[merged_annotation_df$CNA!="Malignant",], x = "Cell_type", y = "log_counts", color = "Type",palette =c("chartreuse","darkgrey"),add = "jitter", add.params = list(size=0.1))+ rotate_x_text(45)+ ylab("log (# detected counts of GFP transcript)")+xlab("Cells")+ stat_compare_means(aes(group = Type), label = "p.signif")
dev.off()



sum(merged_annotation_df[merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_POS",]$log_counts>1) ## GFPpos nonmalignant cells with detected GFP

nrow(merged_annotation_df[merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_POS",]) ## total GFPpos nonmalignant cells


sum(merged_annotation_df[merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_NEG",]$log_counts>1) ## GFPneg nonmalignant cells with detected GFP

nrow(merged_annotation_df[merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_NEG",]) ## total GFPneg nonmalignant cells

## removing malignant cells + td tomato samples
merged_annotation_df_nm<-merged_annotation_df[!merged_annotation_df$Sample %in% c("MGG87","MGG123"),]
merged_annotation_df_nm<- merged_annotation_df_nm[merged_annotation_df_nm$CNA!="Malignant",]


prop_df<- data.frame(row.names = c("aRG",  "oRG.Astroglia" ,"Cycling.oRG.Astroglia",  "Immature.interneurons", "IP", "CPN" ,"CFuPN"))


for(Cell_type in rownames(prop_df))
{
  temp_df<- merged_annotation_df_nm[merged_annotation_df_nm$Cell_type==Cell_type,]

  prop_df[Cell_type,"GFP_pos"]<- sum(temp_df$Cell_type==Cell_type & temp_df$Type=="GFP_POS")/sum(temp_df$Cell_type==Cell_type)
  prop_df[Cell_type,"GFP_neg"]<- sum(temp_df$Cell_type==Cell_type & temp_df$Type=="GFP_NEG")/sum(temp_df$Cell_type==Cell_type)
  
  group2<- sum(merged_annotation_df_nm$Type=="GFP_POS") #GFP+ cells
  group1<- sum(merged_annotation_df_nm$Cell_type==Cell_type) # cell type cells (oRG, IP,etc)
  overlap<- sum(merged_annotation_df_nm$Type=="GFP_POS"&merged_annotation_df_nm$Cell_type==Cell_type) # both
  total<- nrow(merged_annotation_df_nm) # all non malig cells
  
  print(Cell_type)
  print(phyper(q=overlap, m=group2, n = total-group2, k=group1,lower.tail= FALSE))

}

pdf("GFP_nm_prop_barplot.pdf",width = 6.5, height = 5,onefile=F)
ggbarplot(melt(as.matrix(prop_df)),y="value",x="Var1",fill="Var2", palette = c("chartreuse","lightgrey"),add="mean")+ rotate_x_text(45)+ theme(axis.title=element_blank())
dev.off()

# per sample
prop_df_sample<-list()
plot_df_sample<-list()


for(sample in unique(merged_annotation_df$Sample))
{
  temp_prop_df<- data.frame(row.names = levels(merged_annotation_df$Cell_type))
  temp_prop_df<- temp_prop_df[rownames(temp_prop_df)%in% as.character(merged_annotation_df[merged_annotation_df$Sample==sample,]$Cell_type),]
  
  for(Cell_type in rownames(temp_prop_df))
  {
    temp_df<- merged_annotation_df[merged_annotation_df$Cell_type==Cell_type & merged_annotation_df$Sample==sample,]
    temp_df<- temp_df[complete.cases(temp_df),]
    
    temp_prop_df[Cell_type,"GFP_pos"]<- sum(temp_df$Cell_type==Cell_type & temp_df$Type=="GFP_POS")/sum(temp_df$Cell_type==Cell_type)
    temp_prop_df[Cell_type,"GFP_neg"]<- sum(temp_df$Cell_type==Cell_type & temp_df$Type=="GFP_NEG")/sum(temp_df$Cell_type==Cell_type)
    
  }
  temp_prop_df<- temp_prop_df[complete.cases(temp_prop_df),]
  
  
  
  melt_pie<- melt(as.matrix(temp_prop_df))
  prop_df_sample[[sample]]<-  melt_pie
  
  melt_pie$Var1<- factor(melt_pie$Var1, levels=c("aRG",  "oRG.Astroglia" ,"Cycling.oRG.Astroglia",  "Immature.interneurons", "IP", "CPN" ,"CFuPN"))
  
  
  if(sample %in% c("MGG87","MGG123")) #tdtomato 
    plot_df_sample[[sample]]<-   print(ggbarplot(melt_pie,y="value",x="Var1",fill="Var2", palette = c("red","lightgrey"))+ ggtitle(sample)+ rotate_x_text(45)+ theme(axis.title=element_blank())+NoLegend())
  else # gfp
    plot_df_sample[[sample]]<-   print(ggbarplot(melt_pie,y="value",x="Var1",fill="Var2", palette = c("chartreuse","lightgrey"))+ ggtitle(sample)+ rotate_x_text(45)+ theme(axis.title=element_blank())+NoLegend())
  
  
  
  
}

# sample without replacement random (size of celltype) cells

for(Cell_type in unique(merged_annotation_df_nm$Cell_type))
{
  size= sum(merged_annotation_df_nm$Cell_type==Cell_type)
  prop<- sum(merged_annotation_df_nm$Cell_type==Cell_type&merged_annotation_df_nm$Type=="GFP_POS" )/(size)
  vector_prop<-c()
  for(i in 1:1000)
  {
    sampled_cells<- sample(rownames(merged_annotation_df_nm),size)
    vector_prop[i]<- sum(merged_annotation_df_nm[sampled_cells,]$Type=="GFP_POS" )/length(sampled_cells)
  }
  
  prop_df[Cell_type,"pval"]<-abs(1000-sum(prop>(0.05+vector_prop)))/1000
  
}


pdf("Final_paper_plots/S3_GFPcelltypebarplot_all.pdf",width = 15, height = 15,onefile=F)
egg::ggarrange(plot_df_sample$MGG101,plot_df_sample$MGG125,plot_df_sample$MGG23,plot_df_sample$MGG65,plot_df_sample$MGG70,plot_df_sample$MGG75,plot_df_sample$MGH143,plot_df_sample$MGG87,plot_df_sample$MGG123)
dev.off()


######### Figure 4 ##########

# mouse
m005_tpm_matrix<- read.table("M005_TPM.txt")
m005_annotation_df<- read.table("M005_METADATA.txt")

colors<-c("#ffa90c","#006400","red","grey")
names(colors)<- c("Human","Mouse","High_for_both","Low QC")

pdf("Human_mouse_counts.pdf",width = 8, height = 8,onefile=F)
ggscatter(m005_annotation_df, x="Mouse",y="Human",color="Classification", palette = colors)
dev.off()


m005_tpm_human_reads<- (m005_tpm_matrix[grep("GRCh38",rownames(m005_tpm_matrix)),]) ## human genes only
rownames(m005_tpm_human_reads)<- substr(rownames(m005_tpm_human_reads),8,30)

mouse_human_genes_ordered<- names(sort(rowMeans(m005_tpm_human_reads),decreasing=T))

mouse_cells_ordered_human<- rownames(m005_annotation_df)[m005_annotation_df$Classification=="Mouse"]
#mouse_cells_ordered_human<- names(sort(colMeans(m005_tpm_human_reads[mouse_human_genes_ordered[1:50],mouse_cells_ordered_human]),decreasing = T))

p1<-ggboxplot(m005_annotation_df[m005_annotation_df$Classification!="High_for_both" & !is.na(m005_annotation_df$Classification),], facet.by="Classification",x="Type", y="log_counts", fill="Type",palette = c("lightgrey","chartreuse"))+NoLegend()+theme(axis.text.x = element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank())+ stat_compare_means(method = "t.test",label = "p.signif")+ylab("log(1+ transcripts)")+ ggtitle("GFP transcripts")

p2<- ggboxplot(m005_annotation_df[m005_annotation_df$Classification!="High_for_both" & !is.na(m005_annotation_df$Classification),], facet.by="Classification",x="Type", y="log_mouse", fill="Type",palette = c("lightgrey","chartreuse"))+theme(axis.text.x = element_blank(),axis.title=element_blank(),axis.ticks.x=element_blank(),legend.position = "right")+ stat_compare_means(method = "t.test",label = "p.signif")+ylab("log(1+ transcripts)")+ ggtitle("mm10 detected genes")


pdf("gfp_mm10_transcript_counts.pdf",width = 8, height = 6,onefile=F)
egg::ggarrange(p1,p2,ncol=2)
dev.off()


# # compare to control
# 
# colnames(merged_tpm_matrix)<- str_replace_all(colnames(merged_tpm_matrix), "[.]","-")
# 
# mgg23_11A_tpm<- merged_tpm_matrix[,rownames(merged_annotation_df[merged_annotation_df$Sample=="MGG23"& merged_annotation_df$Organoid_line=="11A",])]
# 
# mgg23_df_pseudo_bulk<- data.frame(row.names = rownames(mgg23_11A_tpm))
# 
# mgg23_df_pseudo_bulk[,"GFP_pos"]<- log2(1+rowMeans(mgg23_11A_tpm[,colnames(mgg23_11A_tpm) %in% rownames(merged_annotation_df)[merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_POS"]]))
# 
# mgg23_df_pseudo_bulk[,"GFP_neg"]<- log2(1+rowMeans(mgg23_11A_tpm[,colnames(mgg23_11A_tpm) %in% rownames(merged_annotation_df)[merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_NEG"]]))
# 
# mgg23_df_pseudo_bulk[,"Avg_malignant"]<-log2(1+rowMeans(mgg23_11A_tpm[,colnames(mgg23_11A_tpm)%in% rownames(merged_annotation_df)[merged_annotation_df$CNA=="Malignant" ]]))


GFP_pos_only_df<- m005_annotation_df[m005_annotation_df$Type=="GFP_POS",]


GFP_pos_only_df$Cell_type<- factor(GFP_pos_only_df$Cell_type,levels=c("Malignant","IP","CPN","CFuPN" ,"Immature.interneurons", "oRG.Astroglia", "Cycling.oRG.Astroglia"))

p1<-ggboxplot(GFP_pos_only_df, x="Cell_type",y="log_mouse",add="jitter",add.params = list(size=0.2))+rotate_x_text(angle = 45, hjust = 1)+xlab("")+ylab("mouse genes detected (GFP+ - median GFP-)")  

p2<-ggboxplot(GFP_pos_only_df, x="Cell_type",y="log_counts",add="jitter",add.params = list(size=0.2))+rotate_x_text(angle = 45, hjust = 1)+xlab("")+ylab("GFP transcripts detected (GFP+ - median GFP-)")  


pdf("M005_gfp_mouse.pdf",width = 6, height = 10,onefile=F)
egg::ggarrange(p2+NoLegend()+theme(axis.text.x = element_blank()),p1,nrow=2)
dev.off()

## highest expressed malignant genes in GFP+ vs GFP- 

mouse_cols<-rownames(m005_annotation_df)[m005_annotation_df$Classification=="Mouse"]
human_cols<-rownames(m005_annotation_df)[m005_annotation_df$Classification=="Human"]

colnames(m005_tpm_matrix)<- str_replace_all(colnames(m005_tpm_matrix),"[.]","-")

m005_tpm_mouse_reads<- m005_tpm_matrix[,mouse_cols]
m005_tpm_mouse_reads<- count_to_tpm(m005_tpm_mouse_reads[grep("mm10",rownames(m005_tpm_mouse_reads)),]) ## mouse cells with mouse genes only
rownames(m005_tpm_mouse_reads)<- substr(rownames(m005_tpm_mouse_reads),8,30)

mouse_genes_ordered<- names(sort(rowMeans(m005_tpm_mouse_reads),decreasing=T))

GFP_POS_m005<- rownames(m005_annotation_df)[m005_annotation_df$Type=="GFP_POS"]
GFP_NEG_m005<- rownames(m005_annotation_df)[m005_annotation_df$Type=="GFP_NEG"]


m005_tpm_mouse_reads_human_cells<- m005_tpm_matrix[,human_cols] ## human cells with mouse genes
m005_tpm_mouse_reads_human_cells<- m005_tpm_mouse_reads_human_cells[grep("mm10",rownames(m005_tpm_mouse_reads_human_cells)),]
m005_tpm_mouse_reads_human_cells<- (m005_tpm_mouse_reads_human_cells)

m005_tpm_mouse_reads_human_cells<-m005_tpm_mouse_reads_human_cells[,!is.na(colMeans(m005_tpm_mouse_reads_human_cells))]

rownames(m005_tpm_mouse_reads_human_cells)<- substr(rownames(m005_tpm_mouse_reads_human_cells),8,30)
colnames(m005_tpm_mouse_reads_human_cells)<- str_replace(colnames(m005_tpm_mouse_reads_human_cells),"[.]","-")

human_df_pseudo_bulk<- data.frame(row.names = rownames(m005_tpm_mouse_reads_human_cells)) ##pseudo bulk of human cells with mouse reads

human_df_pseudo_bulk[,"GFP_pos"]<- log2(1+rowMeans(m005_tpm_mouse_reads_human_cells[,colnames(m005_tpm_mouse_reads_human_cells)[colnames(m005_tpm_mouse_reads_human_cells) %in% GFP_POS_m005]]))
human_df_pseudo_bulk[,"GFP_neg"]<- log2(1+rowMeans(m005_tpm_mouse_reads_human_cells[,colnames(m005_tpm_mouse_reads_human_cells)[colnames(m005_tpm_mouse_reads_human_cells) %in% GFP_NEG_m005]]))

human_df_pseudo_bulk[,"Avg_malignant"]<- log2(1+rowMeans(m005_tpm_mouse_reads))

human_df_pseudo_bulk_cell_types<- data.frame(row.names = rownames(m005_tpm_mouse_reads_human_cells)) ## pseudo bulk of human cells with mouse reads per cell type

for(cell_type in levels(m005_annotation_df$Cell_type)[levels(m005_annotation_df$Cell_type)!="Malignant"])
{
  human_df_pseudo_bulk_cell_types[,paste0(cell_type,"_GFP+")]<- log2(1+rowMeans(m005_tpm_mouse_reads_human_cells[,colnames(m005_tpm_mouse_reads_human_cells)[colnames(m005_tpm_mouse_reads_human_cells) %in% rownames(m005_annotation_df)[m005_annotation_df$Cell_type==cell_type & m005_annotation_df$Type=="GFP_POS"]]]))
  human_df_pseudo_bulk_cell_types[,paste0(cell_type,"_GFP-")]<- log2(1+rowMeans(m005_tpm_mouse_reads_human_cells[,colnames(m005_tpm_mouse_reads_human_cells)[colnames(m005_tpm_mouse_reads_human_cells) %in% rownames(m005_annotation_df)[m005_annotation_df$Cell_type==cell_type & m005_annotation_df$Type=="GFP_NEG"]]]))
}


p1<-basic_heatmap_wrapper(human_df_pseudo_bulk_cell_types[mouse_genes_ordered[1:100],]-rowMeans(human_df_pseudo_bulk_cell_types[mouse_genes_ordered[1:100],]),showAllGenes = F,showAllCells = T,limits = 1,ylab="rank 1-100")

p2<-basic_heatmap_wrapper(human_df_pseudo_bulk_cell_types[mouse_genes_ordered[101:1000],]-rowMeans(human_df_pseudo_bulk_cell_types[mouse_genes_ordered[101:1000],]),showAllGenes = F,showAllCells = T,limits = 1,ylab="rank 101-1000")


pdf("m005_heatmap_ranked_mouse_genes_bulk_celltypes.pdf",width = 8, height = 10,onefile=F)
egg::ggarrange(p2+theme(axis.text.x = element_blank(),axis.title.x=element_blank())+NoLegend(),p1+xlab(""),left = "           ranked malignant genes (mm10)")
dev.off()


p1<-basic_heatmap_wrapper(human_df_pseudo_bulk[mouse_genes_ordered[1:100],1:2]-rowMeans(human_df_pseudo_bulk[mouse_genes_ordered[1:100],1:2]),showAllGenes = F,showAllCells = T,limits = 1,ylab="1-100")

p2<-basic_heatmap_wrapper(human_df_pseudo_bulk[mouse_genes_ordered[101:1000],1:2]-rowMeans(human_df_pseudo_bulk[mouse_genes_ordered[101:1000],1:2]),showAllGenes = F,showAllCells = T,limits = 1,ylab="101-1000")


pdf("m005_heatmap_ranked_mouse_genes_bulk.pdf",width = 4, height = 10,onefile=F)
egg::ggarrange(p2+theme(axis.text.x = element_blank(),axis.title.x=element_blank())+NoLegend(),p1+xlab(""),left="           ranked malignant genes (mm10)")
dev.off()

# sc
centered_m005_tpm_mouse_reads_human_cells<- m005_tpm_mouse_reads_human_cells-rowMeans(m005_tpm_mouse_reads_human_cells)

GFP_POS_m005_o<- GFP_POS_m005[GFP_POS_m005 %in% colnames(centered_m005_tpm_mouse_reads_human_cells)]
GFP_POS_m005_o<- GFP_POS_m005_o[order(colMeans(centered_m005_tpm_mouse_reads_human_cells[mouse_genes_ordered[1:100],GFP_POS_m005_o]),decreasing=T)]


GFP_NEG_m005_o<- GFP_NEG_m005[GFP_NEG_m005 %in% colnames(centered_m005_tpm_mouse_reads_human_cells)]
GFP_NEG_m005_o<- GFP_NEG_m005_o[order(colMeans(centered_m005_tpm_mouse_reads_human_cells[mouse_genes_ordered[1:100],GFP_NEG_m005_o]),decreasing=T)]

mouse_genes_ordered<- mouse_genes_ordered[!mouse_genes_ordered %in%  c("Gm42418","Hsp90ab1")]

p1<-basic_heatmap_wrapper(centered_m005_tpm_mouse_reads_human_cells[mouse_genes_ordered[1:100],c(GFP_POS_m005_o,GFP_NEG_m005_o)],limits = 1,ylab="1-500")+geom_vline(xintercept = length(GFP_POS_m005_o)+0.5)

p2<- basic_heatmap_wrapper(centered_m005_tpm_mouse_reads_human_cells[mouse_genes_ordered[101:1000],c(GFP_POS_m005_o,GFP_NEG_m005_o)],limits = 1,ylab="501-2000")+geom_vline(xintercept = length(GFP_POS_m005_o)+0.5)

pdf("m005_heatmap_ranked_mouse_genes_sc.pdf",width = 4, height = 10,onefile=F)
egg::ggarrange(p2+theme(axis.text.x = element_blank(),axis.title.x=element_blank())+NoLegend(),p1+xlab(""),left="           ranked malignant genes (mm10)")
dev.off()


plot_scatter<-melt(as.matrix(human_df_pseudo_bulk[,1:2]))

plot_scatter$Avg_malignant<- human_df_pseudo_bulk[as.character(plot_scatter$Var1),]$Avg_malignant
plot_scatter$genes_expression<-""
plot_scatter[(plot_scatter$Var1) %in% mouse_highest_genes,]$genes_expression<-"high"


genes_to_label<- unique(as.character(plot_scatter[plot_scatter$Avg_malignant>4 & plot_scatter$value>0.3,]$Var1))

genes_to_label<- genes_to_label[-grep("Rp",genes_to_label)]
genes_to_label<- genes_to_label[-grep("mt",genes_to_label)]



plot_scatter<- plot_scatter[!(plot_scatter$Var1) %in%  c("Gm42418","Hsp90ab1"),]

reg1<-lm(value ~ Avg_malignant, data = plot_scatter[plot_scatter$Var2=="GFP_pos" & plot_scatter$genes_expression=="high",])
coeff1=coefficients(reg1)


reg2<-lm(value ~ Avg_malignant, data = plot_scatter[plot_scatter$Var2=="GFP_neg" & plot_scatter$genes_expression=="high",])
coeff2=coefficients(reg2)

sp<- ggplot(data=plot_scatter, aes(x=Avg_malignant, y=value,color=Var2)) + geom_point()+ scale_color_manual(values=c("chartreuse","darkgrey"))

pdf("m005_dotplot_mouse_genes.pdf",width = 8, height = 7,onefile=F)
sp+ geom_abline(intercept = coeff1[1], slope = coeff1[2],color="#6BC71D", linetype="dashed", size=1)+ geom_abline(intercept = coeff2[1], slope = coeff2[2],color="#777777", linetype="dashed", size=1)+ theme_classic() + ylab("Avg Exp. nonmalignant")+ xlab("Avg Exp. malignant")+
  stat_cor(data= plot_scatter[plot_scatter$genes_expression=="high",], aes(color = Var2),method = "pearson")
dev.off()


top.genes<- rownames(m005_tpm_mouse_reads)[log2(1+rowMeans(m005_tpm_mouse_reads))>=4]
log_m005_tpm_mouse_reads<- log2(1+m005_tpm_mouse_reads[top.genes,]/10)
colnames(log_m005_tpm_mouse_reads)<- str_replace_all(colnames(log_m005_tpm_mouse_reads),"[.]","-")


#mouse_malig<- readRDS( "mouse_malig.rds")
mouse_malig<- read.table("TableS4_DeNovoProgramsMouse.txt")



NMF_heatmap1<-NMF_heatmap(log_expression_matrix = log_m005_tpm_mouse_reads, NMF_programs = mouse_malig,showAllGenes1 = T)

mouse_GBM<-lapply(Signatures_GBM, function(x){human_to_mouse(x)$mouse}) #converting neftel states to mouse
mouse_GBM<- lapply(mouse_GBM, function(x)x[x %in% rownames(log_m005_tpm_mouse_reads)]) 

GBM_scores<-sigScores(log_m005_tpm_mouse_reads, mouse_GBM) # scoring m005 mouse cells for mouse neftel states

scores_malig_p<-basic_heatmap_wrapper(t(GBM_scores[NMF_heatmap1$cells_order,c("OPC","NPC1","NPC2",'AC',"MES1","MES2","G1S","G2M")]),showAllGenes = T,showAllCells = F,limits = 2)

png("Mouse_malignant_denovostates.png",width = 12, height = 10,units = 'in', res = 300,type="cairo")
egg::ggarrange(NMF_heatmap1$heatmap+xlab(""),scores_malig_p,nrow=2, heights = c(50,10) )
dev.off()


### scoring for canonical programs human non malignant cells from m005
temp_seu<- seu_mouse_model[,seu_mouse_model$Classification=="Human"& seu_mouse_model$Cell_type!="Malignant"]

table(temp_seu@meta.data[,c("Cell_type","Type")])

m005_nonmalignant_tpm<- data.frame(m005_tpm_matrix[,rownames(m005_annotation_df)[m005_annotation_df$Classification=="Human"&m005_annotation_df$Cell_type!="Malignant"]])
colnames(m005_nonmalignant_tpm)<- str_replace_all(colnames(m005_nonmalignant_tpm),"[.]","-")

GFP_POS_nonmal<-colnames(m005_nonmalignant_tpm)[colnames(m005_nonmalignant_tpm) %in% GFP_POS_m005]
GFP_NEG_nonmal<-colnames(m005_nonmalignant_tpm)[colnames(m005_nonmalignant_tpm) %in% GFP_NEG_m005]

#human genes
m005_nonmalignant_tpm<- m005_nonmalignant_tpm[grep("GR",rownames(m005_nonmalignant_tpm)),]
rownames(m005_nonmalignant_tpm)<- substr(rownames(m005_nonmalignant_tpm),8,30)

log_m005_nonmalignant_tpm<- log2(1+m005_nonmalignant_tpm/10)

top.genes<- rownames(m005_nonmalignant_tpm)[log2(1+rowMeans(m005_nonmalignant_tpm))>=4]

Signatures_GBM<- lapply(scalop::Signatures_GBM, function(x){x[x %in% top.genes]})

neftel_states_nonmalig<-sigScores(m = log_m005_nonmalignant_tpm,sigs = Signatures_GBM[c(1,4:8)])

melt_test_df<-melt(as.matrix(t(neftel_states_nonmalig)))

melt_test_df$Var2<- as.character(melt_test_df$Var2)
melt_test_df$Var2<- str_replace_all(melt_test_df$Var2,"[.]",'-')
melt_test_df[melt_test_df$Var2 %in% GFP_NEG_nonmal,"GFP"]<-"NEG"
melt_test_df[melt_test_df$Var2 %in% GFP_POS_nonmal,"GFP"]<-"POS"

melt_test_df$GFP<- factor(melt_test_df$GFP,levels=c("POS","NEG"))
melt_test_df[,"Cell_type"]<- m005_annotation_df[as.character(melt_test_df$Var2),]$Cell_type

pdf("m005_canonical_human_nonmalignant.pdf",width = 8, height = 6,onefile=F)
(ggboxplot(melt_test_df,x="Var1",y="value",fill = "GFP",add.params = list(size=0.5),palette = c("chartreuse","darkgrey"))+ stat_compare_means(aes(group = GFP),method = "t.test",label = "p.signif", method.args = list(adjust.method = "bonferroni"),hide.ns=T)+xlab("GBM state")+ylab("score"))
dev.off()


melt_test_df$Cell_type<- factor(melt_test_df$Cell_type, levels=c ("oRG.Astroglia" ,"Cycling.oRG.Astroglia", "IP","Immature.interneurons","CPN","CFuPN"))

pdf("m005_canonical_human_nonmalignant_celltype.pdf",width = 10, height = 6,onefile=F)
(ggboxplot(melt_test_df,x="Var1",y="value",fill = "GFP",add.params = list(size=0.5),palette = c("chartreuse","darkgrey"),facet.by = "Cell_type")+ stat_compare_means(aes(group = GFP),method = "t.test",label = "p.signif", method.args = list(alternative = "greater",adjust.method = "bonferroni"))+xlab("GBM state")+ylab("score"))
dev.off()

## genes detected from each organism
m005_annotation_df$Type_2<-paste(m005_annotation_df$Classification,m005_annotation_df$Type)

m005_annotation_df2<- m005_annotation_df[m005_annotation_df$Classification!="High_for_both",]

m005_annotation_df2$Type_2<- factor(m005_annotation_df2$Type_2, levels=c("Human GFP_NEG","Human GFP_POS","Mouse GFP_POS"))

p1<-ggboxplot(m005_annotation_df2,y="Mouse",x="Type_2",add="jitter")+ylab("Mouse genes detected")+xlab("")

p2<-ggboxplot(m005_annotation_df2,y="Human",x="Type_2",add="jitter")+ylab("Human genes detected")+xlab("")

pdf("Final_paper_plots/S4_genesdetected.pdf",width = 8, height = 6,onefile=F)
egg::ggarrange(p1,p2,ncol=2)
dev.off()



# human samples

for(model in unique(merged_annotation_df$Sample))
{
  temp_seu<- merged_tpm_matrix[,rownames(merged_annotation_df)[merged_annotation_df$Sample==model]]
  
  mal.cells<-rownames(merged_annotation_df)[merged_annotation_df$Sample==model & merged_annotation_df$CNA=="Malignant"]
  nonmal.cells<-rownames(merged_annotation_df)[merged_annotation_df$Sample==model & merged_annotation_df$CNA=="Non-malignant"]
  
  GFP_POS_nonmal<-rownames(merged_annotation_df)[merged_annotation_df$Sample==model & merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_POS"]
  GFP_NEG_nonmal<-rownames(merged_annotation_df)[merged_annotation_df$Sample==model & merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_NEG"]
  
  temp_pseudo_bulk_df<- data.frame(row.names = rownames(temp_seu))
  temp_pseudo_bulk_df[,"Avg_GFP_POS"]<- log2(1+rowMeans(temp_seu[,GFP_POS_nonmal]))
  temp_pseudo_bulk_df[,"Avg_GFP_NEG"]<- log2(1+rowMeans(temp_seu[,GFP_NEG_nonmal]))
  temp_pseudo_bulk_df[,"Avg_malignant"]<- log2(1+rowMeans(temp_seu[,mal.cells]))
  temp_pseudo_bulk_df[,"Avg_nonmalignant"]<- log2(1+rowMeans(temp_seu[,nonmal.cells]))
  
  top_malignant_genes<- get_top_genes(temp_seu[,mal.cells])
  
  M2 <- lm(data = temp_pseudo_bulk_df[,], formula = Avg_malignant~Avg_nonmalignant)
  
  
  expressed_in_GFP_pos  <- rownames(temp_pseudo_bulk_df)[temp_pseudo_bulk_df$Avg_malignant>=4][(M2$residuals[  rownames(temp_pseudo_bulk_df)[temp_pseudo_bulk_df$Avg_malignant>=4]]>=quantile(M2$residuals[  rownames(temp_pseudo_bulk_df)[temp_pseudo_bulk_df$Avg_malignant>=4]],0.9))]
  
  
  
  temp_pseudo_bulk_df$gene<-rownames(temp_pseudo_bulk_df)
  temp_pseudo_bulk_df[,"color"]<-"none"
  temp_pseudo_bulk_df[expressed_in_GFP_pos,"color"]<-"residual"
  
  p1<-ggscatter(temp_pseudo_bulk_df[,], y="Avg_malignant",x="Avg_nonmalignant",color="color",add = "reg.line",conf.int = TRUE, palette = c("grey","red"), add.params = list(color = "black", fill = "lightgray"))+ xlab("Avg expression organoid cells")+ylab("Avg expresion in malignant cells")+ NoLegend() + geom_point(data=temp_pseudo_bulk_df[temp_pseudo_bulk_df$color=="residual",], x=temp_pseudo_bulk_df[temp_pseudo_bulk_df$color=="residual",]$Avg_nonmalignant, y=temp_pseudo_bulk_df[temp_pseudo_bulk_df$color=="residual",]$Avg_malignant,color="red")+ggtitle("")
  
  melted_matrix<- melt(as.matrix(temp_pseudo_bulk_df[,1:2]))
  melted_matrix$color<-"All genes"
  melted_matrix[melted_matrix$Var1 %in% expressed_in_GFP_pos,"color"]="Malignant-enriched genes"
  
  melted_matrix$color<- factor(melted_matrix$color,levels=c("Malignant-enriched genes","All genes"))
  
  melted_matrix$Var2<- str_remove_all(melted_matrix$Var2,"Avg_")
  
  p2<- ggboxplot(melted_matrix, x="Var2",y="value",color="Var2",facet.by = "color",add="jitter",palette = c("chartreuse","darkgrey"))+xlab("")+ylab("Avg expression organoid cells")+ stat_compare_means(method = "t.test",label = "p.signif",paired=T,method.args = list(alternative="less"))+NoLegend()
  
  png(paste0(model,"_malig_enrich_genes.png"),width = 12, height = 6,units = 'in', res = 300,type="cairo")
  print(egg::ggarrange(p1,p2,ncol=2,top = model))
  dev.off()
  
  
} 


##combine all models




mal.cells<-rownames(merged_annotation_df)[ merged_annotation_df$CNA=="Malignant"]
nonmal.cells<-rownames(merged_annotation_df)[ merged_annotation_df$CNA=="Non-malignant"]

GFP_POS_nonmal<-rownames(merged_annotation_df)[merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_POS"]
GFP_NEG_nonmal<-rownames(merged_annotation_df)[merged_annotation_df$CNA=="Non-malignant" & merged_annotation_df$Type=="GFP_NEG"]


temp_pseudo_bulk_df<- data.frame(row.names = rownames(merged_tpm_matrix))
temp_pseudo_bulk_df[,"Avg_GFP_POS"]<- log2(1+rowMeans(merged_tpm_matrix[,GFP_POS_nonmal]))
temp_pseudo_bulk_df[,"Avg_GFP_NEG"]<- log2(1+rowMeans(merged_tpm_matrix[,GFP_NEG_nonmal]))
temp_pseudo_bulk_df[,"Avg_malignant"]<- log2(1+rowMeans(merged_tpm_matrix[,mal.cells]))
temp_pseudo_bulk_df[,"Avg_nonmalignant"]<- log2(1+rowMeans(merged_tpm_matrix[,nonmal.cells]))

top_malignant_genes<- get_top_genes(merged_tpm_matrix[,mal.cells])

M2 <- lm(data = temp_pseudo_bulk_df[,], formula = Avg_malignant~Avg_nonmalignant)

expressed_in_GFP_pos  <- rownames(temp_pseudo_bulk_df)[temp_pseudo_bulk_df$Avg_malignant>=4][(M2$residuals[  rownames(temp_pseudo_bulk_df)[temp_pseudo_bulk_df$Avg_malignant>=4]]>=quantile(M2$residuals[  rownames(temp_pseudo_bulk_df)[temp_pseudo_bulk_df$Avg_malignant>=4]],0.9))]



temp_pseudo_bulk_df$gene<-rownames(temp_pseudo_bulk_df)
temp_pseudo_bulk_df[,"color"]<-"none"
temp_pseudo_bulk_df[expressed_in_GFP_pos,"color"]<-"residual"

p1<-ggscatter(temp_pseudo_bulk_df[,], y="Avg_malignant",x="Avg_nonmalignant",color="color",add = "reg.line",conf.int = TRUE, palette = c("grey","red"), add.params = list(color = "black", fill = "lightgray"))+ xlab("Avg expression organoid cells")+ylab("Avg expresion malignant cells")+ NoLegend() + geom_point(data=temp_pseudo_bulk_df[temp_pseudo_bulk_df$color=="residual",], x=temp_pseudo_bulk_df[temp_pseudo_bulk_df$color=="residual",]$Avg_nonmalignant, y=temp_pseudo_bulk_df[temp_pseudo_bulk_df$color=="residual",]$Avg_malignant,color="red")+ggtitle(paste(""))

melted_matrix<- melt(as.matrix(temp_pseudo_bulk_df[,1:2]))
melted_matrix$color<-"All genes"
melted_matrix[melted_matrix$Var1 %in% expressed_in_GFP_pos,"color"]="Malignant-enriched genes"

melted_matrix$color<- factor(melted_matrix$color,levels=c("Malignant-enriched genes","All genes"))

p2<- ggboxplot(melted_matrix, x="color",y="value",fill="Var2",add.params = list(size=0.5),palette = c("chartreuse","darkgrey"))+xlab("")+ylab("Avg expression organoid cells")+ stat_compare_means(method = "t.test",label = "p.signif")

pdf("malig_enrich_genes_GFPpos.pdf",width = 10, height = 5,onefile=F)
print(egg::ggarrange(p1,p2,ncol=2))
dev.off()

celltype_pseudo_bulk_df<- data.frame(row.names=rownames(temp_pseudo_bulk_df))



for(cell_type in  unique(merged_annotation_df$Cell_type)[1:6] )
{
  
 
  celltype_pseudo_bulk_df[,paste0(cell_type,"_GFP_POS")]<- log2(1+rowMeans(merged_tpm_matrix[,rownames(merged_annotation_df)[merged_annotation_df$Type=="GFP_POS" & merged_annotation_df$Cell_type==cell_type]]))
  
  celltype_pseudo_bulk_df[,paste0(cell_type,"_GFP_NEG")]<- log2(1+rowMeans(merged_tpm_matrix[,rownames(merged_annotation_df)[merged_annotation_df$Type=="GFP_NEG" & merged_annotation_df$Cell_type==cell_type]]))
}


melted_matrix<- melt(as.matrix(celltype_pseudo_bulk_df))
melted_matrix$color<-"All genes"
melted_matrix[melted_matrix$Var1 %in% expressed_in_GFP_pos,"color"]="Malignant-enriched genes"

melted_matrix[grep("POS",melted_matrix$Var2),"GFP"]<-"POS"
melted_matrix[grep("NEG",melted_matrix$Var2),"GFP"]<-"NEG"
melted_matrix$GFP<- factor(melted_matrix$GFP,levels=c("POS","NEG"))

melted_matrix$cell_type<- split_names(string = as.character(melted_matrix$Var2),char = "_",position = 1)


melted_matrix$color<- factor(melted_matrix$color,levels=c("Malignant-enriched genes","All genes"))


pdf("malig_enrich_genes_celltype.pdf",width = 10, height = 8,onefile=F)
print(ggboxplot(melted_matrix, x="color",y="value",fill="GFP",add.params=list(size=0.5),palette = c("chartreuse","darkgrey"),facet.by = "cell_type")+xlab("")+ylab("Avg expression in organoid cells")+ stat_compare_means(method = "t.test",label = "p.signif"))
dev.off()



pseudo_bulk_df_human<- data.frame(row.names = rownames(merged_tpm_matrix))



for(cell_type in unique(merged_annotation_df$Cell_type)[1:7])
{
  temp_GFP_POS_nonmal<- GFP_POS_nonmal[merged_annotation_df[GFP_POS_nonmal,]$Cell_type == cell_type]
  temp_GFP_NEG_nonmal<- GFP_NEG_nonmal[merged_annotation_df[GFP_NEG_nonmal,]$Cell_type == cell_type]
  
  pseudo_bulk_df_human[,paste0(cell_type,"_Avg_GFP+")]<- log2(1+rowMeans(merged_tpm_matrix[,temp_GFP_POS_nonmal])) -  rowMeans(cbind(log2(1+rowMeans(merged_tpm_matrix[,temp_GFP_POS_nonmal])),log2(1+rowMeans(merged_tpm_matrix[,temp_GFP_NEG_nonmal]))))
  
  pseudo_bulk_df_human[,paste0(cell_type,"Avg_GFP-")]<- log2(1+rowMeans(merged_tpm_matrix[,temp_GFP_NEG_nonmal]))-  rowMeans(cbind(log2(1+rowMeans(merged_tpm_matrix[,temp_GFP_POS_nonmal])),log2(1+rowMeans(merged_tpm_matrix[,temp_GFP_NEG_nonmal]))))
}


program_plots_list<- list()

for(program in names(Signatures_GBM))
{
  genes<-scalop::Signatures_GBM[[program]]
  genes<- genes[order(rowMeans(pseudo_bulk_df_human[genes,grep("GFP+",colnames(pseudo_bulk_df_human))])-rowMeans(pseudo_bulk_df_human[genes,grep("GFP-",colnames(pseudo_bulk_df_human))]),decreasing=T)]
  program_plots_list[[program]]<-basic_heatmap_wrapper(pseudo_bulk_df_human[genes,],limits = 0.5, showAllGenes = T,showAllCells = T)+xlab("")+ggtitle(program)+geom_vline(xintercept = 2.5)+geom_vline(xintercept = 4.5)+geom_vline(xintercept = 6.5)+geom_vline(xintercept = 8.5)+geom_vline(xintercept = 10.5)+geom_vline(xintercept = 12.5)+geom_vline(xintercept = 14.5)
}

pdf("canonicalGBM_full_in_human_samples.pdf",width = 15, height = 10,onefile=F)
egg::ggarrange(program_plots_list$AC+NoLegend()+theme(axis.text.x=element_blank()),program_plots_list$MES1+theme(axis.text.x=element_blank())+NoLegend(),program_plots_list$MES2+theme(axis.text.x=element_blank())+NoLegend(),program_plots_list$NPC1+NoLegend(),program_plots_list$NPC2+NoLegend(),program_plots_list$OPC,ncol=3,top="Human cells from cell line experiments + Human genes")
dev.off()




test_df<-sigScores(m = merged_tpm_matrix[,],sigs = Signatures_GBM[c("AC","MES1","MES2","OPC","NPC1","NPC2")])

melt_test_df<-melt(as.matrix(t(test_df)))

melt_test_df$Var2<- as.character(melt_test_df$Var2)
melt_test_df$Var2<- str_replace_all(melt_test_df$Var2,"[.]",'-')
melt_test_df[melt_test_df$Var2 %in% GFP_NEG_nonmal,"GFP"]<-"NEG"
melt_test_df[melt_test_df$Var2 %in% GFP_POS_nonmal,"GFP"]<-"POS"

melt_test_df$GFP<- factor(melt_test_df$GFP,levels=c("POS","NEG"))

melt_test_df<- melt_test_df[!is.na(melt_test_df$GFP),]

pdf("canonicalGBM_in_human_samples.pdf",width = 6.5, height = 5,onefile=F)
ggboxplot(melt_test_df,x="Var1",y="value",fill = "GFP",add.params = list(size=0.5),palette = c("chartreuse","darkgrey"))+ stat_compare_means(aes(group = GFP),method = "t.test",label = "p.signif")+ylab("Score")+xlab("")
dev.off()


melt_test_df$Sample<- split_names(melt_test_df$Var2,position = 1,char = "_")


pdf("canonicalGBM_in_human_samples_split.pdf",width = 10, height = 10,onefile=F)
(ggboxplot(melt_test_df,x="Var1",y="value",fill = "GFP",facet.by = "Sample",add.params = list(size=0.5),palette = c("chartreuse","darkgrey"))+ stat_compare_means(aes(group = GFP),method = "t.test",label = "p.signif")+xlab("GBM state scores")+ylab("log expression"))
dev.off()

melt_test_df[,"Cell_type"]<- merged_annotation_df[(melt_test_df$Var2),]$Cell_type


pdf("canonicalGBM_in_human_samples_celltype.pdf",width = 10, height = 10,onefile=F)
(ggboxplot(melt_test_df,x="Var1",y="value",fill = "GFP",facet.by="Cell_type",palette = c("chartreuse","darkgrey"))+ stat_compare_means(aes(group = GFP),method = "t.test",label = "p.signif",hide.ns=T)+xlab("GBM state scores")+ylab("log expression"))
dev.off()

## EV markers
test_df<-merged_tpm_matrix[c("CD63", "CD81","CD9", "CD82",  "CD151", "SDCBP" , "PDCD6IP"),]

melt_test_df<-melt(as.matrix((test_df)))

melt_test_df$Var2<- as.character(melt_test_df$Var2)
melt_test_df$Var2<- str_replace_all(melt_test_df$Var2,"[.]",'-')
melt_test_df[melt_test_df$Var2 %in% GFP_NEG_nonmal,"GFP"]<-"NEG"
melt_test_df[melt_test_df$Var2 %in% GFP_POS_nonmal,"GFP"]<-"POS"

melt_test_df$GFP<- factor(melt_test_df$GFP,levels=c("POS","NEG"))

pdf("Human_samples_EV.pdf",width = 6.5, height = 5,onefile=F)
(ggboxplot(melt_test_df,x="Var1",y="value",fill = "GFP",add.params = list(size=0.5),palette = c("chartreuse","darkgrey"))+ stat_compare_means(aes(group = GFP),method = "t.test",label = "p.signif")+xlab("Exosome markers")+ylab("log expression"))
dev.off()


melt_test_df[,"Cell_type"]<- merged_annotation_df[(melt_test_df$Var2),]$Cell_type


pdf("Human_samples_EV_celltype.pdf",width = 10, height = 10,onefile=F)
(ggboxplot(melt_test_df,x="Var1",y="value",fill = "GFP",facet.by="Cell_type",palette = c("chartreuse","darkgrey"))+ stat_compare_means(aes(group = GFP),method = "t.test",label = "p.signif")+xlab("Exosome markers")+ylab("log expression"))
dev.off()


### DE gfp + vs -

pseudo_bulk_df<- data.frame(row.names = rownames(merged_tpm_matrix))
GFP_POS_CELLS<- rownames(merged_annotation_df)[merged_annotation_df$Type=="GFP_POS"& merged_annotation_df$Cell_type!="Malignant"]
GFP_NEG_CELLS<- rownames(merged_annotation_df)[merged_annotation_df$Type=="GFP_NEG"& merged_annotation_df$Cell_type!="Malignant"]


pseudo_bulk_df_cell_type<-  data.frame(row.names = rownames(merged_tpm_matrix))

cell_types<- c(unique(merged_annotation_df$Cell_type), "all")
cell_types<- cell_types[cell_types!="Malignant"]

for(sample in unique(merged_annotation_df$Sample))
{
  for(organoid_line in unique(merged_annotation_df$Organoid_line))
  {
    if(sum(merged_annotation_df$Sample==sample & merged_annotation_df$Organoid_line==organoid_line)<100)
      next
    
    temp_expression_df<- data.frame(merged_tpm_matrix[,rownames(merged_annotation_df)[merged_annotation_df$Sample==sample  & merged_annotation_df$Organoid_line==organoid_line]])
    
    colnames(temp_expression_df)<- str_replace(colnames(temp_expression_df),"[.]","-")
    
    
    pseudo_bulk_df[,paste(sample,"all",organoid_line,"GFP_POS",sep="_")]<- log2(1+rowMeans(temp_expression_df[rownames(pseudo_bulk_df),colnames(temp_expression_df)[colnames(temp_expression_df) %in% GFP_POS_CELLS]])) 
    
    pseudo_bulk_df[,paste(sample,"all",organoid_line,"GFP_NEG",sep="_")]<- log2(1+rowMeans(temp_expression_df[rownames(pseudo_bulk_df),colnames(temp_expression_df)[colnames(temp_expression_df) %in% GFP_NEG_CELLS]])) 
    
    temp_mean<- rowMeans(pseudo_bulk_df[,grep(paste(sample,"all",organoid_line,sep="_"),colnames(pseudo_bulk_df))])
    
    pseudo_bulk_df[,paste(sample,"all",organoid_line,"GFP_POS",sep="_")]<- pseudo_bulk_df[,paste(sample,"all",organoid_line,"GFP_POS",sep="_")]- temp_mean
    pseudo_bulk_df[,paste(sample,"all",organoid_line,"GFP_NEG",sep="_")]<- pseudo_bulk_df[,paste(sample,"all",organoid_line,"GFP_NEG",sep="_")]- temp_mean
    
    
    for(cell_type in cell_types)
    {
      
      
      if(sum(merged_annotation_df$Sample==sample &merged_annotation_df$Cell_type==cell_type & merged_annotation_df$Organoid_line==organoid_line)<100)
        next
      
      temp_expression_df<- data.frame(merged_tpm_matrix[,rownames(merged_annotation_df)[merged_annotation_df$Sample==sample  & merged_annotation_df$Organoid_line==organoid_line & merged_annotation_df$Cell_type==cell_type]])
      
      colnames(temp_expression_df)<- str_replace(colnames(temp_expression_df),"[.]","-")
      pseudo_bulk_df[,paste(sample,cell_type,organoid_line,"GFP_POS",sep="_")]<- log2(1+rowMeans(temp_expression_df[rownames(pseudo_bulk_df),colnames(temp_expression_df)[colnames(temp_expression_df) %in% GFP_POS_CELLS]])) 
      
      pseudo_bulk_df[,paste(sample,cell_type,organoid_line,"GFP_NEG",sep="_")]<- log2(1+rowMeans(temp_expression_df[rownames(pseudo_bulk_df),colnames(temp_expression_df)[colnames(temp_expression_df) %in% GFP_NEG_CELLS]])) 
      
      temp_mean<- rowMeans(pseudo_bulk_df[,grep(paste(sample,cell_type,organoid_line,sep="_"),colnames(pseudo_bulk_df))])
      
      pseudo_bulk_df[,paste(sample,cell_type,organoid_line,"GFP_POS",sep="_")]<- pseudo_bulk_df[,paste(sample,cell_type,organoid_line,"GFP_POS",sep="_")]- temp_mean
      pseudo_bulk_df[,paste(sample,cell_type,organoid_line,"GFP_NEG",sep="_")]<- pseudo_bulk_df[,paste(sample,cell_type,organoid_line,"GFP_NEG",sep="_")]- temp_mean
      
      pseudo_bulk_df_cell_type[,paste(sample,cell_type,organoid_line,sep="_")]<- log2(1+rowMeans(temp_expression_df[rownames(pseudo_bulk_df),])) 
      
    }
    
    
    
  }
  
}


top.genes<- get_top_genes(merged_tpm_matrix)
sig_genes_list<- list()

volc_list<- list()

for( cell_type in  cell_types)
{
  temp_pseudo_bulk_df<- pseudo_bulk_df[,grep(paste0("_",cell_type),colnames(pseudo_bulk_df))]
  if(ncol(temp_pseudo_bulk_df)<3)
    next
  
  t_test <- apply(temp_pseudo_bulk_df, 1, function(x) tryCatch({
    t.test(x[grep("GFP_POS",colnames(temp_pseudo_bulk_df))], x[grep("GFP_NEG",colnames(temp_pseudo_bulk_df))],paired=TRUE)$p.value
  }, error = function(e) {NA}))
  
  fc <- rowMeans(temp_pseudo_bulk_df[,grep("GFP_POS",colnames(temp_pseudo_bulk_df))]) - rowMeans(temp_pseudo_bulk_df[,grep("GFP_NEG",colnames(temp_pseudo_bulk_df))])
  
  
  volc <- data.frame(cbind(t_test, fc))
  volc$Gene<- rownames(volc)
  colnames(volc)<- c("pvalue","log2FoldChange")
  
  volc_top<- volc[top.genes,]
  volc_list[[cell_type]]<-volc_top
  
  volc_top<- volc_top[order(volc_top$log2FoldChange, decreasing=T),]
  upgenes<- rownames(volc_top)[volc_top$log2FoldChange>.2 & volc_top$pvalue<0.05][1:50]
  
  volc_top<- volc_top[order(volc_top$log2FoldChange, decreasing=F),]
  downgenes<- rownames(volc_top)[volc_top$log2FoldChange<(-.2) & volc_top$pvalue<0.05][1:50]
  
  sig_genes_list[[paste(cell_type,"UP")]]<- upgenes[!is.na(upgenes)]
  sig_genes_list[[paste(cell_type,"DOWN")]]<- downgenes[!is.na(downgenes)]

  
  
}

volc_plot_list<- list()

for(name in names(volc_list))
{
  volc_all<- volc_list[[name]]
  volc_all<- volc_all[order(volc_all$log2FoldChange,decreasing = T),]
  volc_all$gene<-rownames(volc_all)
  options(ggrepel.max.overlaps = Inf)
  
  volc_all[unlist(scalop::Signatures_GBM[c("MES1","MES2")]),"color"]<-"MES/AC"
  volc_all[unlist(scalop::Signatures_GBM[c("NPC1","NPC2")]),"color"]<-"NPC/OPC"
  volc_all[unlist(scalop::Signatures_GBM[c("OPC")]),"color"]<-"NPC/OPC"
  volc_all[unlist(scalop::Signatures_GBM[c("AC")]),"color"]<-"MES/AC"
  
  
  volc_all$color<- factor(volc_all$color,levels=c("MES/AC","NPC/OPC"))
  

  
  sig_genes<- rownames(volc_all[volc_all$log2FoldChange>0.5 & -log10(volc_all$pvalue)>2,])
  sig_genes<- sig_genes[!is.na(sig_genes)]
  
  
  sig_genes<- sig_genes[1:10]
  
  sig_genes<- unique(c("NES", "S100A6","S100A10", "GAP43", sig_genes))
  
  sig_genes<- sig_genes[volc_all[sig_genes,]$pvalue<0.05 & volc_all[sig_genes,]$log2FoldChange>0.5 ]
  

  options(ggrepel.max.overlaps = Inf)
  volc_plot_list[[name]]<- ggplot(volc_all, aes(x=log2FoldChange, y=-log10(pvalue)))+ labs(x="log2(FoldChange)", y="-log10(p-value)") + geom_point(alpha = 0.6, size=.8,color="lightgrey") + ggtitle(paste(name,"GFP+ vs GFP-"))+ theme(axis.text=element_text(size=10),axis.title=element_text(size=11), plot.title = element_text(size=13), legend.title = element_text(size = 10),  legend.text = element_text(size = 10),legend.key=element_blank(), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), aspect.ratio = 1)+ geom_point(data=volc_all[c(scalop::Signatures_GBM$NPC1,scalop::Signatures_GBM$NPC2,scalop::Signatures_GBM$OPC),], colour="blue")+ geom_point(data=volc_all[c(scalop::Signatures_GBM$MES1,scalop::Signatures_GBM$MES2,scalop::Signatures_GBM$AC),], colour="red")#+ geom_text_repel(data=volc_all[sig_genes,],label=volc_all[sig_genes,]$gene)+geom_point(data=volc_all[sig_genes,],color="black")
  
  
  
}

pdf("GFPpos_neg_volc_plot_celltype.pdf",width = 10, height = 15,onefile=F)
volc_plot_list[[1]]+NoLegend()+volc_plot_list[[2]]+NoLegend()+volc_plot_list[[3]]+volc_plot_list[[4]]+NoLegend()+volc_plot_list[[5]]+NoLegend()+volc_plot_list[[6]]+NoLegend()+volc_plot_list[[7]]+NoLegend()
dev.off()

##ALL:

labels<- volc_all[volc_all$log2FoldChange>1 & volc_all$pvalue<0.05,]$gene

labels<- labels[labels %in% unlist(Signatures_GBM)]


pdf("GFPpos_neg_volc_plot.pdf",width = 8, height = 6,onefile=F)
ggplot(volc_all, aes(x=log2FoldChange, y=-log10(pvalue)))+ labs(x="log2(FoldChange)", y="-log10(p-value)") + geom_point(alpha = 0.6, size=1,color="grey")+ theme(axis.text=element_text(size=10),axis.title=element_text(size=11), plot.title = element_text(size=13), legend.title = element_text(size = 10),  legend.text = element_text(size = 10),legend.key=element_blank(), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), aspect.ratio = 1) + geom_point(data=volc_all[c(scalop::Signatures_GBM$NPC1,scalop::Signatures_GBM$NPC2,scalop::Signatures_GBM$OPC),], fill="blue",shape=21,size=2)+ geom_point(data=volc_all[c(scalop::Signatures_GBM$MES1,scalop::Signatures_GBM$MES2,scalop::Signatures_GBM$AC),], fill="red",shape=21,size=2)+ geom_text_repel(data=volc_all[labels,],label=labels,size=3.5,min.segment.length = 0.1) #+geom_point(data=volc_all[volc_all$log2FoldChange>1 & volc_all$pvalue<0.05,], fill="lightgrey",shape=21,size=2.5)
dev.off()


reactive_astrocyte_sig<- read.csv("reactive_astrocyte_sig.csv")

reactive_astrocyte_sig$MCAO_Reactive_astrocytes<- mouse_to_human(reactive_astrocyte_sig$MCAO_Reactive_astrocytes)$human

reactive_astrocyte_sig$LPS_Reactive_astrocytes<- mouse_to_human(reactive_astrocyte_sig$LPS_Reactive_astrocytes)$human


ast_list<-list("MCAO"=reactive_astrocyte_sig$MCAO_Reactive_astrocytes[!is.na(reactive_astrocyte_sig$MCAO_Reactive_astrocytes)],"LPS"=reactive_astrocyte_sig$LPS_Reactive_astrocytes[!is.na(reactive_astrocyte_sig$LPS_Reactive_astrocytes)])


reactive_ast_scores<- sigScores(data.frame(merged_tpm_matrix[,rownames(merged_annotation_df)[merged_annotation_df$CNA=="Non-malignant"]]), ast_list)


melt_reactive_ast_scores<- melt(as.matrix(reactive_ast_scores))

melt_reactive_ast_scores$Cell_type<- merged_annotation_df[melt_reactive_ast_scores$Var1,]$Cell_type
melt_reactive_ast_scores$Type<-  merged_annotation_df[melt_reactive_ast_scores$Var1,]$Type

melt_reactive_ast_scores$Cell_type<- factor(melt_reactive_ast_scores$Cell_type, levels=c( "aRG"  , "oRG.Astroglia" ,"Cycling.oRG.Astroglia",  "Immature.interneurons", "IP", "CPN","CFuPN" ))

melt_reactive_ast_scores<- melt_reactive_ast_scores[!is.na(melt_reactive_ast_scores$Cell_type),]

v1<-ggboxplot(melt_reactive_ast_scores, x = "Var2", y = "value",
              facet.by = "Cell_type", fill = "Type",
              palette = c("darkgrey", "chartreuse"), add.params = list(size = 0.3)) +
  stat_compare_means(aes(group = Type), label.y = 1.8,label = "p.signif", method = "t.test",
                     method.args = list(alternative = "greater", adjust.method = "bonferroni"), hide.ns = T)+xlab("")


pdf("reactive_astrocyte_boxplotV2.pdf",width = 8, height = 10,onefile=F)
v1+xlab("Reactive astrocyte sig")+ylab("Score")
dev.off()


v1<-ggboxplot(melt_reactive_ast_scores, x = "Var2", y = "value",
              facet.by = "Cell_type", fill  = "Type",
              palette = c("darkgrey", "chartreuse"), add.params = list(size = 0.3)) +
  stat_compare_means(aes(group = Type), label.y = 1.8,label = "p.signif", method = "t.test",
                     method.args = list(alternative = "greater", adjust.method = "bonferroni"), hide.ns = T)+xlab("")


pdf("reactive_astrocyte_boxplot_celltype.pdf",width = 8, height = 6,onefile=F)
v1+xlab("Reactive astrocyte sig")+ylab("Score")
dev.off()



##### #####



