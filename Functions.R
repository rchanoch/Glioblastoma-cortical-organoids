to_cpm <- function(m) {
  m <- as.matrix(m)
  count_sum <- apply(m, 2, sum)
  m <- (t(t(m)/count_sum)) * 1e+06
  m <- as(m, "Matrix")
  return(m)
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

unlog_matrix<- function(matrix=matrix)
{
  matrix<- (2^(matrix)-1)*10
  return(matrix)
}

get_top_genes<- function(matrix=matrix, isLog=F, cutoff=4)
{
  if(isLog)
    matrix<- (2^(matrix)-1)*10
  
  return(rownames(matrix)[log2(1+rowMeans(matrix))>=cutoff])
  
}



count_to_tpm<- function(matrix)
{
  matrix<- data.frame(matrix)
  count_sum <- apply(matrix, 2, sum) # sums all UMI counts for each sample
  matrix <- t(t(matrix)/count_sum) # divides all genes in each sample by the total sum
  matrix <- matrix*1000000 # "per million"
  return(matrix)
}



create_state_piechart<- function(cell_type_list= cell_type_list, order_by=NULL){
  
  # calculate the start and end angles for each pie
  dat<-  get_state_percentages(cell_type_list)
  colnames(dat)<- c("Cell_Type","Cnt","Channel")
  
  dat_pies <- dplyr::left_join(dat,
                               dat %>% 
                                 group_by(Channel) %>%
                                 dplyr::summarize(Cnt_total = sum(Cnt))) %>%
    group_by(Channel) %>%
    mutate(end_angle = 2*pi*cumsum(Cnt)/Cnt_total,      # ending angle for each pie slice
           start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
           mid_angle = 0.5*(start_angle + end_angle))   # middle of each pie slice, for the text label
  
  rpie = 1 # pie radius
  rlabel = 0.6 * rpie # radius of the labels; a number slightly larger than 0.5 seems to work better,
  # but 0.5 would place it exactly in the middle as the question asks for.
  
  # draw the pies
  
  dat_pies$Cell_Type<- as.factor(dat_pies$Cell_Type)
  
  data <- dat_pies %>% 
    dplyr::arrange(desc(Cell_Type)) %>%
    dplyr::mutate(prop = Cnt) %>%
    dplyr::mutate(ypos = cumsum(prop)- 0.5*prop )
  
  
  ggplot(data, aes(x="", y=Cnt, fill=Cell_Type)) +
    geom_bar(stat="identity", width=1,color="black") +
    coord_polar("y", start=0)+facet_wrap( ~ Channel, ncol = length(unique(dat_pies$Channel)))+theme_void()+ theme(strip.text = element_text(size = 10))+geom_text(aes(y = ypos, label = Cnt), color = "black", size=4)
  
}

split_names<- function(string, position, char)
{
  return(unlist(lapply(strsplit(string,char,1),"[[",position)))
}

volcano_plot<- function(centered_tpm, cluster_cells=cluster_cells)
{
  t_test <- apply(centered_tpm, 1, function(x){t.test(x[is.element(colnames(centered_tpm), cluster_cells)], x[!is.element(colnames(centered_tpm), cluster_cells)])$p.value})
  
  fc <- rowMeans(centered_tpm[,is.element(colnames(centered_tpm), cluster_cells)]) - rowMeans(centered_tpm[,!is.element(colnames(centered_tpm), cluster_cells)])
  
  volc <- data.frame(cbind(t_test, fc))
  volc$Gene<- rownames(volc)
  colnames(volc)<- c("pvalue","log2FoldChange")
  volc<- volc[order(volc$log2FoldChange, decreasing=T),]
  
  return(volc)
  
}

plot_volc<- function(volc, title=NULL, upgenes=NULL,downgenes=NULL)
{
  
  print(ggplot(volc, aes(x=log2FoldChange, y=-log10(pvalue)))+ labs(x="log2FC", y="-log10(p value)") + geom_point(alpha = 0.6, size=2,color="grey") + ggtitle("")+ theme(axis.text=element_text(size=20),axis.title=element_text(size=20), plot.title = element_text(size=24), legend.text = element_text(size=16), legend.key=element_blank(), panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), aspect.ratio = 1) + geom_point(data=volc[downgenes,], colour="red")+ geom_point(data=volc[upgenes,], colour="red")+ ggrepel::geom_text_repel(data=volc[c(upgenes,downgenes),], label=c(upgenes,downgenes))+geom_vline(xintercept = 0, linetype="dashed")+ ggtitle(title))
  
}


basic_heatmap_wrapper<- function(matrix_hc, limits=3, showAllGenes=F, genes_list=NULL, showAllCells=FALSE,legend="Relative\nExpression",xlab="Cells",ylab="")
{
  plot_test <- melt(as.matrix(matrix_hc))
  
  break_points= signif(ncol(matrix_hc),1)/5
  
  genes= c(genes_list)
  if(showAllGenes)
    genes<- c(genes, rownames(matrix_hc))
  
  
  if(showAllCells)
    ggplot(plot_test, aes(Var2, Var1, fill = value))+ geom_tile()+ ggtitle("")+scale_fill_gradient2(low="steelblue", mid="white", high="darkred", oob=squish, midpoint=0, limit = c(-limits,limits), space = "Lab", name=legend,guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + theme(axis.text.y = element_text(size=7),panel.border = element_rect(colour = "black", fill=NA, size=1) ,  axis.ticks = element_blank(), axis.line = element_blank(), axis.text.x = element_text(hjust=1, angle=45),legend.title = element_text(size = 8),  legend.text = element_text(size = 8),legend.key.size = unit(.8,"line"))+ labs(x = xlab, y= ylab)+ scale_y_discrete(breaks=genes)+ scale_x_discrete(breaks=colnames(matrix_hc))
  else
    ggplot(plot_test, aes(Var2, Var1, fill = value))+ geom_tile()+ ggtitle("")+scale_fill_gradient2(low="steelblue", mid="white", high="darkred", oob=squish, midpoint=0, limit = c(-limits,limits), space = "Lab", name=legend,guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + theme(axis.text.y = element_text(size=7),panel.border = element_rect(colour = "black", fill=NA, size=1) ,  axis.ticks = element_blank(), axis.line = element_blank(),legend.title = element_text(size = 8),  legend.text = element_text(size = 8),legend.key.size = unit(.8,"line"))+ labs(x = xlab, y= ylab)+ scale_y_discrete(breaks=genes)+scale_x_discrete(breaks=colnames(matrix_hc)[seq(break_points, ncol(matrix_hc), by=break_points)], labels= seq(break_points, ncol(matrix_hc), by=break_points))
  
  
  
}


NMF_heatmap<- function(NMF_programs=NMF_programs,log_expression_matrix=log_expression_matrix, order_by=NULL,num_genes=15,showAllGenes1=TRUE,genes_list=NULL, limits=3, maxcol=0.3,center=F)
{
  cell_type_list=list()
  NMF_programs<- lapply(NMF_programs, function(x){x[x %in% rownames(log_expression_matrix)]})
  NMF_program_scores<- sigScores(log_expression_matrix, NMF_programs,expr.center = center)
  iden<- maxcol_strict(NMF_program_scores, min =maxcol)
  names<- c(names(NMF_programs))
  
  if(maxcol!=0)
  {
    no_id<- rownames(NMF_program_scores)[!rownames(NMF_program_scores) %in% iden]
    names(no_id)<- rep("unclassified", length(no_id))
    iden<- c(iden,no_id)
    names<- c(names(NMF_programs), "unclassified")
    
  }
  
  
  cells_order<-""
  
  
  for(cell_type in names)
  {
    temp_cells<- iden[names(iden)==cell_type]
    if(cell_type != "unclassified")
    {
      temp_cells<- temp_cells[order(NMF_program_scores[temp_cells,paste(cell_type)], decreasing=T)]
      
    }
    cells_order<- c(cells_order, temp_cells)
    cell_type_list[[cell_type]]<- temp_cells
  }
  
  cells_order<- c(cells_order, cells_order[!cells_order %in% colnames(log_expression_matrix)])
  
  if(length(order_by)>0)
  {
    new_order<-""
    for(i in 1:length(order_by))
    {
      temp_order1<- cells_order[grep(order_by[i], cells_order)]
      new_order<- c(new_order,temp_order1)
    }
    cells_order<- new_order[-1]
  }
  
  marker_programs<- unlist(lapply(NMF_programs, function(x){x[1:num_genes]}))
  cells_order<- cells_order[cells_order %in% colnames(log_expression_matrix)]
  
  if(center==FALSE) centered_log_tpm <- log_expression_matrix-rowMeans(log_expression_matrix)
  else centered_log_tpm <- log_expression_matrix
  
  matrix_hc <- as.matrix(centered_log_tpm[marker_programs,cells_order])
  
  #top_marker_programs<- unlist(lapply(NMF_programs, function(x){x[1:5]}))
  
  #heatmap_NMF<- basic_heatmap_wrapper(matrix_hc, genes_list = c(unlist(Signatures_GBM),top_marker_programs))
  heatmap_NMF<- basic_heatmap_wrapper(matrix_hc,showAllGenes = showAllGenes1, limits = limits, genes_list = genes_list)
  
  return_list<- list("heatmap"=heatmap_NMF,"cell_type_list"=cell_type_list,"cells_order"=cells_order)
  return(return_list)
}



mouse_to_human<- function(mouse_genes=mouse_genes)
{
  
  human_to_mouse<-read.table("HMD_HumanPhenotype.rpt.txt",sep="\t",header=F)
  
  human_to_mouse<- human_to_mouse[,c("V1","V3")]
  colnames(human_to_mouse)<- c("human","mouse")
  
  human_to_mouse<-human_to_mouse[match(mouse_genes, human_to_mouse$mouse),]
  rownames(human_to_mouse)<- mouse_genes
  
  return(human_to_mouse[mouse_genes,c("mouse","human")])
  
  
  
}


human_to_mouse<- function(human_genes=human_genes)
{
  
  human_to_mouse<-read.table("HMD_HumanPhenotype.rpt.txt",sep="\t",header=F)
  human_to_mouse<- human_to_mouse[,c("V1","V3")]
  colnames(human_to_mouse)<- c("human","mouse")
  
  
  human_to_mouse<-human_to_mouse[match(human_genes, human_to_mouse$human),]
  rownames(human_to_mouse)<- human_genes
  
  return(human_to_mouse[human_genes,c("mouse","human")])
  
  
  
}

