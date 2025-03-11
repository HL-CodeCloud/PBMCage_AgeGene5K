# 
# ### Plot heatmap of the age-correlated genes in all the cells regardless of celltypes
# #####################################
# ###
# 
# ### Load the pseudobulk by donor
# pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
# pseudobulkobj=AggregateExpression(pseudobulk_obj, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
# pseudobulkobj$age=as.numeric(unlist(lapply(strsplit(colnames(pseudobulkobj), split="_"), function(x) x[1])))
# pseudobulkobj$sex=unlist(lapply(strsplit(colnames(pseudobulkobj), split="_"), function(x) x[2]))
# pseudobulkobj$donor_id=unlist(lapply(strsplit(colnames(pseudobulkobj), split="_"), function(x) x[3]))
# 
# ### Load the genes
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_donorid.csv.gz", sep="\t")
# pos_genes=Cor_analysis_DF %>% subset(rho_pval<0.05 & analysis=="all" & rho>0) %>% 
#   slice_max(rho, n=15) %>% select(gene) %>% unlist() %>% unname()d
# neg_genes=Cor_analysis_DF %>% subset(rho_pval<0.05 & analysis=="all" & rho<0) %>% 
#   slice_min(rho, n=15) %>% select(gene) %>% unlist() %>% unname()
# genes_selected=c(pos_genes, neg_genes)
# # check the list with DEMAGALHAES_AGING_UP.v2023.2.Hs.gmt
# AgingDB=clusterProfiler::read.gmt("~/Project_PBMCage/DEMAGALHAES_AGING_UP.v2023.2.Hs.gmt")[["gene"]]
# AgingDB[match(genes_selected, AgingDB)]
# 
# ### Make the matrix
# gene_expr=Seurat::FetchData(pseudobulkobj, vars=c("donor_id","sex","age",genes_selected), layer="data")
# gene_expr=gene_expr %>% group_by(age) %>% summarize_at(colnames(.)[!grepl("donor_id|age|sex",colnames(.))], mean)
# # convert to matrix
# gene_expr_mtx=gene_expr %>% column_to_rownames("age") %>% data.table::transpose(.)
# rownames(gene_expr_mtx)=colnames(gene_expr)[!grepl("age", colnames(gene_expr))]
# colnames(gene_expr_mtx)=gene_expr$age
# # scale
# gene_expr_mtx=t(scale(t(gene_expr_mtx)))
# 
# ### Visualize the age effects
# library(ComplexHeatmap)
# p=
#   Heatmap(gene_expr_mtx,
#           name=" ",
#           show_heatmap_legend=FALSE, col=c("dodgerblue3","white","brown3"),
#           cluster_columns=F,
#           row_names_gp=gpar(fontsize=9),
#           column_names_side="top",
#           column_names_rot=90,
#           column_names_gp=gpar(fontsize=9), 
#           height=0.1*nrow(gene_expr_mtx), 
#           width=0.15*ncol(gene_expr_mtx))
# p=draw(p)
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Heatmap.pdf", height=0.15*nrow(gene_expr_mtx), width=0.15*ncol(gene_expr_mtx))
# plot(p)
# dev.off()
# 
# #####################################
# 
# 
# 
# ### Analyze the sig genes found across ages in all cells
# #####################################
# ###
# 
# ### Enrichment analysis
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_donorid.csv.gz", sep="\t")
# sig_genes=Cor_analysis_DF %>% subset(rho_pval<0.05 & analysis=="all") %>% 
#   mutate(orientation=ifelse(rho>0, "upreg.", "downreg.")) %>%
#   select(gene, orientation)
# sig_genes_list=split(sig_genes$gene, sig_genes$orientation)
# 
# library(clusterProfiler)
# ego_all=compareCluster(sig_genes_list, 
#                        fun=enrichGO,
#                        keyType="SYMBOL",
#                        OrgDb="org.Hs.eg.db",
#                        ont="BP",
#                        pvalueCutoff=0.05)
# 
# plot_=
#   clusterProfiler::dotplot(ego_all, showCategory=10,
#                            label_format=50) +
#   scale_colour_gradient(limits=quantile(ego_all@compareClusterResult$p.adjust)[c(1,10)],
#                         low="brown4", high="pink",
#                         labels= ~sprintf(fmt="%0.01e", .)) +
#   scale_size_continuous(range=c(2,6)) +
#   labs(x=NULL) +
#   theme(axis.text.x=element_text(size=10),
#         axis.text.y=element_text(size=10),
#         title=element_text(size=11)) +
#   scale_x_discrete(expand=expansion(mult=c(0, 0), add=c(0.25, 0.25)))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Enrichment_age.cor.genes.pdf", width=5, height=6)
# plot_
# dev.off()
# 
# #####################################
# 
# 
# 
# ### Score the T cell recombination genes
# #####################################
# ### 
# 
# ### Load
# THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
# CD8T=subset(THEObj, Annot.inter=="CD8T.naive")
# # add the agecut metadata
# AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
# AGECUT_name=c()
# for (i in 1:length(AGECUT)) {
#   AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
# }
# agecut_df=CD8T[[]]
# agecut_df$age=as.factor(agecut_df$age)
# agecut_df$agecut=agecut_df$age
# for (i in 1:length(AGECUT)) {
#   agecut_df=agecut_df %>%
#     mutate(agecut=ifelse(age %in% AGECUT[[i]], AGECUT_name[i], agecut))
# }
# table(rownames(agecut_df)==colnames(CD8T)) # check
# CD8T=AddMetaData(CD8T, metadata=agecut_df$agecut, col.name="agecut")
# 
# ### Make a function for plotting on ages
# geneset_scoring_forAge=function(dataframe, celltype_as_LegendName, scale) {
#   df_full_list=list()
#   CD8T_NK_subsets_celltypes=unique(dataframe$Annot.detailed)
#   for (i in 1:length(CD8T_NK_subsets_celltypes)) {
#     matrix_=subset(dataframe, Annot.detailed==CD8T_NK_subsets_celltypes[i])
#     age_not_there=unique(dataframe$age)[!(unique(dataframe$age) %in% unique(matrix_$age))]
#     if (length(age_not_there)!=0) {
#       matrix_foragenotthere=data.frame(Annot.detailed=unique(matrix_$Annot.detailed),
#                                        Annot.inter=unique(matrix_$Annot.inter),
#                                        Annot.rough=unique(matrix_$Annot.rough),
#                                        age=age_not_there,
#                                        BCL11B=rep(NA, length(age_not_there)),
#                                        LEF1=rep(NA, length(age_not_there)),
#                                        LIG4=rep(NA, length(age_not_there)),
#                                        PRKDC=rep(NA, length(age_not_there)),
#                                        TCF7=rep(NA, length(age_not_there)))
#     } else {
#       matrix_foragenotthere=data.frame()
#     }
#     matrix_=data.table::rbindlist(list(matrix_, matrix_foragenotthere))
#     matrix_=matrix_ %>% arrange(age)
#     age_=matrix_$age
#     matrix_=t(matrix_[,c("BCL11B", "LEF1", "LIG4", "PRKDC", "TCF7")])
#     matrix_=as.data.frame(matrix_)
#     colnames(matrix_)=age_
#     rownames(matrix_)=c("BCL11B", "LEF1", "LIG4", "PRKDC", "TCF7")
#     df_full_list[[i]]=matrix_
#   }
#   df_full_combined=data.table::rbindlist(df_full_list)
#   
#   if (scale==T) {
#     df_full_combined=t(scale(t(df_full_combined)))
#     quantiles=quantile(df_full_combined, na.rm=TRUE)
#     col_fun=circlize::colorRamp2(c(-2,0,2), c("#00007E", "white", "#7E0000"))}
#   if (scale==F) {
#     df_full_combined=df_full_combined-rowMeans(df_full_combined, na.rm=T)
#     quantiles=quantile(df_full_combined, na.rm=TRUE)
#     col_fun=circlize::colorRamp2(c(-0.5,0,0.5), c("#00007E", "white", "#7E0000"))}
#   
#   p_age=list()
#   codes=""
#   for (i in 1:length(df_full_list)) {
#     df_to_draw=as.data.frame(df_full_list[[i]])
#     if (scale==T) {df_to_draw=t(scale(t(df_to_draw)))} else {df_to_draw=df_to_draw-rowMeans(df_to_draw, na.rm=T)}
#     
#     if (i==1) {
#       p_age[[i]]=
#         Heatmap(df_to_draw,
#                 na_col="gray99",
#                 col=col_fun,
#                 name=paste0("Mean Expr. of\n", celltype_as_LegendName),
#                 heatmap_legend_param=list(title_position="leftcenter-rot",
#                                           title_gp=gpar(fontsize=8, fontface="bold")),
#                 column_title=NULL,
#                 row_title=CD8T_NK_subsets_celltypes[i],
#                 row_title_gp=gpar(fontsize=8),
#                 row_title_rot=0,
#                 show_column_names=T,
#                 cluster_columns=F,
#                 cluster_rows=F,
#                 show_row_dend=F,
#                 row_names_gp=gpar(fontsize=8),
#                 column_names_gp=gpar(fontsize=8),
#                 column_names_side="top",
#                 width=ncol(df_to_draw)*unit(3,"mm"), 
#                 height=nrow(df_to_draw)*unit(3,"mm"))}
#     if (i>=2) {
#       p_age[[i]]=
#         Heatmap(df_to_draw,
#                 na_col="gray99",
#                 col=col_fun,
#                 show_heatmap_legend=FALSE,
#                 column_title=NULL,
#                 row_title=CD8T_NK_subsets_celltypes[i],
#                 row_title_gp=gpar(fontsize=8),
#                 row_title_rot=0,
#                 show_column_names=F,
#                 cluster_columns=F,
#                 cluster_rows=F,
#                 show_row_dend=F,
#                 row_names_gp=gpar(fontsize=8),
#                 column_names_gp=gpar(fontsize=8),
#                 column_names_side="top",
#                 width=ncol(df_to_draw)*unit(3,"mm"), 
#                 height=nrow(df_to_draw)*unit(3,"mm"))
#     }
#     codes=paste0(codes, " %v% p_age[[", i, "]]")
#   }
#   codes=gsub("^ \\%v\\% ","",codes)
#   eval(parse(text=paste0("ht_list=", codes)))
#   p_age_total=draw(ht_list, ht_gap=unit(1, "mm"))
#   
#   return(p_age_total)
# }
# 
# ### Score the recombination genes in CD8 T cell subtypes
# Idents(CD8T)="age"
# Recombination_genes=FetchData(CD8T, 
#                               vars=c("Annot.detailed","Annot.inter","Annot.rough",
#                                      "age","agecut","sex",
#                                      "BCL11B", "LEF1", "LIG4", "PRKDC", "TCF7"), 
#                               layer="data")
# Recombination_genes_group=Recombination_genes %>%
#   group_by(Annot.detailed, Annot.inter, Annot.rough, age) %>%
#   summarize(BCL11B=mean(BCL11B),
#             LEF1=mean(LEF1),
#             LIG4=mean(LIG4),
#             PRKDC=mean(PRKDC),
#             TCF7=mean(TCF7)) %>%
#   as.data.frame()
# CD8T_recombination_plots=geneset_scoring_forAge(dataframe=Recombination_genes_group, celltype_as_LegendName="CD8T.naive", scale=F)
# 
# ### Plot
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_VDJrecombination_CD8T.pdf", width=12, height=2.5)
# draw(CD8T_recombination_plots)
# dev.off()
# 
# #####################################
# 
# 
# 
# ### Map hub genes from WGCNA modules at the rough level to marker genes of celltypes at various Annot levels
# #####################################
# ### 
# 
# library(Seurat)
# library(hdWGCNA)
# library(dplyr)
# library(ggplot2)
# library(clusterProfiler)
# 
# ### Get the hub genes
# WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Rough.rds")
# hub_genes=GetHubGenes(WGCNAobj, n_hubs=length(WGCNAobj@misc$WGCNA$wgcna_genes)) %>% 
#   subset(module!="grey") %>%
#   group_by(module) %>% 
#   arrange(desc(kME), .by_group=T) %>% as.data.frame() %>% 
#   split(.$module) %>%
#   lapply(., function(df) df %>% select(gene_name, kME) %>% tibble::deframe())
# hub_genes=hub_genes[names(hub_genes)!="grey"]
# 
# ### Load the marker genes
# merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels.csv.gz", sep="\t")
# merged_results=subset(merged_results, Annot.inter!="Other.Eryth" & Annot.inter!="Other.Mast")
# 
# ### Map to rough level
# merged_results_rough=merged_results %>% subset(celltype.level=="Annot.rough") %>% 
#   select(Annot.rough, gene)
# gse_rough=hub_genes %>% 
#   lapply(., function(x) GSEA(x, scoreType="pos", TERM2GENE=merged_results_rough, pvalueCutoff=1))
# 
# gse_rough
# 
# ### Map to inter level
# merged_results_inter=merged_results %>% subset(celltype.level=="Annot.inter") %>% 
#   select(Annot.inter, gene)
# gse_inter=hub_genes %>% 
#   lapply(., function(x) GSEA(x, scoreType="pos", TERM2GENE=merged_results_inter, pvalueCutoff=1))
# 
# gse_inter
# 
# gse_inter_names=names(gse_inter)
# gse_inter_resultDF=lapply(1:length(gse_inter), 
#                           function(idx) gse_inter[[idx]]@result %>% mutate(module=paste0("M",idx))) %>%
#   data.table::rbindlist(.) %>% 
#   select(ID, NES, module) %>% 
#   tidyr::pivot_wider(names_from="module", values_from="NES") %>%
#   tibble::column_to_rownames("ID")
# gse_inter_resultDF=gse_inter_resultDF[rownames(gse_inter_resultDF)[order(rownames(gse_inter_resultDF))],]
# gse_inter_resultDF[is.na(gse_inter_resultDF)]=0
# 
# gse_inter_resultDF_pval=
#   lapply(1:length(gse_inter), 
#          function(idx) gse_inter[[idx]]@result %>% mutate(module=paste0("M",idx))) %>%
#   data.table::rbindlist(.) %>% 
#   select(ID, p.adjust, module) %>% 
#   tidyr::pivot_wider(names_from="module", values_from="p.adjust") %>%
#   tibble::column_to_rownames("ID")
# gse_inter_resultDF_pval=apply(gse_inter_resultDF_pval, c(1,2), function(x) ifelse(x<0.0001, "****", 
#                                                                                   ifelse(x>0.0001 & x<0.001, "***",
#                                                                                          ifelse(x>0.001 & x<0.01, "**", 
#                                                                                                 ifelse(x>0.01 & x<0.05, "*", NA)))))
# gse_inter_resultDF_pval=gse_inter_resultDF_pval[rownames(gse_inter_resultDF_pval)[order(rownames(gse_inter_resultDF_pval))],]
# gse_inter_resultDF_pval[is.na(gse_inter_resultDF_pval)]=""
# 
# 
# library(ComplexHeatmap)
# col_fun=circlize::colorRamp2(c(0, max(gse_inter_resultDF, na.rm=T)), c("white", "brown3"))
# Heatmap(gse_inter_resultDF, col=col_fun, cluster_rows=F, cluster_columns=F)
# ht=
#   Heatmap(gse_inter_resultDF, 
#           name="mat", show_heatmap_legend=F,
#           col=col_fun,
#           cluster_rows=F, cluster_columns=F,
#           column_names_rot=0, column_names_centered=TRUE, column_names_side="top",
#           column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
#           height=5*nrow(gse_inter_resultDF), 
#           width=5*ncol(gse_inter_resultDF),
#           
#           layer_fun=function(j, i, x, y, width, height, fill) {
#             grid.text(as.matrix(gse_inter_resultDF_pval), x, y, 
#                       gp=gpar(fontsize=10))}
#           
#   )
# lgd=Legend(col_fun=col_fun, 
#            title=NULL,
#            legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
#            direction="horizontal", title_position="lefttop", labels_gp=gpar(fontsize=9),
# )
# ht_final=draw(ht, annotation_legend_list=lgd, annotation_legend_side="bottom")
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation_rough_ggplot.pdf", width=4, height=1)
# ht_final
# dev.off()
# 
# 
# ### Score the enrichment terms from aPEAR analysis on WGCNA modules at the inter level
# #####################################
# ### 
# 
# library(clusterProfiler)
# library(hdWGCNA)
# library(aPEAR)
# library(dplyr)
# library(ggplot2)
# 
# ### Extract the modules
# WGCNA_obj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# mt_cor=GetModuleTraitCorrelation(WGCNA_obj)
# Module_genes=GetModules(WGCNA_obj) %>% subset(module!="grey")
# Module_genes$module=as.character(Module_genes$module)
# genes_per_module=split(Module_genes$gene_name, Module_genes$module)
# 
# ### aPEAR enrichment
# EGO=aPEAR_cluster=list()
# for (i in 1:length(genes_per_module)) {
#   genes_=genes_per_module[i][[1]]
#   EGO[[i]]=enrichGO(gene=genes_,
#                     OrgDb="org.Hs.eg.db",
#                     keyType="SYMBOL",
#                     ont="BP",
#                     pvalueCutoff=0.05,
#                     readable=TRUE)
#   result=gofilter(EGO[[i]], level=3)@result
#   message(paste0("EGO ",i," of ", length(genes_per_module), " has been enriched."))
#   
#   aPEAR_cluster[[i]]=findPathClusters(result)$clusters
# }
# 
# ### Get the genes from the aPEAR terms of my interest
# TOP_Genes_per_module=list()
# for (i in 1:length(aPEAR_cluster)) {
#   the_apear_result=aPEAR_cluster[[i]]
#   # take the cluster >3
#   top_cluster=the_apear_result %>% count(Cluster) %>% subset(n>=3) %>% select(Cluster) %>% unlist() %>% unname()
#   top_pathway=the_apear_result %>% subset(Cluster %in% top_cluster)
#   top_pathway=split(top_pathway$Pathway, top_pathway$Cluster)
#   
#   # get the genes enriched in the top pathways
#   TOP_GENEs_per_pathway=list()
#   TOP_GENEs_per_pathway_NAME=c()
#   for (j in 1:length(top_pathway)) {
#     the_ego_result=gofilter(EGO[[i]], level=3)@result
#     top_genes=the_ego_result %>% subset(Description %in% top_pathway[[j]]) %>% select(geneID) %>% unlist() %>% unname()
#     top_genes=unlist(unname(sapply(top_genes, function(x) strsplit(x, "/")[[1]])))
#     top_genes=top_genes[!duplicated(top_genes)]
#     
#     TOP_GENEs_per_pathway[[j]]=top_genes
#     TOP_GENEs_per_pathway_NAME[j]=names(top_pathway)[j]
#   }
#   names(TOP_GENEs_per_pathway)=TOP_GENEs_per_pathway_NAME
#   
#   TOP_Genes_per_module[[i]]=TOP_GENEs_per_pathway
# }
# 
# ### Score the terms
# library(Seurat)
# library(UCell)
# the_object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
# SCORE_DFs=list()
# for (i in 1:length(TOP_Genes_per_module)) {
#   the_module=TOP_Genes_per_module[[i]]
#   # score
#   Ucell_scored_obj=AddModuleScore_UCell(the_object, features=the_module)
#   # extract the ucell scores
#   score_df=Ucell_scored_obj[[]] %>% select(-any_of(c("orig.ident")))
#   score_df=score_df %>% group_by(Annot.inter, age, sex) %>% summarize_at(vars(colnames(.)[grepl("_UCell", colnames(.))]), mean)
#   # save the score df
#   SCORE_DFs[[i]]=score_df
# }
# names(SCORE_DFs)=paste0("M",c(1:5))
# saveRDS(SCORE_DFs, "~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter_Module_aPEAR.top.terms_genes_ExprUcell.rds")
# 
# ### Plot the heatmap
# library(ComplexHeatmap)
# # select the terms of interest
# terms_list=list(M1=list(`immune system development`=c("immune system development","T cell selection")),
#                 M2=list(`immune defense`=c("viral budding","regulation of viral process","positive regulation by symbiont of entry into host"),
#                         `hemostasis`=c("regulation of hemostasis")),
#                 M3=list(`immune defense`=c("regulation of viral process","viral transcription","viral budding"),
#                         `development and reproduction`=c("developmental pigmentation","meiotic cell cycle process","ovulation cycle","multi-multicellular organism process")),
#                 M4=list(`cellular pigmentation`=c("cellular pigmentation"),
#                         `biosynthesis`=c("production of molecular mediator of immune response"),
#                         `aging`=c("anatomical structure regression","muscle adaptation","aging")),
#                 M5=list(`development and reproduction`=c("post-embryonic development","maintenance of cell number","ovulation cycle","embryo implantation"),
#                         `hemostasis`=c("regulation of hemostasis")))
# # process the data
# SCORE_DFs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter_Module_aPEAR.top.terms_genes_ExprUcell.rds")
# 
# ht=list()
# for (i in 1:length(SCORE_DFs)) {
#   SCORE_DFs_merge=SCORE_DFs[[i]]
#   # regardless of sex
#   SCORE_DFs_merge=SCORE_DFs_merge %>% group_by(age, Annot.inter) %>% summarize_at(vars(colnames(.)[grepl("_UCell", colnames(.))]), mean)
#   # filter the terms of interest
#   score_sum=list()
#   for (j in 1:length(terms_list[[i]])) {
#     SCORE_DFs_merge_subset=SCORE_DFs_merge %>% select(any_of(c("Annot.inter","age",paste0(terms_list[[i]][[j]], "_UCell"))))
#     score_sum[[j]]=apply(SCORE_DFs_merge_subset[,!grepl("Annot\\.inter|age",colnames(SCORE_DFs_merge_subset))], 1, sum)
#   }
#   names(score_sum)=names(terms_list[[i]])
#   SCORE_DFs_merge_overallterm=cbind(SCORE_DFs_merge[,c("age","Annot.inter")], as.data.frame(score_sum))
#   colnames(SCORE_DFs_merge_overallterm)=c("age","Annot.inter",names(score_sum))
#   
#   # calculate the pearson cor
#   pearson_results=SCORE_DFs_merge_overallterm %>% split(~Annot.inter) %>% lapply(., function(df) {apply(df[,c(1,3:ncol(df))], 2, function(x) cor.test(df$age, x, method="pearson")$estimate)})
#   pearson_df=as.data.frame(pearson_results) %>% .[2:nrow(.),]
#   rownames(pearson_df)=gsub("_UCell","",rownames(pearson_df))
#   pearson_df.t=data.table::transpose(pearson_df)
#   colnames(pearson_df.t)=rownames(pearson_df); rownames(pearson_df.t)=colnames(pearson_df)
#   pearson_df.t=apply(pearson_df.t, c(1,2), function(x) ifelse(is.na(x), 0, x))
#   
#   col_fun=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
#   ht[[i]]=
#     Heatmap(pearson_df.t,
#             name="mat", show_heatmap_legend=F,
#             col=col_fun,
#             cluster_rows=T, cluster_columns=T, show_column_dend=F,
#             column_names_rot=45, column_names_side="bottom", column_names_centered=F,
#             column_labels=stringr::str_wrap(colnames(pearson_df.t), width=15),
#             column_names_gp=gpar(fontsize=10.5), 
#             row_names_gp=gpar(fontsize=10),
#             width=unit(10, "mm")*ncol(pearson_df.t), 
#             height=unit(4.2, "mm")*nrow(pearson_df.t))
# }
# ht_list=ht[[1]]+ht[[2]]+ht[[3]]+ht[[4]]+ht[[5]]
# 
# lgd=Legend(col_fun=col_fun, 
#            title="Pearson correlation",
#            legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
#            direction="horizontal", title_position="topcenter", title_gp=gpar(font=1, fontsize=11), labels_gp=gpar(fontsize=9),
# )
# ht_final=draw(ht_list, annotation_legend_list=lgd, annotation_legend_side="top")
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleAPEAR_inter_EnrichedGeneExprHeatmap.pdf", height=7, width=6.5)
# draw(ht_final)
# dev.off()
# 
# #####################################
# 
# 
# 
# ### Analyze the sex effect on the enrichment terms from aPEAR analysis on WGCNA modules at the inter level
# #####################################
# ### 
# 
# library(dplyr)
# library(ggplot2)
# 
# ### Load the scoring data of the enrichment terms
# SCORE_DFs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter_Module_aPEAR.top.terms_genes_ExprUcell.rds")
# # select the terms of interest
# terms_list=list(M1=list(`immune system development`=c("immune system development","T cell selection")),
#                 M2=list(`immune defense`=c("viral budding","regulation of viral process","positive regulation by symbiont of entry into host"),
#                         `hemostasis`=c("regulation of hemostasis")),
#                 M3=list(`immune defense`=c("regulation of viral process","viral transcription","viral budding"),
#                         `development and reproduction`=c("developmental pigmentation","meiotic cell cycle process","ovulation cycle","multi-multicellular organism process")),
#                 M4=list(`cellular pigmentation`=c("cellular pigmentation"),
#                         `biosynthesis`=c("production of molecular mediator of immune response"),
#                         `aging`=c("anatomical structure regression","muscle adaptation","aging")),
#                 M5=list(`development and reproduction`=c("post-embryonic development","maintenance of cell number","ovulation cycle","embryo implantation"),
#                         `hemostasis`=c("regulation of hemostasis")))
# 
# PEARSON_DFs=list()
# for (i in 1:length(SCORE_DFs)) {
#   SCORE_DFs_merge=SCORE_DFs[[i]]
#   pearson_df.t_SexCombined=list()
#   
#   # for female
#   SCORE_DFs_merge_F=SCORE_DFs_merge %>% subset(sex=="female") %>% group_by(age, Annot.inter) %>% summarize_at(vars(colnames(.)[grepl("_UCell", colnames(.))]), mean)
#   # filter the terms of interest
#   score_sum=list()
#   for (j in 1:length(terms_list[[i]])) {
#     SCORE_DFs_merge_subset=SCORE_DFs_merge_F %>% select(any_of(c("Annot.inter","age",paste0(terms_list[[i]][[j]], "_UCell"))))
#     score_sum[[j]]=apply(SCORE_DFs_merge_subset[,!grepl("Annot\\.inter|age",colnames(SCORE_DFs_merge_subset))], 1, sum)
#   }
#   names(score_sum)=names(terms_list[[i]])
#   SCORE_DFs_merge_overallterm=cbind(SCORE_DFs_merge_F[,c("age","Annot.inter")], as.data.frame(score_sum))
#   colnames(SCORE_DFs_merge_overallterm)=c("age","Annot.inter",names(score_sum))
#   # pearson_results=SCORE_DFs_merge_overallterm %>% split(~Annot.inter) %>% lapply(., function(df) {apply(df[,c(1,3:ncol(df))], 2, function(x) cor.test(df$age, x, method="pearson")$estimate)})
#   pearson_results=SCORE_DFs_merge_overallterm %>% split(~Annot.inter) %>% lapply(., function(df) {apply(df[,c(1,3:ncol(df))], 2, function(x) abs(lm(x~df$age)$coefficients[2]))})
#   pearson_df=as.data.frame(pearson_results) %>% .[2:nrow(.),]
#   rownames(pearson_df)=gsub("_UCell","",rownames(pearson_df))
#   pearson_df.t=data.table::transpose(pearson_df)
#   colnames(pearson_df.t)=rownames(pearson_df); rownames(pearson_df.t)=colnames(pearson_df)
#   pearson_df.t=apply(pearson_df.t, c(1,2), function(x) ifelse(is.na(x), 0, x))
#   pearson_df.t=as.data.frame(pearson_df.t); pearson_df.t$sex="F"
#   
#   pearson_df.t_SexCombined[["F"]]=pearson_df.t %>% mutate(celltype=rownames(.)) %>% 
#     tidyr::pivot_longer(cols=colnames(.)[!grepl("sex|celltype",colnames(.))], names_to="term", values_to="cor")
#   
#   
#   # for male
#   SCORE_DFs_merge_M=SCORE_DFs_merge %>% subset(sex=="male") %>% group_by(age, Annot.inter) %>% summarize_at(vars(colnames(.)[grepl("_UCell", colnames(.))]), mean)
#   # filter the terms of interest
#   score_sum=list()
#   for (j in 1:length(terms_list[[i]])) {
#     SCORE_DFs_merge_subset=SCORE_DFs_merge_M %>% select(any_of(c("Annot.inter","age",paste0(terms_list[[i]][[j]], "_UCell"))))
#     score_sum[[j]]=apply(SCORE_DFs_merge_subset[,!grepl("Annot\\.inter|age",colnames(SCORE_DFs_merge_subset))], 1, sum)
#   }
#   names(score_sum)=names(terms_list[[i]])
#   SCORE_DFs_merge_overallterm=cbind(SCORE_DFs_merge_M[,c("age","Annot.inter")], as.data.frame(score_sum))
#   colnames(SCORE_DFs_merge_overallterm)=c("age","Annot.inter",names(score_sum))
#   pearson_results=SCORE_DFs_merge_overallterm %>% split(~Annot.inter) %>% lapply(., function(df) {apply(df[,c(1,3:ncol(df))], 2, function(x) abs(lm(x~df$age)$coefficients[2]))})
#   pearson_df=as.data.frame(pearson_results) %>% .[2:nrow(.),]
#   rownames(pearson_df)=gsub("_UCell","",rownames(pearson_df))
#   pearson_df.t=data.table::transpose(pearson_df)
#   colnames(pearson_df.t)=rownames(pearson_df); rownames(pearson_df.t)=colnames(pearson_df)
#   pearson_df.t=apply(pearson_df.t, c(1,2), function(x) ifelse(is.na(x), 0, x))
#   pearson_df.t=as.data.frame(pearson_df.t); pearson_df.t$sex="M"
#   
#   pearson_df.t_SexCombined[["M"]]=pearson_df.t %>% mutate(celltype=rownames(.)) %>% 
#     tidyr::pivot_longer(cols=colnames(.)[!grepl("sex|celltype",colnames(.))], names_to="term", values_to="cor")
#   
#   
#   # calculate the diff between the two genders
#   pearson_df.t_Sex=data.table::rbindlist(pearson_df.t_SexCombined) %>% tidyr::pivot_wider(names_from="sex", values_from="cor") %>%
#     mutate(diff=`F`-`M`)
#   
#   # save to the list
#   PEARSON_DFs[[i]]=pearson_df.t_Sex %>% mutate(module=paste0("M",i))
# }
# 
# ### Combine the results from all the modules
# PEARSON_DFs_combo=data.table::rbindlist(PEARSON_DFs)
# 
# ### Plot
# plot_=
#   ggplot(PEARSON_DFs_combo, aes(x=celltype, y=diff, color=term, fill=term)) +
#   facet_wrap(~module, nrow=1) +
#   geom_bar(stat="identity", position=position_dodge2(preserve="single", width=1), alpha=0.3, width=0.8) +
#   theme_minimal() +
#   labs(x=NULL, y="Sex difference",title=" \n\n\n") +
#   theme(axis.text.x=element_text(size=18, angle=90, hjust=1, vjust=0.5),
#         axis.text.y=element_text(size=16),
#         axis.title.y=element_text(size=20),
#         title=element_text(size=24),
#         panel.grid.minor=element_blank(),
#         # panel.grid.major.x=element_blank(),
#         legend.title=element_blank(),
#         legend.text=element_text(size=20),
#         strip.text=element_text(size=20)) +
#   scale_fill_manual(values=paletteer::paletteer_d("ggsci::category20c_d3")) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20c_d3"))
# 
# # get legend
# legend_p=cowplot::get_legend(plot_)
# # remove legend
# plot_updated=plot_+
#   theme(legend.position="none")
# 
# # cowplot
# plot_final_=
#   cowplot::plot_grid(cowplot::plot_grid(NULL, NULL, NULL, NULL, legend_p, axis="t", nrow=1, align="hv", ncol=5,
#                                         rel_widths=c(1,1,1,1,1)),
#                      cowplot::plot_grid(plot_updated),
#                      axis="l", align="v", nrow=2, rel_heights=c(1,4))
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleAPEAR_inter_EnrichedGeneExpr_FvsM.pdf", height=14, width=40)
# plot(plot_final_)
# dev.off()
# 
# #####################################




### Plot general gene expr in B and T naive cells
#####################################
###

### Load the seurat obj
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all" & grepl("CD4T\\.|CD8T\\.|B\\.",celltypes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% summarize_at(c("rho","rho_pval"), mean))
cor_data_subset=lapply(1:length(cor_data_subset), function(idx) cor_data_subset[[idx]] %>%
                         mutate(Annot.inter=names(cor_data_subset)[idx])) %>%
  data.table::rbindlist()

ggplot(cor_data_subset, aes(x=cor, y=-log10(cor_pval))) +
  ggrastr::geom_point_rast() +
  facet_wrap(~celltypes)










