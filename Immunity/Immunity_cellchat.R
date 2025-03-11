gc()
library(CellChat)
library(tidyverse)
library(Seurat)
library(dplyr)

# ##########################################################################
# ### Analyze on the celltypes at the rough level
# ##########################################################################
# 
# 
# ### Prepare the cellchat objs for the analysis
# #####################################
# ###
# 
# object=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
# object=Seurat::NormalizeData(object)
# 
# ### Prepare Agecut
# object[[]]=object[[]] %>%
#   dplyr::select(-agecut) %>%
#   dplyr::rename(age=Age) %>%
#   mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
#   mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(41:44), "41~44", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
#   dplyr::rename(agecut=agecut.half)
# 
# ### Divide the obj by agecut and take only assigned cells
# AGECUT_name=names(table(object$agecut))
# object_list=list()
# for (i in 1:length(AGECUT_name)) {
#   object_list[[i]]=subset(object, agecut==AGECUT_name[i])
# }
# names(object_list)=AGECUT_name
# 
# ### Create the cellchat objects
# cellchat=list()
# for (i in 1:length(object_list)) {
#   DefaultAssay(object_list[[i]])="RNA"
#   cellchat[[i]]=createCellChat(object=object_list[[i]], meta=object_list[[i]]@meta.data, group.by="Annot.rough")
# }
# 
# ### load the cellchat database
# CellChatDB=CellChatDB.human
# for (i in 1:length(cellchat)) {
#   cellchat[[i]]@DB=CellChatDB
# }
# gc()
# saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough.rds")
# 
# ### Identify overexpression
# processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
# for (i in 1:length(cellchat)) {
#   cellchat[[i]]=subsetData(cellchat[[i]])
#   cellchat[[i]]=identifyOverExpressedGenes(cellchat[[i]])
#   cellchat[[i]]=identifyOverExpressedInteractions(cellchat[[i]])
#   setTxtProgressBar(processbar, i)
# }
# close(processbar)
# gc()
# 
# ### Evaluate the communications
# options(future.globals.maxSize=1*1024^3) # first try computeCommunProb, if error then increase the memory: 850MiB * 1024^2 = 891289600
# processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
# for (i in 1:length(cellchat)) {
#   cellchat[[i]]=computeCommunProb(cellchat[[i]], raw.use=TRUE, population.size=TRUE)
#   cellchat[[i]]=filterCommunication(cellchat[[i]], min.cells=10)
#   setTxtProgressBar(processbar, i)
# }
# close(processbar)
# gc()
# saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough.rds")
# 
# ### Calculate the communication networks
# processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
# for (i in 1:length(cellchat)) {
#   cellchat[[i]]=computeCommunProbPathway(cellchat[[i]])
#   cellchat[[i]]=aggregateNet(cellchat[[i]])
#   setTxtProgressBar(processbar, i)
# }
# close(processbar)
# gc()
# saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough.rds")
# 
# ### Calculate the communication networks
# for (i in 1:length(cellchat)) {
#   cellchat[[i]]=netAnalysis_computeCentrality(cellchat[[i]], slot.name="netP")
# }
# names(cellchat)=AGECUT_name
# gc()
# saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough.rds")
# 
# #####################################
# # 
# # 
# # 
# # ### Visualize the chats among every pair of celltypes in each cellchat objs
# # #####################################
# # ### 
# # 
# # ###
# # cellchat=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_rough.rds")
# # AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
# # AGECUT_name=c()
# # for (i in 1:length(AGECUT)) {
# #   AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
# # }
# # 
# # ### Plot the netVisual_circle (communications from source cells to receiver cells)
# # pdf("~/Project_PBMCage/Plots/Cellchat_perobj_netVisualcircle.pdf", width=10, height=60)
# # par(mfrow=c(length(AGECUT),2), xpd=TRUE, mar=c(1,1,1,1))
# # for (i in 1:length(cellchat)) {
# #   netVisual_circle(cellchat[[i]]@net$count,
# #                    vertex.weight=as.numeric(table(cellchat[[i]]@idents)), 
# #                    weight.scale=T,
# #                    label.edge=F, 
# #                    title.name=AGECUT_name[i],
# #                    color.use=scater:::.get_palette("tableau10medium")[1:length(levels(cellchat[[i]]@idents))],
# #                    alpha.edge=0.6)
# #   netVisual_circle(cellchat[[i]]@net$weight,
# #                    vertex.weight=as.numeric(table(cellchat[[i]]@idents)), 
# #                    weight.scale=T,
# #                    label.edge=F, 
# #                    title.name=AGECUT_name[i],
# #                    color.use=scater:::.get_palette("tableau10medium")[1:length(levels(cellchat[[i]]@idents))],
# #                    alpha.edge=0.6)
# # }
# # dev.off()
# # 
# # ### Plot the netVisual_heatmap (communications from source cells to receiver cells)
# # source("~/Rscripts/Cellchat_functions.R")
# # library(ComplexHeatmap)
# # # determine the color scale
# # mat_count=mat_weight=list()
# # for (i in 1:length(cellchat)) {
# #   mat_count[[i]]=netVisual_heatmap2_getMAT(cellchat[[i]], measure=c("count"))
# #   mat_weight[[i]]=netVisual_heatmap2_getMAT(cellchat[[i]], measure=c("weight"))
# # }
# # mat_count_max=max(sapply(mat_count, max))
# # mat_count_min=min(sapply(mat_count, min))
# # mat_weight_max=max(sapply(mat_weight, max))
# # mat_weight_min=min(sapply(mat_weight, min))
# # print(list(mat_count_max, mat_count_min, mat_weight_max, mat_weight_min))
# # # plot; remember to modify the color scale!
# # plot_ht_count=plot_ht_weight=list()
# # for (i in 1:length(cellchat)) {
# #   plot_ht_count[[i]]=
# #     netVisual_heatmap2(cellchat[[i]], color.heatmap="Reds", measure="count", title.name=paste0("Number of interactions\n[",AGECUT_name[i],"]"),
# #                       color.heatmap.use_forcomb=circlize::colorRamp2(c(0, 25), c("#FFF5F0", "#67000D")))
# #   plot_ht_weight[[i]]=
# #     netVisual_heatmap2(cellchat[[i]], color.heatmap="Reds", measure="weight", title.name=paste0("Interaction strength\n[",AGECUT_name[i],"]"),
# #                       color.heatmap.use_forcomb=circlize::colorRamp2(c(0, 0.1), c("#FFF5F0", "#67000D")))
# # }
# # codes1=codes2=""
# # for (i in 1:length(cellchat)) {
# #   codes1=paste0(codes1, " + plot_ht_count[[",i,"]]")
# #   codes2=paste0(codes2, " + plot_ht_weight[[",i,"]]")
# # }
# # codes1=gsub("^ \\+ ", "", codes1)
# # codes2=gsub("^ \\+ ", "", codes2)
# # 
# # lgd_count=Legend(col_fun=circlize::colorRamp2(c(0, 25), c("#FFF5F0", "#67000D")),
# #                  title="Number of interactions",
# #                  title_gp=gpar(fontsize=8, fontface="plain"),
# #                  title_position="leftcenter-rot",
# #                  border=NA,
# #                  legend_height=unit(20, "mm"),
# #                  labels_gp=gpar(fontsize=8),
# #                  grid_width=unit(2, "mm"))
# # eval(parse(text=paste0("plot_ht_count_all=draw(",codes1,", ht_gap=unit(0.5, 'cm'), annotation_legend_list=lgd_count)")))
# # 
# # lgd_weight=Legend(col_fun=circlize::colorRamp2(c(0, 0.1), c("#FFF5F0", "#67000D")),
# #                   title="Number of interactions",
# #                   title_gp=gpar(fontsize=8, fontface="plain"),
# #                   title_position="leftcenter-rot",
# #                   border=NA,
# #                   legend_height=unit(20, "mm"),
# #                   labels_gp=gpar(fontsize=8),
# #                   grid_width=unit(2, "mm"))
# # eval(parse(text=paste0("plot_ht_weight_all=draw(",codes2,", ht_gap=unit(0.5, 'cm'), annotation_legend_list=lgd_weight)")))
# # 
# # pdf("~/Project_PBMCage/Plots/Cellchat_perobj_netVisualheatmap.pdf", width=20, height=3)
# # plot(plot_ht_count_all)
# # plot(plot_ht_weight_all)
# # dev.off()
# # 
# # ### Plot the netVisual_aggregate (communications from source cells to receiver cells)
# # levels(cellchat[[1]]@idents)
# #   # list all the celltypes and use vertex.receiver to set the target celltypes in the left plot
# #   # the remaining celltypes will be taken as the target celltypes in the right plot
# # vertex.receiver=c(1,2,3,8)
# # plot_aggr=list()
# # for (i in 1:length(cellchat)) {
# #   plot_aggr_each=
# #     netVisual_aggregate(cellchat[[i]], signaling=cellchat[[i]]@netP$pathways,  
# #                         vertex.receiver=vertex.receiver,
# #                         color.use=scater:::.get_palette("tableau10medium")[1:length(levels(cellchat[[i]]@idents))],
# #                         layout="hierarchy",
# #                         alpha.edge=0.8
# #                         )
# #   plot_aggr_each[[1]][[60]][[2]][[2]]=rep(NA, length(plot_aggr_each[[1]][[60]][[2]][[2]])) # remove the title
# #   plot_aggr_each[[1]][[60]][[2]][[2]][1]=AGECUT_name[i]
# #   plot_aggr[[i]]=cowplot::ggdraw(plot_aggr_each)
# # }
# # pdf("~/Project_PBMCage/Plots/Cellchat_perobj_netVisualaggregate.pdf", height=60, width=8)
# # cowplot::plot_grid(plotlist=plot_aggr, align="hv", nrow=length(AGECUT))
# # dev.off()
# # 
# # ### Plot the netAnalysis_signalingRole_heatmap (signaling pathways in each celltype)
# # pathway.union=list()
# # for (i in 1:(length(cellchat)-1)) {
# #   pathway.union[[i]]=union(cellchat[[i]]@netP$pathways, cellchat[[i+1]]@netP$pathways)
# # }
# # pathway.union=unlist(pathway.union)
# # pathway.union=pathway.union[!duplicated(pathway.union)]
# # 
# # library(ComplexHeatmap)
# # plot_sig=list()
# # for (i in 1:length(cellchat)) {
# #   plot_sig_all=
# #     netAnalysis_signalingRole_heatmap(cellchat[[i]],
# #                                       pattern="all", 
# #                                       signaling=pathway.union,
# #                                       title=names(cellchat)[i],
# #                                       width=8, height=10)
# #   plot_sig_out=
# #     netAnalysis_signalingRole_heatmap(cellchat[[i]],
# #                                       pattern="outgoing", 
# #                                       signaling=pathway.union,
# #                                       title=names(cellchat)[i],
# #                                       width=8, height=10)
# #   plot_sig_in=
# #     netAnalysis_signalingRole_heatmap(cellchat[[i]],
# #                                       pattern="incoming", 
# #                                       signaling=pathway.union,
# #                                       title=names(cellchat)[i],
# #                                       width=8, height=10)
# #   
# #   plot_sig[[i]]=draw(plot_sig_all + plot_sig_out + plot_sig_in, ht_gap=unit(0.5, "cm"))
# # }
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cellchat_perobj_netAnalysissignalingRoleheatmap.pdf", width=15, height=7)
# # for (i in 1:length(plot_sig)) plot(plot_sig[[i]])
# # dev.off()
# # 
# # ### Plot the netVisual_bubble (ligand-receptor pairs expr on source cells to receiver cells)
# # plot_bubble=list()
# # for (i in 1:length(cellchat)) {
# #   plot_bubble[[i]]=
# #     netVisual_bubble(cellchat[[i]],
# #                      sources.use=vertex.receiver,
# #                      targets.use=levels(cellchat[[i]]@idents)[sapply(1:length(levels(cellchat[[i]]@idents)), function(x) !(x %in% vertex.receiver))], 
# #                      remove.isolate=FALSE,
# #                      title.name=AGECUT_name[i]
# #                      )
# # }
# # pdf("~/Project_PBMCage/Plots/Cellchat_perobj_netVisualbubble.pdf")
# # for (i in 1:length(cellchat)) {
# #   plot(plot_bubble[[i]])
# # }
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### Merge all the cellchat objs and plot the interaction comparison
# # #####################################
# # ### 
# # ###
# # 
# # ### Merge the cellchat objs all in one
# # cellchat=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_rough.rds")
# # AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
# # AGECUT_name=c()
# # for (i in 1:length(AGECUT)) {
# #   AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
# # }
# # names(cellchat)=AGECUT_name
# # cellchat.list=mergeCellChat(cellchat, add.names=names(cellchat), cell.prefix=TRUE)
# # plot_count=compareInteractions(cellchat.list, show.legend=F, measure="count")
# # plot_weight=compareInteractions(cellchat.list, show.legend=F, measure="weight")
# # plot_interactions=plot_count+plot_weight
# # 
# # pdf("~/Project_PBMCage/Plots/Cellchat_compareall_compareInteractions.pdf", width=9, height=7)
# # plot(plot_interactions)
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### Merge young and old cellchat objs for comparison
# # #####################################
# # ### 
# # ### (i.e., AGECUT 1-7: "19~24" "25~30" "31~42" "43~51" "52~58" "59~64" "65~73", vs., AGECUT 8-11: "74~81" "82~88" "89~91" "92~97")
# # 
# # library(Seurat)
# # object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
# # 
# # ### Prepare the cellchat objs for the analysis
# # # Add Agecut information from the determined agecuts by "WGCNA_withoutSampleFilter.R"
# # agecut_df=object[[]]
# # agecut_df$age=as.factor(agecut_df$age)
# # agecut_df$agecut=agecut_df$age
# # AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
# # AGECUT_name=c()
# # for (i in 1:length(AGECUT)) {
# #   AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
# # }
# # for (i in 1:length(AGECUT)) {
# #   agecut_df=agecut_df %>%
# #     mutate(agecut=ifelse(age %in% AGECUT[[i]], AGECUT_name[i], agecut))
# # }
# # agecut_df=agecut_df %>%
# #   mutate(young_old=ifelse(agecut %in% AGECUT_name[1:7], "young", 
# #                           ifelse(agecut %in% AGECUT_name[8:11], "old", NA)))
# # table(rownames(agecut_df)==colnames(object)) # check
# # object=AddMetaData(object, metadata=agecut_df$young_old, col.name="young_old")
# # 
# # ### Divide the obj by young_old and take only assigned cells
# # AGECUTs=c("young","old")
# # object_list=list()
# # for (i in 1:length(AGECUTs)) {
# #   object_list[[i]]=subset(object, young_old==AGECUTs[i])
# # }
# # names(object_list)=AGECUTs
# # 
# # ### Create the cellchat objects
# # cellchat=list()
# # for (i in 1:length(object_list)) {
# #   DefaultAssay(object_list[[i]])="RNA"
# #   cellchat[[i]]=createCellChat(object=object_list[[i]], meta=object_list[[i]]@meta.data, group.by="Annot.rough")
# # }
# # 
# # ### load the cellchat database
# # CellChatDB=CellChatDB.human
# # for (i in 1:length(cellchat)) {
# #   cellchat[[i]]@DB=CellChatDB
# # }
# # saveRDS(cellchat, "~/Project_PBMCage/Tempt_RDS/Cellchat_rough_youngAndold.rds")
# # 
# # ### Identify overexpression
# # processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
# # for (i in 1:length(cellchat)) {
# #   cellchat[[i]]=subsetData(cellchat[[i]])
# #   cellchat[[i]]=identifyOverExpressedGenes(cellchat[[i]])
# #   cellchat[[i]]=identifyOverExpressedInteractions(cellchat[[i]])
# #   setTxtProgressBar(processbar, i)
# # }
# # close(processbar)
# # saveRDS(cellchat, "~/Project_PBMCage/Tempt_RDS/Cellchat_rough_youngAndold.rds")
# # 
# # ### Evaluate the communications
# # options(future.globals.maxSize=1*1024^3) # first try computeCommunProb, if error then increase the memory: 850MiB * 1024^2 = 891289600
# # processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
# # for (i in 1:length(cellchat)) {
# #   cellchat[[i]]=computeCommunProb(cellchat[[i]], raw.use=TRUE, population.size=TRUE)
# #   cellchat[[i]]=filterCommunication(cellchat[[i]], min.cells=10)
# #   setTxtProgressBar(processbar, i)
# # }
# # close(processbar)
# # saveRDS(cellchat, "~/Project_PBMCage/Tempt_RDS/Cellchat_rough_youngAndold.rds")
# # 
# # ### Calculate the communication networks
# # processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
# # for (i in 1:length(cellchat)) {
# #   cellchat[[i]]=computeCommunProbPathway(cellchat[[i]])
# #   cellchat[[i]]=aggregateNet(cellchat[[i]])
# #   setTxtProgressBar(processbar, i)
# # }
# # close(processbar)
# # saveRDS(cellchat, "~/Project_PBMCage/Tempt_RDS/Cellchat_rough_youngAndold.rds")
# # 
# # ### Calculate the communication networks
# # for (i in 1:length(cellchat)) {
# #   cellchat[[i]]=netAnalysis_computeCentrality(cellchat[[i]], slot.name="netP")
# # }
# # saveRDS(cellchat, "~/Project_PBMCage/Tempt_RDS/Cellchat_rough_youngAndold.rds")
# # 
# # ### Merge the cellchat objs
# # cellchat=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_rough_youngAndold.rds")
# # names(cellchat)=c("young","old")
# # cellchat.list=mergeCellChat(cellchat, add.names=c("young","old"), cell.prefix=TRUE)
# # plot_count=compareInteractions(cellchat.list, show.legend=F, measure="count")
# # plot_weight=compareInteractions(cellchat.list, show.legend=F, measure="weight")
# # plot_interactions=plot_count+plot_weight
# # 
# # ### Compute the similarity between young and old
# # cellchat.list=computeNetSimilarityPairwise(cellchat.list, type="functional")
# # cellchat.list=netEmbedding(cellchat.list, umap.method="uwot", type="functional")
# # cellchat.list=netClustering(cellchat.list, do.parallel=FALSE, type="functional")
# # # netVisual_embeddingPairwise(cellchat.list, type="functional", label.size=3.5)
# # # netVisual_embeddingPairwiseZoomIn(cellchat.list, type="functional", nCol=2)
# # # cellchat.list=computeNetSimilarityPairwise(cellchat.list, type="structural")
# # # cellchat.list=netEmbedding(cellchat.list, umap.method="uwot", type="structural")
# # # cellchat.list=netClustering(cellchat.list, do.parallel=FALSE, type="structural")
# # # netVisual_embeddingPairwise(cellchat.list, type="structural", label.size=3.5)
# # # netVisual_embeddingPairwiseZoomIn(cellchat.list, type="structural", nCol=2)
# # saveRDS(cellchat.list, "~/Project_PBMCage/Tempt_RDS/Cellchat_rough_youngAndold_mergedObj.rds")
# # 
# # ### Plot the netVisual_circle (compare the communications from source cells to receiver cells between young and old)
# # # trace(netVisual_diffInteraction, edit=TRUE)
# # plot_compar_count=
# #   netVisual_diffInteraction(cellchat.list, weight.scale=TRUE, measure="count", title.name="age<74 vs. age>=74")
# # plot_compar_weight=
# #   netVisual_diffInteraction(cellchat.list, weight.scale=TRUE, measure="weight", title.name="age<74 vs. age>=74")
# # 
# # ### Plot the rankSimilarity (compare the functional similarity of pathways between each two comparison pairs)
# # plot_ranksim=rankSimilarity(cellchat.list, type="functional") + 
# #   ggtitle("Functional similarity of pathway [age<74 vs. age>=74]")
# # # rankSimilarity(cellchat.list, type="structural") + ggtitle("Structural similarity of pathway")
# # 
# # ### Plot the signaling comparison between each compare_pair by merging their pathways
# # library(ComplexHeatmap)
# # ht_all=ht_out=ht_in=list()
# # cellchat1=cellchat[[1]]
# # cellchat2=cellchat[[2]]
# # pathway.union=union(cellchat1@netP$pathways, cellchat2@netP$pathways)
# # ht_all1=netAnalysis_signalingRole_heatmap(cellchat1, pattern="all", signaling=pathway.union, title="age<74", width=8, height=10)
# # ht_all2=netAnalysis_signalingRole_heatmap(cellchat2, pattern="all", signaling=pathway.union, title="age>=74", width=8, height=10)
# # ht_out1=netAnalysis_signalingRole_heatmap(cellchat1, pattern="outgoing", signaling=pathway.union, title="age<74", width=8, height=10)
# # ht_out2=netAnalysis_signalingRole_heatmap(cellchat2, pattern="outgoing", signaling=pathway.union, title="age>=74", width=8, height=10)
# # ht_in1=netAnalysis_signalingRole_heatmap(cellchat1, pattern="incoming", signaling=pathway.union, title="age<74", width=8, height=10)
# # ht_in2=netAnalysis_signalingRole_heatmap(cellchat2, pattern="incoming", signaling=pathway.union, title="age>=74", width=8, height=10)
# # ht_all=draw(ht_all1 + ht_all2, ht_gap=unit(0.5, "cm"))
# # ht_out=draw(ht_out1 + ht_out2, ht_gap=unit(0.5, "cm"))
# # ht_in=draw(ht_in1 + ht_in2, ht_gap=unit(0.5, "cm"))
# # 
# # ### Plot the rankNet (compare the signaling pathways between young and old and do statistics)
# # plot_ranknet=
# #   rankNet(cellchat.list, mode="comparison", stacked=T, do.stat=TRUE)
# # 
# # ### Plot the netVisual_bubble (compare the expr of ligand-receptor pairs on source cells to receiver cells between young and old)
# # levels(cellchat[[1]]@idents)
# # # list all the celltypes and use vertex.receiver to set the target celltypes in the left plot
# # # the remaining celltypes will be taken as the target celltypes in the right plot
# # vertex.receiver=c(1,2,3,8)
# # plot_compare_bubble=
# #   netVisual_bubble(cellchat.list, 
# #                    comparison=c(1,2),
# #                    sources.use=vertex.receiver,
# #                    targets.use=c(1:length(levels(cellchat.list@idents$joint)))[sapply(1:length(levels(cellchat.list@idents$joint)), function(x) !(x %in% vertex.receiver))], 
# #                    remove.isolate=FALSE,
# #                    angle.x=45,
# #                    title.name="age<74 vs. age>=74")
# # 
# # ### Plot all the figs above
# # pdf("~/Project_PBMCage/Plots/Cellchat_compareYoungAndOld_allplots.pdf", width=10, height=7)
# # plot(plot_interactions)
# # plot_compar_weight
# # plot(plot_ranksim)
# # plot(ht_all)
# # plot(ht_out)
# # plot(ht_in)
# # plot(plot_ranknet)
# # plot(plot_compare_bubble)
# # dev.off()
# # 
# # #####################################
# 
# 
# 
### Extract the results
#####################################
###

library(CellChat)
library(dplyr)

cellchat=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough.rds")
agecuts=names(cellchat)

### Count and Weight
cellchat.count=cellchat.weight=list()
for (i in 1:length(cellchat)) {
  cellchat.count[[i]]=cellchat[[i]]@net$count
  cellchat.count[[i]]=as.data.frame(cellchat.count[[i]]) %>%
    tibble::rownames_to_column(var="Source") %>%
    tidyr::pivot_longer(!Source, names_to="Target", values_to="Count")
  cellchat.count[[i]]$Agecut=agecuts[i]

  cellchat.weight[[i]]=cellchat[[i]]@net$weight
  cellchat.weight[[i]]=as.data.frame(cellchat.weight[[i]]) %>%
    tibble::rownames_to_column(var="Source") %>%
    tidyr::pivot_longer(!Source, names_to="Target", values_to="Weight")
  cellchat.weight[[i]]$Agecut=agecuts[i]
}
CELLCHAT.COUNT=data.table::rbindlist(cellchat.count)
CELLCHAT.WEIGHT=data.table::rbindlist(cellchat.weight)

### Signals
cellchat.signal=list()
for (i in 1:length(cellchat)) {
  tempt=cellchat[[i]]@netP$prob
  cellchat.signal_j=list()
  for (j in 1:dim(tempt)[3]) {
    cellchat.signal_j[[j]]=as.data.frame(tempt[,,j]) %>%
      tibble::rownames_to_column(var="Source") %>%
      tidyr::pivot_longer(!Source, names_to="Target", values_to="Prob")
    cellchat.signal_j[[j]]$Signal=dimnames(tempt)[[3]][j]
  }
  cellchat.signal[[i]]=data.table::rbindlist(cellchat.signal_j)
  cellchat.signal[[i]]$Agecut=agecuts[i]
}
CELLCHAT.SIGNAL=data.table::rbindlist(cellchat.signal)

### Pairs
cellchat.pair=list()
for (i in 1:length(cellchat)) {
  tempt=cellchat[[i]]@net$prob
  pair_list_names=names(tempt[1,1,])
  pair_list=lapply(1:dim(tempt)[3], function(idx) sum(tempt[,,idx])) %>% unlist(); names(pair_list)=pair_list_names
  pair_list=pair_list[pair_list!=0]
  match_pair_signal=cellchat[[i]]@LR$LRsig %>% .[pair_list_names,] %>% dplyr::rename(Pair=interaction_name)
  
  cellchat.pair_j=list()
  for (j in 1:length(pair_list)) {
    cellchat.pair_j[[j]]=as.data.frame(tempt[,,names(pair_list)[j]]) %>%
      tibble::rownames_to_column(var="Source") %>%
      tidyr::pivot_longer(!Source, names_to="Target", values_to="Prob")
    cellchat.pair_j[[j]]$Pair=names(pair_list)[j]
  }
  cellchat.pair[[i]]=data.table::rbindlist(cellchat.pair_j)
  cellchat.pair[[i]]$Agecut=agecuts[i]
  cellchat.pair[[i]]=cellchat.pair[[i]] %>% left_join(match_pair_signal, by=c("Pair"))
}
CELLCHAT.PAIR=data.table::rbindlist(cellchat.pair)

saveRDS(list(CELLCHAT.COUNT=CELLCHAT.COUNT, CELLCHAT.WEIGHT=CELLCHAT.WEIGHT, CELLCHAT.SIGNAL=CELLCHAT.SIGNAL, CELLCHAT.PAIR=CELLCHAT.PAIR),
        "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")

#####################################
# 
# 
# 
# ### Plot the flow-in/out among rough cell types, regardless of signals
# #####################################
# ###
# 
# library(ggplot2)
# library(dplyr)
# 
# CELLCHAT.WEIGHT=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.WEIGHT"]]
# 
# ### Plot
# library(ggalluvial)
# celltype_flow.in.out=
#   ggplot(CELLCHAT.WEIGHT,
#          aes(y=Weight, axis1=Source, axis2=Target)) +
#   geom_alluvium(aes(fill=Source),
#                 width=3/8, knot.pos=0, reverse=FALSE) +
#   geom_stratum(color="grey50", fill="white", linewidth=0.1, alpha=.25, width=3/8, reverse=FALSE) +
#   scale_fill_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
#   guides(fill="none") +
#   geom_text(stat="stratum", aes(label=after_stat(stratum)), reverse=FALSE, size=3.5, angle=90) +
#   scale_x_continuous(breaks=1:2, labels=c("Source","Target")) +
#   coord_flip() +
#   theme_void() +
#   theme(legend.position="bottom",
#         plot.background=element_rect(color="transparent", fill="white"),
#         axis.text.y=element_text(size=11, angle=90, hjust=0.5))
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_rough_celltype_alluvial.pdf", height=3.5, width=8)
# celltype_flow.in.out
# dev.off()
# 
# #####################################
# # 
# # axis_x_data=data.frame(x=c(unique(CELLCHAT.SIGNAL_NK_source$Agecut), "999"),
# #                        y=rep(0, length(unique(CELLCHAT.SIGNAL_NK_source$Agecut))+1))
# # axis_plot=
# #   ggplot(axis_x_data, aes(x=x, y=y)) +
# #   geom_point(data=. %>% subset(x!="999"), size=2, shape=1, color="black") +
# #   geom_segment(aes(x=unique(CELLCHAT.SIGNAL_NK_source$Agecut)[1], 
# #                    xend="999", 
# #                    y=0, yend=0),
# #                arrow=arrow(length=unit(0.3, "cm")),
# #                color="grey50", linewidth=unit(0.1, "cm")) +
# #   theme_void() +
# #   theme(axis.text.x=element_text(size=9)) +
# #   coord_cartesian(ylim=c(-0.1,0.1)) +
# #   scale_x_discrete(labels=c(unique(CELLCHAT.SIGNAL_NK_source$Agecut), " "))
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cellchat_x_axis_bar.pdf", height=0.32, width=10)
# # plot(axis_plot)
# # dev.off()
# # 
# #   
# #   
# # 
# # 
# # 
# ### Plot the interaction counts/weights among cell types across age
# #####################################
# ###
# 
# ###
# library(ggplot2)
# library(dplyr)
# 
# CELLCHAT.COUNT=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.COUNT"]]
# CELLCHAT.WEIGHT=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.WEIGHT"]]
# 
# # # plot counts
# # count_plot=
# #   ggplot(CELLCHAT.COUNT, aes(x=Agecut, y=Count, color=Target)) +
# #   geom_point(size=0.5) +
# #   geom_line(aes(group=Target, color=Target), linetype="dashed", linewidth=0.5) +
# #   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
# #   facet_grid(cols=vars(Source))
# 
# # plot weight
# weight_plot_cut=
#   ggplot(CELLCHAT.WEIGHT %>% subset(!((Source=="CD4T cells" | Source=="CD8T cells") & Target=="CD8T cells")),
#          aes(x=Agecut, y=Weight, color=Target)) +
#   facet_grid(cols=vars(Source)) +
#   geom_point(size=0.5, alpha=0.25) +
#   geom_line(aes(group=Target), linetype="solid", linewidth=0.25) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
#   labs(title="Source", y="Interaction strength", x=NULL) +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=9),
#         strip.text=element_text(size=10),
#         # axis.text.x=element_text(size=7.5, angle=90, hjust=0, vjust=0.5),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2, alpha=1)))
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_Weight.pdf", width=10, height=3)
# plot(weight_plot_cut)
# dev.off()
# 
# # plot weight from CD4T/CD8T to CD8T
# weight_plot_CD4CD8=
#   ggplot(CELLCHAT.WEIGHT %>% subset((Source=="CD4T cells" | Source=="CD8T cells") & Target=="CD8T cells"),
#          aes(x=Agecut, y=Weight, color=Target)) +
#   facet_grid(cols=vars(Source)) +
#   geom_point(size=0.5) +
#   geom_line(aes(group=Target), linetype="dashed", linewidth=0.5) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[3]) +
#   labs(title="Signals to CD8T cells", y="Interaction strength", x=NULL) +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=9),
#         strip.text=element_text(size=10),
#         axis.text.x=element_text(size=7.5, angle=90, hjust=0, vjust=0.5),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color="none")
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_Weight_CD4CD8T.pdf", width=4, height=3.5)
# plot(weight_plot_CD4CD8)
# dev.off()
# 
# #####################################
# 
# 
# 
# ### Plot the signals in and out NK cells
# #####################################
# ###
# 
# library(ggplot2)
# library(dplyr)
# 
# CELLCHAT.SIGNAL=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.SIGNAL"]]
# 
# # take the NK cells only
# CELLCHAT.SIGNAL_NK_source=CELLCHAT.SIGNAL %>% 
#   subset((Source=="NK cells" | Target=="NK cells") & !grepl("DCs|Other cells",Target) & !grepl("DCs|Other cells",Source))
# 
# # mark the top signals
# marking_target=CELLCHAT.SIGNAL_NK_source %>%
#   subset(Prob!=0) %>%
#   group_by(Source, Target, Signal) %>%
#   summarize_at("Prob", mean) %>%
#   slice_max(Prob, n=3) %>%
#   mutate(marker=Signal) %>%
#   dplyr::select(-Prob)
# 
# plot_per_signal_NK_target=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Source=="NK cells"), aes(x=Agecut, y=Prob, color=Target)) +
#   facet_wrap(~Target, nrow=1) +
#   geom_point(size=0.5, alpha=0.25) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.6) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(1,2,3,5,6,8)]) +
#   labs(y="Interaction strength", x=NULL, title="Signals from NK cells", subtitle="Target") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# plot_per_signal_NK_source=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Target=="NK cells"), aes(x=Agecut, y=Prob)) +
#   facet_wrap(~Source, nrow=1) +
#   geom_point(size=0.5, alpha=0.25, 
#              color=paletteer::paletteer_d("ggsci::category20_d3")[6]) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5,
#             color=paletteer::paletteer_d("ggsci::category20_d3")[6]) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.6) +
#   labs(y="Interaction strength", x=NULL, title="Signals to NK cells", subtitle="Source") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_Signals_NK.pdf", height=4.5, width=8.5)
# # cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2, rel_heights=c(0.7,1))
# cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2)
# dev.off()
# 
# #####################################
# 
# 
# 
# ### Plot the signals in and out CD4T cells
# #####################################
# ###
# 
# library(ggplot2)
# library(dplyr)
# 
# CELLCHAT.SIGNAL=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.SIGNAL"]]
# 
# # take the NK cells only
# CELLCHAT.SIGNAL_NK_source=CELLCHAT.SIGNAL %>% 
#   subset((Source=="CD4T cells" | Target=="CD4T cells") & !grepl("DCs|Other cells",Target) & !grepl("DCs|Other cells",Source))
# 
# # mark the top signals
# marking_target=CELLCHAT.SIGNAL_NK_source %>%
#   subset(Prob!=0) %>%
#   group_by(Source, Target, Signal) %>%
#   summarize_at("Prob", mean) %>%
#   slice_max(Prob, n=3) %>%
#   mutate(marker=Signal) %>%
#   dplyr::select(-Prob)
# 
# plot_per_signal_NK_target=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Source=="CD4T cells"), aes(x=Agecut, y=Prob, color=Target)) +
#   facet_wrap(~Target, nrow=1) +
#   geom_point(size=0.5, alpha=0.25) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.5) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(1,2,3,5,6,8)]) +
#   labs(y="Interaction strength", x=NULL, title="Signals from CD4T cells", subtitle="Target") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# plot_per_signal_NK_source=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Target=="CD4T cells"), aes(x=Agecut, y=Prob, color=Source)) +
#   facet_wrap(~Source, nrow=1) +
#   geom_point(size=0.5, alpha=0.25, 
#              color=paletteer::paletteer_d("ggsci::category20_d3")[2]) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5,
#             color=paletteer::paletteer_d("ggsci::category20_d3")[2]) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.5) +
#   scale_color_manual(values=) +
#   labs(y="Interaction strength", x=NULL, title="Signals to CD4T cells", subtitle="Source") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_Signals_CD4T.pdf", height=4.5, width=8.5)
# # cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2, rel_heights=c(0.7,1))
# cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2)
# dev.off()
# 
# #####################################
# 
# 
# 
# ### Plot the signals in and out B cells
# #####################################
# ###
# 
# library(ggplot2)
# library(dplyr)
# 
# CELLCHAT.SIGNAL=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.SIGNAL"]]
# 
# # take the NK cells only
# CELLCHAT.SIGNAL_NK_source=CELLCHAT.SIGNAL %>% 
#   subset((Source=="B cells" | Target=="B cells") & !grepl("DCs|Other cells",Target) & !grepl("DCs|Other cells",Source))
# 
# # mark the top signals
# marking_target=CELLCHAT.SIGNAL_NK_source %>%
#   subset(Prob!=0) %>%
#   group_by(Source, Target, Signal) %>%
#   summarize_at("Prob", mean) %>%
#   slice_max(Prob, n=3) %>%
#   mutate(marker=Signal) %>%
#   dplyr::select(-Prob)
# 
# plot_per_signal_NK_target=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Source=="B cells"), aes(x=Agecut, y=Prob, color=Target)) +
#   facet_wrap(~Target, nrow=1) +
#   geom_point(size=0.5, alpha=0.25) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.5) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(1,2,3,5,6,8)]) +
#   labs(y="Interaction strength", x=NULL, title="Signals from B cells", subtitle="Target") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# plot_per_signal_NK_source=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Target=="B cells"), aes(x=Agecut, y=Prob, color=Source)) +
#   facet_wrap(~Source, nrow=1) +
#   geom_point(size=0.5, alpha=0.25, 
#              color=paletteer::paletteer_d("ggsci::category20_d3")[1]) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5,
#             color=paletteer::paletteer_d("ggsci::category20_d3")[1]) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.5) +
#   scale_color_manual(values=) +
#   labs(y="Interaction strength", x=NULL, title="Signals to B cells", subtitle="Source") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_Signals_B.pdf", height=4.5, width=8.5)
# # cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2, rel_heights=c(0.7,1))
# cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2)
# dev.off()
# 
# #####################################
# 
# 
# 
# ### Plot the signals in and out Monocytes
# #####################################
# ###
# 
# library(ggplot2)
# library(dplyr)
# 
# CELLCHAT.SIGNAL=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.SIGNAL"]]
# 
# # take the NK cells only
# CELLCHAT.SIGNAL_NK_source=CELLCHAT.SIGNAL %>% 
#   subset((Source=="Monocytes" | Target=="Monocytes") & !grepl("DCs|Other cells",Target) & !grepl("DCs|Other cells",Source))
# 
# # mark the top signals
# marking_target=CELLCHAT.SIGNAL_NK_source %>%
#   subset(Prob!=0) %>%
#   group_by(Source, Target, Signal) %>%
#   summarize_at("Prob", mean) %>%
#   slice_max(Prob, n=3) %>%
#   mutate(marker=Signal) %>%
#   dplyr::select(-Prob)
# 
# plot_per_signal_NK_target=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Source=="Monocytes"), aes(x=Agecut, y=Prob, color=Target)) +
#   facet_wrap(~Target, nrow=1) +
#   geom_point(size=0.5, alpha=0.25) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.5) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(1,2,3,5,6,8)]) +
#   labs(y="Interaction strength", x=NULL, title="Signals from Monocytes", subtitle="Target") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# plot_per_signal_NK_source=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Target=="Monocytes"), aes(x=Agecut, y=Prob, color=Source)) +
#   facet_wrap(~Source, nrow=1) +
#   geom_point(size=0.5, alpha=0.25, 
#              color=paletteer::paletteer_d("ggsci::category20_d3")[5]) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5,
#             color=paletteer::paletteer_d("ggsci::category20_d3")[5]) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.5) +
#   scale_color_manual(values=) +
#   labs(y="Interaction strength", x=NULL, title="Signals to Monocytes", subtitle="Source") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_Signals_Monocytes.pdf", height=4.5, width=8.5)
# # cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2, rel_heights=c(0.7,1))
# cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2)
# dev.off()
# 
# #####################################
# 
# 
# 
# ### Plot the signals in and out CD8T cells
# #####################################
# ###
# 
# library(ggplot2)
# library(dplyr)
# 
# CELLCHAT.SIGNAL=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.SIGNAL"]]
# 
# # take the NK cells only
# CELLCHAT.SIGNAL_NK_source=CELLCHAT.SIGNAL %>% 
#   subset((Source=="CD8T cells" | Target=="CD8T cells") & !grepl("DCs|Other cells",Target) & !grepl("DCs|Other cells",Source))
# 
# # mark the top signals
# marking_target=CELLCHAT.SIGNAL_NK_source %>%
#   subset(Prob!=0) %>%
#   group_by(Source, Target, Signal) %>%
#   summarize_at("Prob", mean) %>%
#   slice_max(Prob, n=3) %>%
#   mutate(marker=Signal) %>%
#   dplyr::select(-Prob)
# 
# plot_per_signal_NK_target=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Source=="CD8T cells"), aes(x=Agecut, y=Prob, color=Target)) +
#   facet_wrap(~Target, nrow=1) +
#   geom_point(size=0.5, alpha=0.25) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.5) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(1,2,3,5,6,8)]) +
#   labs(y="Interaction strength", x=NULL, title="Signals from CD8T cells", subtitle="Target") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# plot_per_signal_NK_source=
#   ggplot(CELLCHAT.SIGNAL_NK_source %>% subset(Target=="CD8T cells"), aes(x=Agecut, y=Prob, color=Source)) +
#   facet_wrap(~Source, nrow=1) +
#   geom_point(size=0.5, alpha=0.25, 
#              color=paletteer::paletteer_d("ggsci::category20_d3")[3]) +
#   geom_line(aes(group=Signal), linetype="solid", linewidth=0.25, alpha=0.5,
#             color=paletteer::paletteer_d("ggsci::category20_d3")[3]) +
#   ggrepel::geom_text_repel(data=. %>% left_join(marking_target, by=c("Source","Target","Signal")) %>%
#                              group_by(Source, Target, Signal) %>% filter(Prob==max(Prob)), 
#                            aes(label=marker), size=3.5, color="black", min.segment.length=0, box.padding=0.5) +
#   scale_color_manual(values=) +
#   labs(y="Interaction strength", x=NULL, title="Signals to CD8T cells", subtitle="Source") +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=10.5), 
#         plot.subtitle=element_text(size=10),
#         strip.text=element_text(size=10),
#         # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
#         axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.position="none",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_Signals_CD8T.pdf", height=4.5, width=8.5)
# # cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2, rel_heights=c(0.7,1))
# cowplot::plot_grid(plotlist=list(plot_per_signal_NK_source, plot_per_signal_NK_target), nrow=2)
# dev.off()
# 
# #####################################





### Plot the categories of receptor/ligand pairs across age
#####################################
###

library(ggplot2)
library(dplyr)

CELLCHAT.PAIR=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.PAIR"]]

CELLCHAT.PAIR_subset=CELLCHAT.PAIR %>%
  subset(Prob!=0) %>%
  group_by(Source, Target, Agecut) %>%
  add_count()

plot_signalN=
  ggplot(CELLCHAT.PAIR_subset, aes(x=Agecut, y=n, color=Target)) +
  facet_wrap(~Source, nrow=2) +
  geom_point(size=0.5, alpha=0.25) +
  geom_line(aes(group=Target), linetype="solid", linewidth=0.25, alpha=0.5) +
  labs(y="SignalN", x=NULL, title=NULL) +
  theme_minimal() +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        title=element_text(size=10.5), 
        plot.subtitle=element_text(size=10),
        strip.text=element_text(size=10),
        # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), ncol=1))

plot(plot_signalN)

#####################################



### Plot the receptor/ligand pairs underlying the altered signals: MHC-I, CLEC, CD99, MIF
#####################################
###

library(ggplot2)
library(dplyr)

CELLCHAT.PAIR=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_rough_allresults.rds")[["CELLCHAT.PAIR"]]
CELLCHAT.PAIR=CELLCHAT.PAIR %>% 
  subset(!grepl("DCs|Other cells",Target) & !grepl("DCs|Other cells",Source))

# B to CD8T
CELLCHAT.PAIR_subset=CELLCHAT.PAIR %>%
  subset(pathway_name=="MHC-I" & Source=="B cells" & Target=="CD8T cells") %>%
  subset(Prob!=0)
plot_pair=
  ggplot(CELLCHAT.PAIR_subset, aes(x=Agecut, y=Prob, color=Pair)) +
  geom_point(size=0.5) +
  geom_line(aes(group=Pair), linewidth=0.25) +
  ggrepel::geom_label_repel(data=. %>% subset(grepl("[ABCE]_CD8B",Pair)) %>%
                              group_by(Source, Target, Pair) %>% filter(Agecut==Agecut[2]),
                            aes(label=Pair), color="black", size=3.5, box.padding=0.3, min.segment.length=0,
                            show.legend=F) +
  labs(y="Interaction strength", x=NULL, title="Signals from B cells to CD8T cells", subtitle="MHC-I") +
  theme_minimal() +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        title=element_text(size=10.5), 
        plot.subtitle=element_text(size=10),
        strip.text=element_text(size=10),
        # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_color_manual(values=paletteer::paletteer_c("ggthemes::Red-Blue Diverging", length(unique(CELLCHAT.PAIR_subset$Pair))))

# all celltypes (B, CD4T, CD8T, Monocytes, NK cells, OtherT) to NK, MHCI
CELLCHAT.PAIR_subset=CELLCHAT.PAIR %>%
  subset(pathway_name=="MHC-I" & Target=="NK cells") %>%
  subset(Prob!=0)
# plot_pair=
  ggplot(CELLCHAT.PAIR_subset, aes(x=Agecut, y=Prob, color=Pair)) +
  geom_point(data=.  %>%
               subset(grepl("NKG2A|KLRC1|KLRK1",Pair)) %>%
               group_by(Target, Pair, Agecut) %>%
               summarize_at("Prob", sum),
             size=0.25) +
  geom_line(data=.  %>%
              group_by(Target, Pair, Agecut) %>%
              summarize_at("Prob", sum),
            aes(group=Pair), linewidth=0.25) +
  ggrepel::geom_label_repel(data=.  %>%
                              group_by(Target, Pair, Agecut) %>%
                              summarize_at("Prob", sum) %>% 
                              subset(grepl("NKG2A|KLRC1|KLRK1",Pair)) %>%
                              group_by(Target, Pair) %>% filter(Agecut=="51~55"),
                            aes(label=Pair), color="black", size=3.5, box.padding=0.3, min.segment.length=0,
                            show.legend=F) +
  labs(y="Interaction strength", x=NULL, title="Signals to NK cells", subtitle="MHC-I") +
  theme_minimal() +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        title=element_text(size=10.5), 
        plot.subtitle=element_text(size=10),
        strip.text=element_text(size=10),
        # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_color_manual(values=paletteer::paletteer_c("ggthemes::Red-Blue Diverging", length(unique(CELLCHAT.PAIR_subset$Pair))))

# all celltypes (B, CD4T, CD8T, Monocytes, NK cells, OtherT) to NK, CLEC
CELLCHAT.PAIR_subset=CELLCHAT.PAIR %>%
  subset(pathway_name=="CLEC" & Target=="NK cells") %>%
  subset(Prob!=0)
# plot_pair=
  ggplot(CELLCHAT.PAIR_subset, aes(x=Agecut, y=Prob, color=Pair)) +
  geom_point(data=.  %>%
               group_by(Target, Pair, Agecut) %>%
               summarize_at("Prob", sum),
             size=0.25) +
  geom_line(data=.  %>%
              group_by(Target, Pair, Agecut) %>%
              summarize_at("Prob", sum),
            aes(group=Pair), linewidth=0.25) +
  # geom_line(data=.  %>%
  #             mutate(Pair_Source=paste0(Pair,";",Source)),
  #           aes(group=Pair_Source), linewidth=0.2, alpha=0.5, linetype="dashed") +
  ggrepel::geom_label_repel(data=.  %>%
                              group_by(Target, Pair, Agecut) %>%
                              summarize_at("Prob", sum) %>% 
                              group_by(Target, Pair) %>% filter(Agecut=="51~55"),
                            aes(label=Pair), color="black", size=3.5, box.padding=0.3, min.segment.length=0,
                            show.legend=F) +
  labs(y="Interaction strength", x=NULL, title="Signals to NK cells", subtitle="CLEC") +
  theme_minimal() +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        title=element_text(size=10.5), 
        plot.subtitle=element_text(size=10),
        strip.text=element_text(size=10),
        # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_color_manual(values=paletteer::paletteer_c("ggthemes::Red-Blue Diverging", length(unique(CELLCHAT.PAIR_subset$Pair))))

# between CD8T/NK and NK, CD99
CELLCHAT.PAIR_subset=CELLCHAT.PAIR %>%
  subset(pathway_name=="CD99" & Source %in% c("CD8T cells","NK cells") & Target=="NK cells") %>%
  subset(Prob!=0)
plot_pair=
  ggplot(CELLCHAT.PAIR_subset, aes(x=Agecut, y=Prob, color=Pair)) +
  geom_point(data=. %>%
               group_by(Target, Pair, Agecut) %>%
               summarize_at("Prob", sum),
             size=0.25) +
  geom_line(data=.  %>%
              group_by(Target, Pair, Agecut) %>%
              summarize_at("Prob", sum),
            aes(group=Pair), linewidth=0.25) +
  ggrepel::geom_label_repel(data=.  %>%
                              group_by(Target, Pair, Agecut) %>%
                              summarize_at("Prob", sum) %>% 
                              group_by(Target, Pair) %>% filter(Agecut=="51~55"),
                            aes(label=Pair), color="black", size=3.5, box.padding=0.3, min.segment.length=0,
                            show.legend=F) +
  labs(y="Interaction strength", x=NULL, title="Signals between CD8T/NK and NK cells", subtitle="CD99") +
  theme_minimal() +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        plot.title.position="plot",
        title=element_text(size=10.5), 
        plot.subtitle=element_text(size=10),
        strip.text=element_text(size=10),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_color_manual(values=paletteer::paletteer_c("ggthemes::Red-Blue Diverging", length(unique(CELLCHAT.PAIR_subset$Pair))))

# NK to CD8T, MHC-I
CELLCHAT.PAIR_subset=CELLCHAT.PAIR %>%
  subset(pathway_name=="MHC-I" & Source=="NK cells" & Target=="CD8T cells") %>%
  subset(Prob!=0)
# plot_pair=
ggplot(CELLCHAT.PAIR_subset, aes(x=Agecut, y=Prob, color=Pair)) +
  geom_point(size=0.5) +
  facet_wrap(~Pair) +
  geom_line(aes(group=Pair), linewidth=0.25) +
  ggrepel::geom_label_repel(data=. %>% subset(grepl("[ABCE]_CD8B",Pair)) %>%
                              group_by(Source, Target, Pair) %>% filter(Agecut==Agecut[2]),
                            aes(label=Pair), color="black", size=3.5, box.padding=0.3, min.segment.length=0,
                            show.legend=F) +
  labs(y="Interaction strength", x=NULL, title="Signals from B cells to CD8T cells", subtitle="MHC-I") +
  theme_minimal() +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        title=element_text(size=10.5), 
        plot.subtitle=element_text(size=10),
        strip.text=element_text(size=10),
        # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_color_manual(values=paletteer::paletteer_c("ggthemes::Red-Blue Diverging", length(unique(CELLCHAT.PAIR_subset$Pair))))

# CD4T/CD8T to CD8T, MIF
CELLCHAT.PAIR_subset=CELLCHAT.PAIR %>%
  subset(pathway_name=="MIF" & Source=="CD4T cells" & Target=="CD8T cells") %>%
  subset(Prob!=0)
# plot_pair=
  ggplot(CELLCHAT.PAIR_subset, aes(x=Agecut, y=Prob, color=Pair)) +
  geom_point(size=0.5) +
  facet_wrap(~Source) +
  geom_line(aes(group=Pair), linewidth=0.25) +
  ggrepel::geom_label_repel(data=. %>% 
                              # subset(grepl("[ABCE]_CD8A",Pair)) %>%
                              group_by(Source, Target, Pair) %>% filter(Agecut=="51~55"),
                            aes(label=Pair), color="black", size=3.5, box.padding=0.3, min.segment.length=0,
                            show.legend=F) +
  labs(y="Interaction strength", x=NULL, title="Signals from CD4T to CD8T cells", subtitle="MHC-I") +
  theme_minimal() +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        title=element_text(size=10.5), 
        plot.subtitle=element_text(size=10),
        strip.text=element_text(size=10),
        # axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_color_manual(values=paletteer::paletteer_c("ggthemes::Red-Blue Diverging", length(unique(CELLCHAT.PAIR_subset$Pair))))

  
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ### Signals
# cellchat.signal=list()
# for (i in 1:length(cellchat)) {
#   tempt=cellchat[[i]]@netP$prob
#   cellchat.signal_j=list()
#   for (j in 1:dim(tempt)[3]) {
#     cellchat.signal_j[[j]]=as.data.frame(tempt[,,j]) %>%
#       tibble::rownames_to_column(var="Source") %>%
#       tidyr::pivot_longer(!Source, names_to="Target", values_to="Prob")
#     cellchat.signal_j[[j]]$Signal=dimnames(tempt)[[3]][j]
#   }
#   cellchat.signal[[i]]=data.table::rbindlist(cellchat.signal_j)
#   cellchat.signal[[i]]$Agecut=agecuts[i]
# }
# CELLCHAT.SIGNAL=data.table::rbindlist(cellchat.signal)
# 
# # take only the most contributing signals
# chosen_signal=c("CD99","CLEC","LCK","MHC-I","MIF")
# 
# plot_per_signal=list()
# for (i in 1:length(chosen_signal)) {
#   tempt_df=CELLCHAT.SIGNAL %>% subset(Source %in% c("B cells","DCs","Monocytes","CD4T cells","CD8T cells") &
#                                         Target %in% c("B cells","DCs","Monocytes","CD4T cells","CD8T cells") &
#                                         Signal==chosen_signal[i])
#   plot_per_signal[[i]]=
#     ggplot(tempt_df, aes(x=Agecut, y=Prob, color=Target)) +
#     facet_grid(cols=vars(Source)) +
#     geom_point(size=0.5) +
#     geom_line(aes(group=Target), linetype="dashed", linewidth=0.5) +
#     scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
#     labs(y=unique(CELLCHAT.SIGNAL$Signal)[i], x=NULL) +
#     theme_minimal() +
#     theme(title=element_text(size=9),
#           strip.text=element_blank(),
#           axis.text.x=element_blank(),
#           axis.text.y=element_text(size=9),
#           axis.title.y=element_text(size=10),
#           legend.position="none",
#           legend.title=element_text(size=10),
#           legend.text=element_text(size=9)) +
#     guides(color=guide_legend(override.aes=list(size=2)))
# }
# plot_per_signal[1]=plot_per_signal[1] %>%
#   lapply(., function(p) p + theme(strip.text=element_text(size=10)) + labs(title="Source"))
# 
# plot_per_signal[length(plot_per_signal)]=plot_per_signal[length(plot_per_signal)] %>%
#   lapply(., function(p) p + theme(axis.text.x=element_text(size=7.5, angle=90, hjust=0, vjust=0.5),
#                                   legend.position="bottom",
#                                   legend.direction="horizontal") +
#                             scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#                             guides(color=guide_legend(nrow=2, byrow=TRUE)))
# 
# Signal_plot_all=cowplot::plot_grid(plotlist=plot_per_signal, align="hv", ncol=1, rel_heights=c(1,1,1,1,1.8))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cellchat_Signals.pdf", height=8, width=5)
# plot(Signal_plot_all)
# dev.off()
# 
# # plot the top signals in source-CD4/CD8 ~ target-CD8
# CELLCHAT.SIGNAL_topSignals=CELLCHAT.SIGNAL %>%
#   subset(Source %in% c("CD4T cells","CD8T cells") & Target=="CD8T cells") %>%
#   group_by(Source, Target, Signal) %>%
#   summarize_at("Prob", mean) %>%
#   slice_max(Prob, n=5) %>%
#   ungroup() %>%
#   dplyr::select(Signal) %>%
#   tibble::deframe() %>% .[!duplicated(.)]
# 
# # weight_plot_all=
#   ggplot(CELLCHAT.SIGNAL %>% subset(Source %in% c("CD4T cells","CD8T cells") & Target=="CD8T cells" & Signal %in% CELLCHAT.SIGNAL_topSignals), 
#          aes(x=Agecut, y=Prob, color=Signal)) +
#   facet_grid(cols=vars(Source)) +
#   geom_point(size=0.5) +
#   geom_line(aes(group=Signal), linetype="dashed", linewidth=0.5) +
#   scale_color_manual(values=paletteer::paletteer_d("ggthemes::few_Medium")) +
#   labs(title="Signals to CD8T cells", y="Interaction strength", x=NULL) +
#   theme_minimal() +
#   theme(plot.background=element_rect(fill="white", color="transparent"),
#         title=element_text(size=9),
#         strip.text=element_text(size=10),
#         axis.text.x=element_text(size=7.5, angle=90, hjust=0, vjust=0.5),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2)) +
#   guides(color=guide_legend(override.aes=list(size=2)))
#   
# 
# #####################################



### Plot the expression of the key ligands/receptors
#####################################
###

### Load correlation df
cor_data=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")

key_ligands=list(`B cells`=c("HLA-A","HLA-B","HLA-C","HLA-E"),
                 `NK cells`=c("KLRK1","CD94","NKG2A","KLRC1","CD99","HLA-A","HLA-B","HLA-C","HLA-E"),
                 `CD8T cells`=c("CD99","CD8A","CD8B","CD74","CXCR4"),
                 `CD4T cells`=c("MIF"),
                 `Monocytes`=c("HLA-E"),
                 `OtherT`=c("HLA-E")
                 )
cor_data_df=cor_data %>% 
  subset(gene %in% unlist(key_ligands) & analysis=="all" & celltypes!="DCs" & celltypes!="Other cells") %>%
  dplyr::select(celltypes, gene, rho) %>%
  tidyr::pivot_wider(names_from="gene", values_from="rho") %>%
  tibble::column_to_rownames("celltypes")
cor_data_df_cor=cor_data %>% 
  subset(gene %in% unlist(key_ligands) & analysis=="all" & celltypes!="DCs" & celltypes!="Other cells") %>%
  dplyr::select(celltypes, gene, rho_pval) %>%
  mutate(rho_pval=ifelse(rho_pval<0.0001,"****",
                         ifelse(rho_pval>=0.0001 & rho_pval<0.001, "***",
                                ifelse(rho_pval>=0.001 & rho_pval<0.01, "**",
                                       ifelse(rho_pval>=0.01 & rho_pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="gene", values_from="rho_pval") %>%
  tibble::column_to_rownames("celltypes") %>%
  as.matrix()
  
### Plot
library(ComplexHeatmap)
range(cor_data_df) # check
col_fun_sub=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))

ht_list=
  Heatmap(cor_data_df,
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(cor_data_df_cor, i, j)
            grid.text(v, x, y, rot=90, just="top",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          name="mat", show_heatmap_legend=F,
          heatmap_legend_param=list(
            title=expression(Spearman~rho),
            legend_height=unit(4.6, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
          ),
          rect_gp=gpar(col="white", lwd=0.5),
          col=col_fun_sub,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="Terekhova et al. (2023) dataset", 
          column_title_gp=gpar(fontsize=11),
          width=ncol(cor_data_df)*unit(4.5, "mm"),
          height=nrow(cor_data_df)*unit(8, "mm"))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_rough_KeyLRpairs_expr.pdf", height=3, width=3)
whole_ht
dev.off()

#####################################
  
  
  
  
  

##########################################################################
### Analyze on the celltypes at the inter level
##########################################################################


### Prepare the cellchat objs for the analysis
#####################################
###

object=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
object=Seurat::NormalizeData(object)

### Prepare Agecut
object[[]]=object[[]] %>%
  dplyr::select(-agecut) %>%
  dplyr::rename(age=Age) %>%
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:44), "41~44", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

### Divide the obj by agecut and take only assigned cells
AGECUT_name=names(table(object$agecut))
object_list=list()
for (i in 1:length(AGECUT_name)) {
  object_list[[i]]=subset(object, agecut==AGECUT_name[i])
}
names(object_list)=AGECUT_name

### Create the cellchat objects
cellchat=list()
for (i in 1:length(object_list)) {
  DefaultAssay(object_list[[i]])="RNA"
  cellchat[[i]]=createCellChat(object=object_list[[i]], meta=object_list[[i]]@meta.data, group.by="Annot.inter")
}

### load the cellchat database
CellChatDB=CellChatDB.human
for (i in 1:length(cellchat)) {
  cellchat[[i]]@DB=CellChatDB
}
saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter.rds")

### Identify overexpression
processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
for (i in 1:length(cellchat)) {
  cellchat[[i]]=subsetData(cellchat[[i]])
  cellchat[[i]]=identifyOverExpressedGenes(cellchat[[i]])
  cellchat[[i]]=identifyOverExpressedInteractions(cellchat[[i]])
  setTxtProgressBar(processbar, i)
}
close(processbar)

### Evaluate the communications
options(future.globals.maxSize=1*1024^3) # first try computeCommunProb, if error then increase the memory: 850MiB * 1024^2 = 891289600
processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
for (i in 1:length(cellchat)) {
  cellchat[[i]]=computeCommunProb(cellchat[[i]], raw.use=TRUE, population.size=TRUE)
  cellchat[[i]]=filterCommunication(cellchat[[i]], min.cells=10)
  setTxtProgressBar(processbar, i)
}
close(processbar)
saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter.rds")

### Calculate the communication networks
processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
for (i in 1:length(cellchat)) {
  cellchat[[i]]=computeCommunProbPathway(cellchat[[i]])
  cellchat[[i]]=aggregateNet(cellchat[[i]])
  setTxtProgressBar(processbar, i)
}
close(processbar)
saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter.rds")

### Calculate the communication networks
for (i in 1:length(cellchat)) {
  cellchat[[i]]=netAnalysis_computeCentrality(cellchat[[i]], slot.name="netP")
}
names(cellchat)=AGECUT_name
saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter.rds")

#####################################



### Merge all the cellchat objs and plot the interaction comparison
#####################################
###

library(CellChat)
### Merge the cellchat objs all in one
cellchat=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter.rds")
cellchat.list=mergeCellChat(cellchat, add.names=names(cellchat), cell.prefix=TRUE)
saveRDS(cellchat.list, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter_mergeCellchatList.rds")

cellchat.list=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter_mergeCellchatList.rds")
plot_count=compareInteractions(cellchat.list, show.legend=F, measure="count")
plot_weight=compareInteractions(cellchat.list, show.legend=F, measure="weight")
plot_interactions=plot_count+plot_weight

### Plot by ggplot2
getAnywhere("compareInteractions")
count_df=as.data.frame(sapply(cellchat.list@net, function(x) sum(x$count)))
weight_df=as.data.frame(sapply(cellchat.list@net, function(x) sum(x$weight)))
df=merge(count_df, weight_df, by=0); colnames(df)=c("agecut","count","weight"); df$agecut=as.character(df$agecut)
df=df %>% tidyr::pivot_longer(cols=c("count","weight"), names_to="measure", values_to="value") %>%
  mutate(measure=ifelse(measure=="count", "Number of interactions","Interaction strength"))
df$measure=forcats::fct_relevel(df$measure, c("Number of interactions","Interaction strength"))

plot_interactions=
  ggplot(df, aes(x=agecut, y=value)) +
  geom_bar(stat="identity", fill="grey60", width=0.8) +
  facet_wrap(~measure, scales="free_y", labeller=labeller(scales::label_wrap(10))) +
  theme_minimal() +
  labs(x=NULL, y=NULL) +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=9),
        title=element_text(size=10),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        strip.text=element_text(size=11),
        legend.position="left",
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cellchat_inter_compareall_compareInteractions.pdf", width=5.5, height=4)
plot(plot_interactions)
dev.off()

#####################################




# ### Compare the signals in the Cellchat obj at inter level between different agecuts, regardless of celltypes
# #####################################
# ###
# 
# library(CellChat)
# 
# ### Load the obj
# cellchat.list=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_inter_mergeCellchatList.rds")
# 
# ### Fetch the contribution data using CellChat plot
# plot_ranknet=
#   rankNet(cellchat.list, mode="comparison", comparison=c(1:length(unique(cellchat.list@meta$agecut))), stacked=T, do.stat=TRUE)
# df=plot_ranknet$data %>% select(name, contribution, group)
# 
# # order based on the contribution and select only the top15 to show
# signal_order=df %>%
#   group_by(name) %>%
#   summarize_at("contribution", sum) %>%
#   arrange(desc(contribution)) %>%
#   ungroup() %>%
#   dplyr::select(name) %>% tibble::deframe() %>% as.character() %>% .[1:15]
# df=df %>% subset(name %in% signal_order)
# df$name=forcats::fct_relevel(df$name, signal_order)
# 
# # plot by ggplot2
# information_flow_plot=
#   ggplot(df, aes(x=name, y=contribution, color=group)) +
#   geom_point(size=1) +
#   geom_hline(yintercept=0, linewidth=0.25, linetype="dotted") +
#   scale_color_manual(values=paletteer::paletteer_c("ggthemes::Red-Blue Diverging", length(unique(df$group)))) +
#   theme_classic() +
#   labs(x=NULL, y="Information flow") +
#   theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
#         axis.text.y=element_text(size=9),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=11),
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   ggbreak::scale_y_break(c(0.06, 0.11), scales=0.2) +
#   guides(color=guide_legend(title="age", override.aes=list(size=2), ncol=3))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cellchat_inter_GeneralInformationFlow.pdf", width=6, height=3, onefile=FALSE)
# information_flow_plot
# dev.off()
# 
# #####################################



### Extract the results for inter
#####################################
###

library(CellChat)
library(dplyr)

cellchat=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter.rds")
agecuts=names(cellchat)

### Count and Weight
cellchat.count=cellchat.weight=list()
for (i in 1:length(cellchat)) {
  cellchat.count[[i]]=cellchat[[i]]@net$count
  cellchat.count[[i]]=as.data.frame(cellchat.count[[i]]) %>%
    tibble::rownames_to_column(var="Source") %>%
    tidyr::pivot_longer(!Source, names_to="Target", values_to="Count")
  cellchat.count[[i]]$Agecut=agecuts[i]
  
  cellchat.weight[[i]]=cellchat[[i]]@net$weight
  cellchat.weight[[i]]=as.data.frame(cellchat.weight[[i]]) %>%
    tibble::rownames_to_column(var="Source") %>%
    tidyr::pivot_longer(!Source, names_to="Target", values_to="Weight")
  cellchat.weight[[i]]$Agecut=agecuts[i]
}
CELLCHAT.COUNT=data.table::rbindlist(cellchat.count)
CELLCHAT.WEIGHT=data.table::rbindlist(cellchat.weight)

### Signals
cellchat.signal=list()
for (i in 1:length(cellchat)) {
  tempt=cellchat[[i]]@netP$prob
  cellchat.signal_j=list()
  for (j in 1:dim(tempt)[3]) {
    cellchat.signal_j[[j]]=as.data.frame(tempt[,,j]) %>%
      tibble::rownames_to_column(var="Source") %>%
      tidyr::pivot_longer(!Source, names_to="Target", values_to="Prob")
    cellchat.signal_j[[j]]$Signal=dimnames(tempt)[[3]][j]
  }
  cellchat.signal[[i]]=data.table::rbindlist(cellchat.signal_j)
  cellchat.signal[[i]]$Agecut=agecuts[i]
}
CELLCHAT.SIGNAL=data.table::rbindlist(cellchat.signal)

### Pairs
cellchat.pair=list()
for (i in 1:length(cellchat)) {
  tempt=cellchat[[i]]@net$prob
  pair_list_names=names(tempt[1,1,])
  pair_list=lapply(1:dim(tempt)[3], function(idx) sum(tempt[,,idx])) %>% unlist(); names(pair_list)=pair_list_names
  pair_list=pair_list[pair_list!=0]
  match_pair_signal=cellchat[[i]]@LR$LRsig %>% .[pair_list_names,] %>% dplyr::rename(Pair=interaction_name)
  
  cellchat.pair_j=list()
  for (j in 1:length(pair_list)) {
    cellchat.pair_j[[j]]=as.data.frame(tempt[,,names(pair_list)[j]]) %>%
      tibble::rownames_to_column(var="Source") %>%
      tidyr::pivot_longer(!Source, names_to="Target", values_to="Prob")
    cellchat.pair_j[[j]]$Pair=names(pair_list)[j]
  }
  cellchat.pair[[i]]=data.table::rbindlist(cellchat.pair_j)
  cellchat.pair[[i]]$Agecut=agecuts[i]
  cellchat.pair[[i]]=cellchat.pair[[i]] %>% left_join(match_pair_signal, by=c("Pair"))
}
CELLCHAT.PAIR=data.table::rbindlist(cellchat.pair)

saveRDS(list(CELLCHAT.COUNT=CELLCHAT.COUNT, CELLCHAT.WEIGHT=CELLCHAT.WEIGHT, CELLCHAT.SIGNAL=CELLCHAT.SIGNAL, CELLCHAT.PAIR=CELLCHAT.PAIR), 
        "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter_allresults.rds")

#####################################



### Plot the interaction counts/weights among inter celltypes across age
#####################################
###

###
library(ggplot2)
library(dplyr)

CELLCHAT.COUNT=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter_allresults.rds")[["CELLCHAT.COUNT"]]
CELLCHAT.WEIGHT=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter_allresults.rds")[["CELLCHAT.WEIGHT"]]
CELLCHAT.PAIR=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_inter_allresults.rds")[["CELLCHAT.PAIR"]]

# B, CD4T, NK to CD8T, all pairs
cellchat_subset=
  CELLCHAT.PAIR %>% 
  subset(grepl("B\\.|NK\\.|CD4T\\.",Source) & grepl("CD8T\\.",Target) & 
           Pair %in% c("HLA-A_CD8B","HLA-B_CD8B","HLA-C_CD8B","HLA-E_CD8B","CD99_CD99","HLA-A_CD8A","HLA-B_CD8A","HLA-C_CD8A","HLA-E_CD8A","MIF_CD74_CXCR4","MIF_CD74_CD44")) %>%
  group_by(Source, Target, Agecut) %>%
  summarize_at("Prob", sum)
# plot weight
# weight_plot_cut=
  ggplot(cellchat_subset, 
         aes(x=Agecut, y=Prob, color=Target)) +
  facet_grid(cols=vars(Source)) +
  geom_point(size=0.5, alpha=0.25) +
  geom_line(aes(group=Target), linetype="solid", linewidth=0.25) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title="Source", y="Interaction strength", x=NULL) +
  theme_minimal() +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        title=element_text(size=9),
        strip.text=element_text(size=10),
        # axis.text.x=element_text(size=7.5, angle=90, hjust=0, vjust=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_x_discrete(guide=guide_axis(n.dodge=2)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1)))

# all to NK cells, all pairs
cellchat_subset=
  CELLCHAT.PAIR %>% 
  subset(grepl("B\\.|CD4T\\.|CD8T\\.[CnT][^r]|Mono\\.|NK\\.|OtherT\\.[^d]",Source) & grepl("NK\\.",Target) & 
           Pair %in% c("HLA-E_KLRK1","HLA-E_CD94:NKG2A","HLA-E_KLRC1","CD99_CD99")) %>%
  group_by(Source, Target, Agecut) %>%
  summarize_at("Prob", sum)
# plot weight
# weight_plot_cut=
  ggplot(cellchat_subset, 
         aes(x=Agecut, y=Prob, color=Target)) +
  facet_grid(cols=vars(Source)) +
  geom_point(size=0.5, alpha=0.25) +
  geom_line(aes(group=Target), linetype="solid", linewidth=0.25) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title="Source", y="Interaction strength", x=NULL) +
  theme_minimal() +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        title=element_text(size=9),
        strip.text=element_text(size=10),
        # axis.text.x=element_text(size=7.5, angle=90, hjust=0, vjust=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_x_discrete(guide=guide_axis(n.dodge=2)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1)))

#####################################





##########################################################################
### Analyze on the celltypes at the detailed level
##########################################################################


### Prepare the cellchat objs for the analysis
#####################################
###

object=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
object=subset(object, Annot.rough %in% c("B cells","CD4T cells","CD8T cells","NK cells","Monocytes","OtherT"))
object=Seurat::NormalizeData(object)

### Prepare Agecut
object[[]]=object[[]] %>%
  dplyr::select(-agecut) %>%
  dplyr::rename(age=Age) %>%
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:44), "41~44", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

### Divide the obj by agecut and take only assigned cells
AGECUT_name=names(table(object$agecut))
object_list=list()
for (i in 1:length(AGECUT_name)) {
  object_list[[i]]=subset(object, agecut==AGECUT_name[i])
}
names(object_list)=AGECUT_name

### Create the cellchat objects
cellchat=list()
for (i in 1:length(object_list)) {
  DefaultAssay(object_list[[i]])="RNA"
  cellchat[[i]]=createCellChat(object=object_list[[i]], meta=object_list[[i]]@meta.data, group.by="Annot.detailed")
}

### load the cellchat database
CellChatDB=CellChatDB.human
for (i in 1:length(cellchat)) {
  cellchat[[i]]@DB=CellChatDB
}
saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_detailed.rds")

### Identify overexpression
processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
for (i in 1:length(cellchat)) {
  cellchat[[i]]=subsetData(cellchat[[i]])
  cellchat[[i]]=identifyOverExpressedGenes(cellchat[[i]])
  cellchat[[i]]=identifyOverExpressedInteractions(cellchat[[i]])
  setTxtProgressBar(processbar, i)
}
close(processbar)

### Evaluate the communications
options(future.globals.maxSize=1*1024^3) # first try computeCommunProb, if error then increase the memory: 850MiB * 1024^2 = 891289600
processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
for (i in 1:length(cellchat)) {
  cellchat[[i]]=computeCommunProb(cellchat[[i]], raw.use=TRUE, population.size=TRUE)
  cellchat[[i]]=filterCommunication(cellchat[[i]], min.cells=10)
  setTxtProgressBar(processbar, i)
}
close(processbar)
saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_detailed.rds")

### Calculate the communication networks
processbar=txtProgressBar(min=0, max=length(cellchat), style=3, width=60, char="=")
for (i in 1:length(cellchat)) {
  cellchat[[i]]=computeCommunProbPathway(cellchat[[i]])
  cellchat[[i]]=aggregateNet(cellchat[[i]])
  setTxtProgressBar(processbar, i)
}
close(processbar)
saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_detailed.rds")

### Calculate the communication networks
for (i in 1:length(cellchat)) {
  cellchat[[i]]=netAnalysis_computeCentrality(cellchat[[i]], slot.name="netP")
}
names(cellchat)=AGECUT_name
saveRDS(cellchat, "~/Project_PBMCage/Immunity/Tempt_RDS/Immunity_Cellchat_detailed.rds")

#####################################


