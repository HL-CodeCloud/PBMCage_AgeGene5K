library(condiments)
library(slingshot)
library(dplyr)
library(tidyr)
library(ggplot2)

### Prepare the obj
#####################################
###

### Load the data
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
CD8T=subset(THEObj, Annot.rough=="DCs" | Annot.inter=="Other.progenitor") # add progenitors so it won't be all terminals

### Run Umap
CD8T=CD8T %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures() %>%
  Seurat::ScaleData() %>%
  Seurat::RunPCA() %>%
  Seurat::FindNeighbors(dims=1:30) %>%
  Seurat::RunUMAP(dims=1:30)

### Check if enough observations in each group (>20)
table(CD8T$Annot.inter, CD8T$agecut)
table(CD8T$Annot.detailed, CD8T$agecut)
meta_data=CD8T[[]]
meta_data=meta_data %>% 
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Other.progenitor.HSC","Other.progenitor.CLP","Other.progenitor.MEP","Other.progenitor.MPP","Other.progenitor.NKP"), 
                               "Other.progenitor.HSPC", Annot.detailed))
table(meta_data$Annot.detailed, meta_data$agecut)
CD8T$Annot.detailed=meta_data$Annot.detailed

### Extract the umap results
df=CD8T@reductions$umap@cell.embeddings
df=df %>% as.data.frame() %>% mutate(cellid=rownames(.)) %>% right_join(CD8T@meta.data %>% mutate(cellid=rownames(.)), by="cellid")
# check if enough observations in each group
table(df$Annot.inter, df$agecut)
table(df$Annot.detailed, df$agecut)
# if not enough observations in some groups (should more than >50), then merge
df=df %>% 
  mutate(agecut_merged=ifelse(age %in% c(19:68), "19~68", "69~97"))
table(df$Annot.inter, df$agecut_merged)
table(df$Annot.detailed, df$agecut_merged)
saveRDS(df, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")

###  Plot the cells across conditions on the umap
conditions=names(table(df$agecut_merged))
ggplot(df, aes(x=umap_1, y=umap_2, color=match(agecut_merged, conditions))) +
  ggrastr::geom_point_rast() +
  scale_color_gradient(low="white",high="brown3")

#####################################



### Analyze topology difference
#####################################
###

df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")

#### Visualize imbalance score
scores=imbalance_score(Object=df %>% select(umap_1, umap_2) %>% as.matrix(),
                       conditions=df$agecut_merged)
df$scores=scores$scores; df$scaled_scores=scores$scaled_scores
saveRDS(df, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")

ggplot(df, aes(x=umap_1, y=umap_2, color=scaled_scores)) +
  ggrastr::geom_point_rast() +
  scale_color_gradient(low="white",high="brown3")

### Run the topology test
rd=as.matrix(df[, c("umap_1", "umap_2")])
sds=slingshot(rd, df$Annot.detailed, start.clus="Other.progenitor.HSPC",
              allow.breaks=FALSE) # Annot.detailed contains inadequate cells per group
saveRDS(sds, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_origin.rds")

top_res=topologyTest(sds=sds, conditions=df$agecut_merged)
saveRDS(top_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_topologyTest.rds")

knitr::kable(top_res)
# turns out not to share common trajectories

#####################################



### Combine multiple trajectories
#####################################
###

### Get slingshot results in each condition
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
sds=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_origin.rds")
top_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_topologyTest.rds")

if (top_res$p.value<0.05) {
  sdss=slingshot_conditions(sds, df$agecut_merged,
                            approx_points=FALSE,
                            extend="n",
                            reweight=FALSE, reassign=FALSE)
  sdss$condition_id=names(sdss)
  n_condition=length(sdss$condition_id)
  n_lineage=length(slingLineages(sdss[[1]]))
  sdss$mapping=matrix(rep(1:n_lineage, each=n_condition),
                      nrow=n_lineage, ncol=n_condition,
                      byrow=TRUE)
  # merge the slingshot results
  sds_merged=do.call(merge_sds, sdss)
  saveRDS(sds_merged, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_merged.rds")
}

#####################################



### Analyze differential progression
# ... * interpretation: difference in cell development/pseudotime in each lineage among different conditions
#####################################
###

### Differential progression analysis
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_merged.rds")
prog_res=progressionTest(sds_merged, conditions=df$agecut_merged, global=TRUE, lineages=TRUE)
knitr::kable(prog_res)
saveRDS(prog_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_progressionTest.rds")
# turns out to be with sig. general difference, but not for any lineage specifically

### Plot the progression
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_merged.rds")

slingLineages(sds_merged)

df_merged=full_join(
  df %>% select(cellid, umap_1, umap_2, Annot.detailed, agecut_merged) %>%
    mutate(cells=paste0("Cell-",rownames(.))),
  slingPseudotime(sds_merged) %>%
    as.data.frame() %>%
    mutate(cells=rownames(.)),
  by="cells"
) %>%
  tidyr::pivot_longer(starts_with("Lineage"), names_to="Curve", values_to="pst")

ggplot(df_merged, aes(x=pst)) +
  geom_density(alpha=0.4, aes(fill=agecut_merged), col="transparent") +
  geom_density(aes(col=agecut_merged), fill="transparent", size=0.5) +
  guides(col="none") +
  scale_fill_brewer(palette="Accent") +
  scale_color_brewer(palette="Accent") +
  labs(x="Pseudotime", fill="Type") +
  facet_wrap(~Curve, scales="free_x")

###
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
sds=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_origin.rds")
conditions=unique(df$agecut_merged)

midpoint=df %>% group_by(Annot.detailed, agecut_merged) %>%
  dplyr::summarise(umap_1=mean(umap_1), umap_2=mean(umap_2), .groups=NULL)

plot_pro=
  ggplot(df, aes(x=umap_1, y=umap_2)) +
  facet_wrap(~agecut_merged) +
  geom_hex(bins=45, aes(fill=after_stat(density))) +
  scale_fill_distiller(palette="Reds", direction=1) +
  geom_point(data=slingMST(sds, as.df=TRUE), inherit.aes=F, aes(x=umap_1, y=umap_2), size=1) +
  geom_path(data=slingMST(sds, as.df=TRUE), inherit.aes=F, aes(x=umap_1, y=umap_2, group=Lineage), size=0.5) +
  ggrepel::geom_text_repel(data=midpoint %>% group_by(Annot.detailed) %>% summarize_at(c("umap_1","umap_2"), mean) %>%
                             mutate(agecut_merged="19~68"),
                           inherit.aes=FALSE, force_pull=100,
                           aes(x=umap_1, y=umap_2, label=stringr::str_wrap(gsub("\\."," ",gsub("^DC\\.","",Annot.detailed)), 10)),
                           alpha=0.8,
                           color="black", lineheight=0.75, min.segment.length=0, box.padding=0.5
  ) +
  guides(fill="none") +
  theme_minimal() +
  theme(strip.text=element_text(size=10), 
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        plot.background=element_rect(color="transparent", fill="white"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Pseudotime_Trajectory_DCs_progression.pdf", height=5, width=7)
plot(plot_pro)
dev.off()

#####################################



### Analyze differential differentiation
# (Use only when common trajectory is found!!!)
# (The equivalent function in the case of different trajectories is the "fateSelectionTest" below!!!)
# ... * interpretation: difference in cell differentiation preference to a certain lineage among different conditions
#####################################
###

top_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_topologyTest.rds")

if (top_res$p.value>0.05) {
  sds=readRDS( "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_origin.rds")
  
  ### Check the number of lineages
  n_lineage=length(unique(slingMST(sds, as.df=TRUE)$Lineage))
  n_lineage
  
  ### Differential differentiation analysis
  if (n_lineage>1) {
    dif_res=differentiationTest(sds, conditions=df$agecut_merged, global=TRUE, pairwise=TRUE)
    knitr::kable(dif_res)
    saveRDS(dif_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_differentiationTest.rds")
  }
}

#####################################



### Analyze differential fate selection
# ... * interpretation: difference in cell differentiation preference to a certain lineage among different conditions
#####################################
###

###
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_merged.rds")

fate_res=fateSelectionTest(sds_merged, conditions=df$agecut_merged, global=TRUE, pairwise=TRUE)
knitr::kable(fate_res)
# turns out to be no difference
saveRDS(fate_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_fateSelectionTest.rds")

### Plot the fates
weights=condiments:::.sling_reassign(sds_merged)
df_difftest=df %>%
  mutate(cells=paste0("Cell-",rownames(.))) %>%
  full_join(
    weights %>%
      as.data.frame() %>%
      mutate(cells=rownames(.)) %>%
      dplyr::rename("Lineage1"=V1, "Lineage2"=V2, "Lineage3"=V3) %>%
      tidyr::pivot_longer(starts_with("Lineage"), names_to="Curve", values_to="weights")
  )
df_w=df_difftest %>%
  group_by(cells) %>%
  mutate(weights=weights/sum(weights)) %>%
  ungroup() %>%
  group_by(agecut_merged, Curve) %>%
  summarise(weights=mean(weights), .groups=NULL)
# p2=
ggplot(df_w, aes(x=Curve, fill=agecut_merged, y=weights)) +
  geom_col(position="dodge") +
  scale_fill_brewer(palette="Accent") +
  labs(x = "", y="Mean weight")

#####################################



### Plot the trajectories
#####################################
###

df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_merged.rds")
# load the tests results
top_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_topologyTest.rds")
prog_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_progressionTest.rds")
fate_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_fateSelectionTest.rds")

df$agecut_merged=forcats::fct_relevel(df$agecut_merged, names(table(df$agecut_merged)))
p=
  ggplot(df, aes(x=umap_1, y=umap_2, col=agecut_merged, alpha=agecut_merged)) +
  ggrastr::geom_point_rast(size=0.01, shape=20,
                           color=paletteer::paletteer_d("ggsci::category20_d3")[4]) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  scale_alpha_manual(values=seq(0.1,0.2,length.out=length(unique(df$agecut_merged)))) +
  theme_light() +
  labs(x="UMAP_1", y="UMAP_2", title="DCs",
       subtitle=paste0("Progression: p=",ifelse(prog_res$p.value[2]<0.001, sprintf("%0.01e",prog_res$p.value[2]), round(prog_res$p.value[2],3)), # general p<0.05 but no lineage-specific p<0.05, therefore use the first specific p
                       ";    Fate selection: p=",ifelse(fate_res$p.value[1]<0.001, sprintf("%0.01e",fate_res$p.value[1]), round(fate_res$p.value[1],3)))) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        plot.subtitle=element_text(size=10),
        title=element_text(size=10),
        legend.position="bottom",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  guides(alpha="none")

midpoint=df %>% group_by(Annot.detailed, agecut_merged) %>%
  dplyr::summarise(umap_1=mean(umap_1), umap_2=mean(umap_2), .groups=NULL)

edges=lapply(slingLineages(sds_merged), function(lin){
  from=lin[1:(length(lin)-1)]
  to=lin[2:length(lin)]
  return(data.frame("from"=from, "to"=to))}
) %>%
  bind_rows()
for (batch in unique(df$agecut_merged)) {
  cl_batch=midpoint %>% filter(agecut_merged==batch)
  edges_batch <- left_join(edges, cl_batch, by = c("from"="Annot.detailed")) %>%
    left_join(cl_batch %>%
                dplyr::rename("umap1_end"="umap_1", "umap2_end"="umap_2") %>%
                select(-agecut_merged),
              by=c("to"="Annot.detailed"))
  p=
    p +
    geom_segment(data=edges_batch, inherit.aes=FALSE,
                 aes(x=umap_1, y=umap_2, xend=umap1_end, yend=umap2_end, col=agecut_merged),
                 alpha=1,
                 lineend="round", linejoin="round",
                 size=0.5)
}

p_withLabel=
  p +
  ggrepel::geom_label_repel(data=midpoint %>% group_by(Annot.detailed) %>% 
                              mutate(Annot.detailed=ifelse(Annot.detailed=="Other.progenitor.HSPC","HSPC",Annot.detailed)) %>% 
                              summarize_at(c("umap_1","umap_2"), mean),
                            inherit.aes=FALSE, force_pull=100,
                            aes(x=umap_1, y=umap_2, label=stringr::str_wrap(gsub("mem\\.","mem ",Annot.detailed), 10)),
                            alpha=0.8, force=2,
                            color="black") +
  guides(color=guide_legend(title="age", title.position="left", title.hjust=0.5, nrow=1))

p_withLabel

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Pseudotime_Trajectory_DCs.pdf", height=5, width=7)
plot(p_withLabel)
dev.off()

#####################################



# ### Analyze differential progression
# # ... * interpretation: difference in cell development/pseudotime in each lineage among different conditions
# #####################################
# ###
# 
# ### Differential progression analysis
# df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
# sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_origin.rds")
# prog_res=progressionTest(sds_merged, conditions=df$agecut_merged, global=TRUE, lineages=TRUE)
# knitr::kable(prog_res)
# saveRDS(prog_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_progressionTest.rds")
# 
# ### Plot the progression
# df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
# sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_origin.rds")
# df_merged=full_join(
#   df %>% select(cellid, umap_1, umap_2, Annot.detailed, agecut_merged) %>%
#     mutate(cells=paste0("Cell-",rownames(.))),
#   slingPseudotime(sds_merged) %>%
#     as.data.frame() %>%
#     mutate(cells=rownames(.)),
#   by="cells"
# ) %>%
#   tidyr::pivot_longer(starts_with("Lineage"), names_to="Curve", values_to="pst")
# 
# ggplot(df_merged, aes(x=pst)) +
#   geom_density(alpha=0.4, aes(fill=agecut_merged), col="transparent") +
#   geom_density(aes(col=agecut_merged), fill="transparent", size=1.5) +
#   guides(col="none") +
#   scale_fill_brewer(palette="Accent") +
#   scale_color_brewer(palette="Accent") +
#   labs(x="Pseudotime", fill="Type") +
#   facet_wrap(~Curve, scales="free_x")
# 
# #####################################
# 
# 
# 
# ### Analyze differential differentiation
# # (Use only when common trajectory is found!!!)
# # (The equivalent function in the case of different trajectories is the "fateSelectionTest" below!!!)
# # ... * interpretation: difference in cell differentiation preference to a certain lineage among different conditions
# #####################################
# ###
# 
# top_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_topologyTest.rds")
# 
# if (top_res$p.value>0.05) {
#   sds=readRDS( "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs.slim_sds_origin.rds")
#   
#   ### Check the number of lineages
#   n_lineage=length(unique(slingMST(sds, as.df=TRUE)$Lineage))
#   n_lineage
#   
#   ### Differential differentiation analysis
#   if (n_lineage>1) {
#     dif_res=differentiationTest(sds, conditions=df$agecut_merged, global=TRUE, pairwise=TRUE)
#     knitr::kable(dif_res)
#     saveRDS(dif_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_differentiationTest.rds")
#   }
# }
# 
# #####################################
# 
# 
# 
# ### Plot the trajectories
# #####################################
# ### Since it shares common trajectories, skip "Combine multiple trajectories" step
# 
# ###
# df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
# sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_origin.rds")
# # load the tests results
# top_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_topologyTest.rds")
# prog_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_progressionTest.rds")
# fate_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_fateSelectionTest.rds")
# 
# df$agecut_merged=forcats::fct_relevel(df$agecut_merged, names(table(df$agecut_merged)))
# p=
#   ggplot(df, aes(x=umap_1, y=umap_2, col=agecut_merged, alpha=agecut_merged)) +
#   ggrastr::geom_point_rast(size=0.01, shape=20,
#                            color=paletteer::paletteer_d("ggsci::category20_d3")[2]) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
#   scale_alpha_manual(values=seq(0.1,0.2,length.out=length(unique(df$agecut_merged)))) +
#   theme_light() +
#   labs(x="UMAP_1", y="UMAP_2", title="CD4T cells",
#        subtitle=paste0("Progression: p=",ifelse(prog_res$p.value[1]<0.001, sprintf("%0.01e",prog_res$p.value[1]), round(prog_res$p.value[1],3)),
#                        ";    Fate selection: p=",ifelse(fate_res$p.value[1]<0.001, sprintf("%0.01e",fate_res$p.value[1]), round(fate_res$p.value[1],3)))) +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         plot.subtitle=element_text(size=10),
#         title=element_text(size=10),
#         legend.position="bottom",
#         legend.title=element_text(size=10),
#         legend.text=element_text(size=9)) +
#   guides(alpha="none")
# 
# midpoint=df %>% group_by(Annot.detailed, agecut_merged) %>%
#   dplyr::summarise(umap_1=mean(umap_1), umap_2=mean(umap_2), .groups=NULL)
# 
# edges=lapply(slingLineages(sds_merged), function(lin){
#   from=lin[1:(length(lin)-1)]
#   to=lin[2:length(lin)]
#   return(data.frame("from"=from, "to"=to))}
# ) %>%
#   bind_rows()
# for (batch in unique(df$agecut_merged)) {
#   cl_batch=midpoint %>% filter(agecut_merged==batch)
#   edges_batch <- left_join(edges, cl_batch, by = c("from"="Annot.detailed")) %>%
#     left_join(cl_batch %>%
#                 dplyr::rename("umap1_end"="umap_1", "umap2_end"="umap_2") %>%
#                 select(-agecut_merged),
#               by=c("to"="Annot.detailed"))
#   p=
#     p +
#     geom_segment(data=edges_batch, inherit.aes=FALSE,
#                  aes(x=umap_1, y=umap_2, xend=umap1_end, yend=umap2_end, col=agecut_merged),
#                  alpha=1,
#                  lineend="round", linejoin="round",
#                  size=0.5)
# }
# 
# p_withLabel=
#   p +
#   ggrepel::geom_label_repel(data=midpoint %>% group_by(Annot.detailed) %>% summarize_at(c("umap_1","umap_2"), mean),
#                             inherit.aes=FALSE, force_pull=100,
#                             aes(x=umap_1, y=umap_2, label=stringr::str_wrap(gsub("CD4T\\.","",gsub("_|\\."," ",Annot.detailed)), 10)),
#                             alpha=0.8,
#                             color="black") +
#   guides(color=guide_legend(title="age", title.position="left", title.hjust=0.5, nrow=1))
# 
# p_withLabel
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Pseudotime_Trajectory_CD4T.pdf", height=5, width=7)
# plot(p_withLabel)
# dev.off()
# 
# #####################################
# 
# 
# 
# ### If different conditions have different trajectories
# #####################################
# ###
# 
# 
# if (F) {
#   ### If different conditions have different trajectories
#   # get slingshot results in each condition
#   sdss=slingshot_conditions(sds, df$agecut_merged, 
#                             approx_points=FALSE,
#                             extend="n", 
#                             reweight=FALSE, reassign=FALSE)
#   sdss$condition_id=names(sdss)
#   n_condition=length(sdss$condition_id)
#   n_lineage=length(slingLineages(sdss[[1]]))
#   sdss$mapping=matrix(rep(1:n_lineage, each=n_condition), 
#                       nrow=n_lineage, ncol=n_condition, 
#                       byrow=TRUE)
#   sdss_origin=sdss
#   # merge the slingshot results
#   sds=do.call(merge_sds, sdss)
#   
#   
#   
#   # analyze differential progression
#   # * interpretation: difference in cell development/pseudotime in each lineage among different conditions
#   progressionTest(sds, conditions=df$agecut_merged, lineages=TRUE)
#   df_merged=full_join(
#     df %>% select(cellid, umap_1, umap_2, Annot.detailed, agecut_merged) %>%
#       mutate(cells=paste0("Cell-",rownames(.))),
#     slingPseudotime(sds) %>% 
#       as.data.frame() %>%
#       mutate(cells=rownames(.)),
#     by="cells"
#   ) %>%
#     tidyr::pivot_longer(starts_with("Lineage"), names_to="Curve", values_to="pst")
#   
#   ggplot(df_merged, aes(x=pst)) +
#     geom_density(alpha=0.4, aes(fill=agecut_merged), col="transparent") +
#     geom_density(aes(col=agecut_merged), fill="transparent", size=1.5) +
#     guides(col="none") +
#     scale_fill_brewer(palette="Accent") +
#     scale_color_brewer(palette="Accent") +
#     labs(x="Pseudotime", fill="Type") +
#     facet_wrap(~Curve, scales="free_x")
#   p5
#   
#   # analyze differential differentiation (i.e., differential fate selection)
#   # * interpretation: difference in cell differentiation preference to a certain lineage among different conditions
#   fateSelectionTest(sds, conditions=df$agecut_merged, pairwise=TRUE)
#   weights=condiments:::.sling_reassign(sds)
#   df_difftest=df %>%
#     mutate(cells=paste0("Cell-",rownames(.))) %>%
#     full_join(
#       weights %>% 
#         as.data.frame() %>%
#         mutate(cells=rownames(.)) %>%
#         dplyr::rename("Lineage1"=V1, "Lineage2"=V2) %>%
#         tidyr::pivot_longer(starts_with("Lineage"), names_to="Curve", values_to="weights")
#     )
#   df_w=df_difftest %>%
#     group_by(cells) %>%
#     mutate(weights=weights/sum(weights)) %>%
#     ungroup() %>%
#     group_by(agecut_merged, Curve) %>%
#     summarise(weights=mean(weights), .groups=NULL)
#   # p2=
#   ggplot(df_w, aes(x=Curve, fill=agecut_merged, y=weights)) +
#     geom_col(position="dodge") +
#     scale_fill_brewer(palette="Accent") +
#     labs(x = "", y="Mean weight")
# }
# 
# #####################################
# 
# 
# 
# ### Visualize the common trajectory
# #####################################
# ###
# 
# ### 
# df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_Umap.rds")
# sds=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_sds_origin.rds")
# conditions=unique(df$agecut_merged)
# ggplot(df, aes(x=umap_1, y=umap_2, col=match(agecut_merged, conditions))) +
#   facet_wrap(~agecut_merged) +
#   geom_point(alpha=0.5) +
#   geom_point(data=slingMST(sds, as.df=TRUE), inherit.aes=F, aes(x=umap_1, y=umap_2), size=2) +
#   geom_path(data=slingMST(sds, as.df=TRUE), inherit.aes=F, aes(x=umap_1, y=umap_2, group=Lineage), size=1.5) +
#   guides(col=FALSE)
# 
# #####################################
# 
# 
# 
# ### Analyze differential progression
# # ... * interpretation: difference in cell development/pseudotime in each lineage among different conditions
# #####################################
# ###
# 
# ### Differential progression analysis
# prog_res=progressionTest(sds, conditions=df$agecut_merged, global=TRUE, lineages=TRUE)
# knitr::kable(prog_res)
# saveRDS(prog_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_progressionTest.rds")
# 
# ### Plot
# psts=slingPseudotime(sds) %>%
#   as.data.frame() %>%
#   mutate(cells=rownames(.),
#          conditions=df$agecut_merged) %>%
#   tidyr::pivot_longer(starts_with("Lineage"), values_to="pseudotime", names_to="lineages")
# 
# ggplot(psts, aes(x=pseudotime, fill=conditions)) +
#   geom_density(alpha=0.5) +
#   facet_wrap(~lineages) +
#   theme(legend.position="bottom")
# 
# 
# 
# ### Analyze differential differentiation
# # ... * interpretation: difference in cell differentiation preference to a certain lineage among different conditions
# #####################################
# ###
# 
# ### Check the number of lineages
# n_lineage=length(unique(slingMST(sds, as.df=TRUE)$Lineage))
# n_lineage
# 
# ### Differential differentiation analysis
# if (n_lineage>1) {
#   dif_res=differentiationTest(sds, conditions=df$agecut_merged, global=TRUE, pairwise=TRUE)
#   knitr::kable(dif_res)
#   saveRDS(dif_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_DCs_differentiationTest.rds")
# }
# 
# 
# df$weight_1=slingCurveWeights(sds, as.probs=TRUE)[, 1] # there should be n-1 weights we can show, where n is the number of lineages
# ggplot(df, aes(x=weight_1, fill=agecut_merged)) +
#   geom_density(alpha=0.5) +
#   labs(x="Curve weight for the first lineage")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# weights=condiments:::.sling_reassign(sds)
# df_difftest=df %>%
#   mutate(cells=paste0("Cell-",rownames(.))) %>%
#   full_join(
#     weights %>% 
#       as.data.frame() %>%
#       mutate(cells=rownames(.)) %>%
#       # dplyr::rename("Lineage1"=V1, "Lineage2"=V2, "Lineage3"=V3) %>%
#       tidyr::pivot_longer(starts_with("Lineage"), names_to="Curve", values_to="weights")
#   )
# df_w=df_difftest %>%
#   group_by(cells) %>%
#   mutate(weights=weights/sum(weights)) %>%
#   ungroup() %>%
#   group_by(agecut_merged, Curve) %>%
#   summarise(weights=mean(weights), .groups=NULL)
# # p2=
# ggplot(df_w, aes(x=Curve, fill=agecut_merged, y=weights)) +
#   geom_col(position="dodge") +
#   scale_fill_brewer(palette="Accent") +
#   labs(x = "", y="Mean weight")
# 
# 
