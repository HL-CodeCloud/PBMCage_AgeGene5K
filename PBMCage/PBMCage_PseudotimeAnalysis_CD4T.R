library(condiments)
library(slingshot)
library(dplyr)
library(tidyr)
library(ggplot2)

# ### Prepare the obj
# #####################################
# ###
# 
# ### Load the data
# THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
# CD8T=subset(THEObj, Annot.rough=="CD4T cells")
# 
# ### Slim the obj otherwise it runs forever
# subset_=CD8T %>%
#   Seurat::NormalizeData() %>%
#   Seurat::FindVariableFeatures()
# subset_=subset_ %>%
#   Seurat::SketchData(.,
#   ncells=50000L,
#   sketched.assay="sketch",
#   method="LeverageScore",
#   var.name="leverage.score",
#   seed=2024L,
#   cast="dgCMatrix",
#   verbose=TRUE)
# assay.v5=Seurat::GetAssay(subset_, assay="sketch")
# subset_=Seurat::CreateSeuratObject(assay.v5)
# subset_[[]]=THEObj[[]][match(colnames(subset_), colnames(THEObj)),]
# 
# ### Run Umap
# CD8T=subset_ %>%
#   Seurat::NormalizeData() %>%
#   Seurat::FindVariableFeatures() %>%
#   Seurat::ScaleData() %>%
#   Seurat::RunPCA() %>%
#   Seurat::FindNeighbors(dims=1:30) %>%
#   Seurat::RunUMAP(dims=1:30)
# 
# ### Check if enough observations in each group (>20)
# table(CD8T$Annot.inter, CD8T$agecut)
# table(CD8T$Annot.detailed, CD8T$agecut)
# meta_data=CD8T[[]]
# meta_data=meta_data %>% 
#   mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.naive.COTL1pos_SOX4neg","CD4T.naive.COTL1pos_SOX4pos"), "CD4T.naive.COTL1pos", Annot.detailed)) %>%
#   mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif","CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
#   mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG","CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
#   mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Treg","CD4T.Treg.FOXP3hi"), "CD4T.Treg", Annot.detailed))
# table(meta_data$Annot.detailed, meta_data$agecut)
# CD8T$Annot.detailed=meta_data$Annot.detailed
# 
# ### Extract the umap results
# df=CD8T@reductions$umap@cell.embeddings
# df=df %>% as.data.frame() %>% mutate(cellid=rownames(.)) %>% right_join(CD8T@meta.data %>% mutate(cellid=rownames(.)), by="cellid")
# 
# # if not enough observations in some groups, then merge
# df=df %>%
#   mutate(agecut_merged=ifelse(age %in% c(19:30), "19~30",
#                               ifelse(age %in% c(31:47), "31~47",
#                                      ifelse(age %in% c(48:68), "48~68",
#                                             ifelse(age %in% c(69:97), "69~97", NA)))))
# 
# table(df$Annot.inter, df$agecut_merged)
# table(df$Annot.detailed, df$agecut_merged)
# saveRDS(df, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_Umap.rds")
# 
# ###  Plot the cells across conditions on the umap
# conditions=names(table(df$agecut_merged))
# ggplot(df, aes(x=umap_1, y=umap_2, color=match(agecut_merged, conditions))) +
#   ggrastr::geom_point_rast() +
#   scale_color_gradient(low="white",high="brown3")
# 
# #####################################



### Analyze topology difference
#####################################
###

df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_Umap.rds")

#### Visualize imbalance score
scores=imbalance_score(Object=df %>% select(umap_1, umap_2) %>% as.matrix(),
                       conditions=df$agecut_merged)
df$scores=scores$scores; df$scaled_scores=scores$scaled_scores
saveRDS(df, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_Umap.rds")

ggplot(df, aes(x=umap_1, y=umap_2, color=scaled_scores)) +
  ggrastr::geom_point_rast() +
  scale_color_gradient(low="white",high="brown3")

### Run the topology test
rd=as.matrix(df[, c("umap_1", "umap_2")])
sds=slingshot(rd, df$Annot.detailed, start.clus="CD4T.naive.HNRNPH1") # automate starts at CD4T.naive.HNRNPH1, good:)
saveRDS(sds, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_sds_origin.rds")

df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_Umap.rds")
sds=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_sds_origin.rds")
top_res=topologyTest(sds=sds, conditions=df$agecut_merged)
saveRDS(top_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_topologyTest.rds")

knitr::kable(top_res)

#####################################



### Combine multiple trajectories
#####################################
###

### Get slingshot results in each condition
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_Umap.rds")
sds=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_sds_origin.rds")
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
saveRDS(sds_merged, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_sds_merged.rds")

### Plot the trajectories
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_Umap.rds")
sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_sds_merged.rds")
# load the tests results
top_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_topologyTest.rds")
prog_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_progressionTest.rds")
fate_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_fateSelectionTest.rds")

df$agecut_merged=forcats::fct_relevel(df$agecut_merged, names(table(df$agecut_merged)))
p=
  ggplot(df, aes(x=umap_1, y=umap_2, col=agecut_merged, alpha=agecut_merged)) +
  ggrastr::geom_point_rast(size=0.01, shape=20,
                           color=paletteer::paletteer_d("ggsci::category20_d3")[2]) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  scale_alpha_manual(values=seq(0.1,0.2,length.out=length(unique(df$agecut_merged)))) +
  theme_light() +
  labs(x="UMAP_1", y="UMAP_2", title="CD4T cells",
       subtitle=paste0("Progression: p=",ifelse(prog_res$p.value[1]<0.001, sprintf("%0.01e",prog_res$p.value[1]), round(prog_res$p.value[1],3)),
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
  ggrepel::geom_label_repel(data=midpoint %>% group_by(Annot.detailed) %>% summarize_at(c("umap_1","umap_2"), mean),
                            inherit.aes=FALSE, force_pull=100,
                            aes(x=umap_1, y=umap_2, label=stringr::str_wrap(gsub("CD4T ","",gsub("_|\\."," ",Annot.detailed)), 10)),
                            alpha=0.8, box.padding=1,
                            color="black", min.segment.length=0) +
  guides(color=guide_legend(title="age", title.position="left", title.hjust=0.5, nrow=1))

p_withLabel

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Pseudotime_Trajectory_CD4T.pdf", height=5, width=7)
plot(p_withLabel)
dev.off()

#####################################



### Analyze differential progression
# ... * interpretation: difference in cell development/pseudotime in each lineage among different conditions
#####################################
###

### Differential progression analysis
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_Umap.rds")
sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_sds_merged.rds")
prog_res=progressionTest(sds_merged, conditions=df$agecut_merged, global=TRUE, lineages=TRUE)
knitr::kable(prog_res)
saveRDS(prog_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_progressionTest.rds")

### Plot the progression
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_Umap.rds")
sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_sds_merged.rds")
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
  geom_density(aes(col=agecut_merged), fill="transparent", size=1.5) +
  guides(col="none") +
  scale_fill_brewer(palette="Accent") +
  scale_color_brewer(palette="Accent") +
  labs(x="Pseudotime", fill="Type") +
  facet_wrap(~Curve, scales="free_x")

#####################################



### Analyze differential differentiation
# (Use only when common trajectory is found!!!)
# (The equivalent function in the case of different trajectories is the "fateSelectionTest" below!!!)
# ... * interpretation: difference in cell differentiation preference to a certain lineage among different conditions
#####################################
###

top_res=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_topologyTest.rds")

if (top_res$p.value>0.05) {
  sds=readRDS( "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_sds_origin.rds")

  ### Check the number of lineages
  n_lineage=length(unique(slingMST(sds, as.df=TRUE)$Lineage))
  n_lineage

  ### Differential differentiation analysis
  if (n_lineage>1) {
    dif_res=differentiationTest(sds, conditions=df$agecut_merged, global=TRUE, pairwise=TRUE)
    knitr::kable(dif_res)
    saveRDS(dif_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_differentiationTest.rds")
  }
}

#####################################



### Analyze differential fate selection
# ... * interpretation: difference in cell differentiation preference to a certain lineage among different conditions
#####################################
###

###
df=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_Umap.rds")
sds_merged=readRDS("~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_sds_merged.rds")

fate_res=fateSelectionTest(sds_merged, conditions=df$agecut_merged, global=TRUE, pairwise=TRUE)
knitr::kable(fate_res)
saveRDS(fate_res, "~/Project_PBMCage/Tempt_RDS/pseudotime_results/PBMCage_CD4T.slim_fateSelectionTest.rds")

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
#
#
#
#
#
#
# ##########################################################################
# ###### Below are for the case with common trajectory, inapplicable for CD8T here
# ##########################################################################
#
#
# ### Visualize the common trajectory
# #####################################
# ###
#
# ###
#
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
# #####################################
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
# ggplot(df, aes(x=weight_1, fill=conditions)) +
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
#       dplyr::rename("Lineage1"=V1, "Lineage2"=V2) %>%
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
