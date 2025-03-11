
### Causual analysis between expression of LCK and MYC in CD8T
#####################################
###

if (!require("CausalImpact")) install.packages("CausalImpact")
library(CausalImpact)

pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough","LCK","MYC"), layer="data") %>%
  subset(Annot.rough=="CD8T cells") %>%
  dplyr::select(age, LCK, MYC)

MYC_on_LCK=expr %>%
  mutate(LCKpos_or_neg=ifelse(LCK<=quantile(expr$LCK, 0.1), "neg", "pos")) %>%
  mutate(MYCpos_or_neg=ifelse(MYC<=quantile(expr$MYC, 0.1), "neg", "pos"))

# check if LCK is the cause
mat_pre=MYC_on_LCK %>% subset(LCKpos_or_neg=="pos") %>% select(MYC, age) %>% as.matrix()
mat_post=MYC_on_LCK %>% subset(LCKpos_or_neg=="neg") %>% select(MYC, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))

impact_LCK=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
impact_LCK$summary
plot_impact_LCK=plot(impact_LCK, "cumulative")

# check if MYC is the cause
mat_pre=MYC_on_LCK %>% subset(MYCpos_or_neg=="pos") %>% select(LCK, age) %>% as.matrix()
mat_post=MYC_on_LCK %>% subset(MYCpos_or_neg=="neg") %>% select(LCK, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))

impact_MYC=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
impact_MYC$summary
plot_impact_MYC=plot(impact_MYC, "cumulative")

# merge the results for plotting
merged_df=data.table::rbindlist(list(plot_impact_LCK$data %>% mutate(intervention="LCK"),
                                     plot_impact_MYC$data %>% mutate(intervention="MYC")))
summary_LCK=data.frame(AbsEffect=sprintf("%.2e",impact_LCK$summary$AbsEffect[1]),
                        RelEffect=sprintf("%.2f%%",impact_LCK$summary$RelEffect[1]),
                        AbsEffect.lower=sprintf("%.2e",impact_LCK$summary$AbsEffect.lower[1]),
                        AbsEffect.upper=sprintf("%.2e",impact_LCK$summary$AbsEffect.upper[1]),
                        RelEffect.lower=sprintf("%.2f%%",impact_LCK$summary$RelEffect.lower[1]),
                        RelEffect.upper=sprintf("%.2f%%",impact_LCK$summary$RelEffect.upper[1]),
                        pval=sprintf("%.3f",impact_LCK$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="LCK")

summary_MYC=data.frame(AbsEffect=sprintf("%.2e",impact_MYC$summary$AbsEffect[1]),
                       RelEffect=sprintf("%.2f%%",impact_MYC$summary$RelEffect[1]),
                       AbsEffect.lower=sprintf("%.2e",impact_MYC$summary$AbsEffect.lower[1]),
                       AbsEffect.upper=sprintf("%.2e",impact_MYC$summary$AbsEffect.upper[1]),
                       RelEffect.lower=sprintf("%.2f%%",impact_MYC$summary$RelEffect.lower[1]),
                       RelEffect.upper=sprintf("%.2f%%",impact_MYC$summary$RelEffect.upper[1]),
                       pval=sprintf("%.3f",impact_MYC$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="MYC")

summary_df_merged=rbind(summary_LCK, summary_MYC)
saveRDS(summary_df_merged, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")

# plot
library(gridExtra)
library(ggplot2)
plot_infer=
  ggplot(merged_df, aes(x=time, y=response)) +
  facet_wrap(~intervention, scales="free", ncol=1) +
  geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
  theme_classic() +
  labs(y="Cumulative effect (%)", x="metacell #") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        strip.text=element_text(size=10))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_LCK.vs.MYC.pdf", width=2.5, height=3.5)
plot_infer
dev.off()

#####################################



### Causual analysis between expression of LCK and NFKB in CD8T
#####################################
###

if (!require("CausalImpact")) install.packages("CausalImpact")
library(CausalImpact)

pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough","NFKB1","NFKB2","MAPK1","MAPK7","MAPK8"), layer="data") %>%
  subset(Annot.rough=="CD8T cells") %>%
  dplyr::select(age, NFKB1, NFKB2, MAPK1, MAPK7, MAPK8) %>%
  mutate(MAPK=(MAPK1+MAPK7+MAPK8)/3, NFKB=(NFKB1+NFKB2)/2)

NFKB_on_MAPK=expr %>%
  mutate(MAPKpos_or_neg=ifelse(MAPK<=quantile(expr$MAPK, 0.1), "neg", "pos")) %>%
  mutate(NFKBpos_or_neg=ifelse(NFKB<=quantile(expr$NFKB, 0.1), "neg", "pos"))

# check if LCK is the cause of NFKB
mat_pre=NFKB_on_MAPK %>% subset(MAPKpos_or_neg=="pos") %>% select(NFKB, age) %>% as.matrix()
mat_post=NFKB_on_MAPK %>% subset(MAPKpos_or_neg=="neg") %>% select(NFKB, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))

impact_MAPK=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
impact_MAPK$summary
plot_impact_MAPK=plot(impact_MAPK, "cumulative")

# check if NFKB is the cause
mat_pre=NFKB_on_MAPK %>% subset(NFKBpos_or_neg=="pos") %>% select(MAPK, age) %>% as.matrix()
mat_post=NFKB_on_MAPK %>% subset(NFKBpos_or_neg=="neg") %>% select(MAPK, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))

impact_NFKB=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
impact_NFKB$summary
plot_impact_NFKB=plot(impact_NFKB, "cumulative")

# merge the results for plotting
merged_df=data.table::rbindlist(list(plot_impact_LCK$data %>% mutate(intervention="LCK"),
                                     plot_impact_NFKB$data %>% mutate(intervention="NFKB")))
summary_LCK=data.frame(AbsEffect=sprintf("%.2e",impact_LCK$summary$AbsEffect[1]),
                       RelEffect=sprintf("%.2f%%",impact_LCK$summary$RelEffect[1]),
                       AbsEffect.lower=sprintf("%.2e",impact_LCK$summary$AbsEffect.lower[1]),
                       AbsEffect.upper=sprintf("%.2e",impact_LCK$summary$AbsEffect.upper[1]),
                       RelEffect.lower=sprintf("%.2f%%",impact_LCK$summary$RelEffect.lower[1]),
                       RelEffect.upper=sprintf("%.2f%%",impact_LCK$summary$RelEffect.upper[1]),
                       pval=sprintf("%.3f",impact_LCK$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="LCK")

summary_NFKB=data.frame(AbsEffect=sprintf("%.2e",impact_NFKB$summary$AbsEffect[1]),
                       RelEffect=sprintf("%.2f%%",impact_NFKB$summary$RelEffect[1]),
                       AbsEffect.lower=sprintf("%.2e",impact_NFKB$summary$AbsEffect.lower[1]),
                       AbsEffect.upper=sprintf("%.2e",impact_NFKB$summary$AbsEffect.upper[1]),
                       RelEffect.lower=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.lower[1]),
                       RelEffect.upper=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.upper[1]),
                       pval=sprintf("%.3f",impact_NFKB$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="NFKB")

summary_df_merged=rbind(summary_LCK, summary_NFKB)
summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")

# plot
library(gridExtra)
library(ggplot2)
# plot_infer=
  ggplot(merged_df, aes(x=time, y=response)) +
  facet_wrap(~intervention, scales="free", ncol=1) +
  geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
  theme_classic() +
  labs(y="Cumulative effect (%)", x="metacell #") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        strip.text=element_text(size=10))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_LCK.vs.NFKB.pdf", width=2.5, height=3.5)
plot_infer
dev.off()

#####################################



### Causual analysis between expression of MYC and cytoplasmic translation
#####################################
###

### Get the hub genes
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="purple") %>% # cytoplasmic translation
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["purple"]]@result # cytoplasmic translation
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hubs=intersect(MEs_hub1, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough","MYC",
                                               MEs_hubs),
                       layer="data") %>%
  subset(Annot.rough %in% c("CD8T cells")) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, MYC)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, MYC, target.genes) %>%
  filter(MYC!=0 & target.genes!=0)

### Plot the association
plot_cor=
  ggplot(expr, aes(x=MYC, y=target.genes)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="left", label.y.npc="top") +
  theme_minimal() +
  labs(x="Expr. of MYC", y="Avg. expr. of the hub genes", subtitle="cytoplasmic translation\n(in CD8T)") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_MYC.vs.CytoplasmicTranslation.pdf", height=2.5, width=2.5)
plot_cor
dev.off()
  
# ### Causual analysis
# if (!require("CausalImpact")) install.packages("CausalImpact")
# library(CausalImpact)
# 
# # determine the neg/pos thershold
# target.genes_on_MYC=expr %>%
#   mutate(MYCpos_or_neg=ifelse(MYC<=quantile(expr$MYC, 0.1), "neg", "pos")) %>%
#   mutate(target.genespos_or_neg=ifelse(target.genes<=quantile(expr$target.genes, 0.1), "neg", "pos"))
# 
# # check if MYC is the cause of target.genes
# mat_pre=target.genes_on_MYC %>% subset(MYCpos_or_neg=="pos") %>% select(MYC, age) %>% as.matrix()
# mat_post=target.genes_on_MYC %>% subset(MYCpos_or_neg=="neg") %>% select(MYC, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_MYC=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_MYC$summary
# plot_impact_MYC=plot(impact_MYC, "cumulative")
# 
# # check if target genes are the cause
# mat_pre=target.genes_on_MYC %>% subset(target.genespos_or_neg=="pos") %>% select(target.genes, age) %>% as.matrix()
# mat_post=target.genes_on_MYC %>% subset(target.genespos_or_neg=="neg") %>% select(target.genes, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_target.genes=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_target.genes$summary
# plot_impact_target.genes=plot(impact_target.genes, "cumulative")
# 
# # merge the results for plotting
# merged_df=data.table::rbindlist(list(plot_impact_MYC$data %>% mutate(intervention="MYC"),
#                                      plot_impact_target.genes$data %>% mutate(intervention="cytoplasmic translation")))
# summary_MYC=data.frame(AbsEffect=sprintf("%.2e",impact_MYC$summary$AbsEffect[1]),
#                        RelEffect=sprintf("%.2f%%",impact_MYC$summary$RelEffect[1]),
#                        AbsEffect.lower=sprintf("%.2e",impact_MYC$summary$AbsEffect.lower[1]),
#                        AbsEffect.upper=sprintf("%.2e",impact_MYC$summary$AbsEffect.upper[1]),
#                        RelEffect.lower=sprintf("%.2f%%",impact_MYC$summary$RelEffect.lower[1]),
#                        RelEffect.upper=sprintf("%.2f%%",impact_MYC$summary$RelEffect.upper[1]),
#                        pval=sprintf("%.3f",impact_MYC$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="MYC")
# 
# summary_target.genes=data.frame(AbsEffect=sprintf("%.2e",impact_target.genes$summary$AbsEffect[1]),
#                                 RelEffect=sprintf("%.2f%%",impact_target.genes$summary$RelEffect[1]),
#                                 AbsEffect.lower=sprintf("%.2e",impact_target.genes$summary$AbsEffect.lower[1]),
#                                 AbsEffect.upper=sprintf("%.2e",impact_target.genes$summary$AbsEffect.upper[1]),
#                                 RelEffect.lower=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.lower[1]),
#                                 RelEffect.upper=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.upper[1]),
#                                 pval=sprintf("%.3f",impact_target.genes$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="cytoplasmic translation")
# 
# summary_df_merged=rbind(summary_MYC, summary_target.genes)
# summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
# saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# 
# # plot
# library(gridExtra)
# library(ggplot2)
# merged_df$intervention=forcats::fct_relevel(merged_df$intervention, c("cytoplasmic translation","MYC"))
# # plot_infer=
# ggplot(merged_df, aes(x=time, y=response)) +
#   facet_wrap(~intervention, scales="free") +
#   geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
#   geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
#   theme_classic() +
#   labs(y="Cumulative effect (%)", x="metacell #") +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=11),
#         strip.text=element_text(size=10))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_MYC.vs.CytoplasmicTranslation.pdf", width=4, height=3)
# plot_infer
# dev.off()

#####################################



### Causual analysis between expression of MYC and tRNA modification in CD8T
#####################################
###

### Get the hub genes
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["greenyellow"]]@result # tRNA modification
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="greenyellow") %>%
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, genes_selected)[1:10] # take the top 30 hub genes from the enriched module

# MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
#   subset(module=="greenyellow") %>%
#   arrange(desc(kME)) %>%
#   dplyr::select(gene_name) %>% tibble::deframe()
# MEs_hub2=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
#   subset(module=="black") %>%
#   arrange(desc(kME)) %>%
#   dplyr::select(gene_name) %>% tibble::deframe()
# MEs_hubs=intersect(MEs_hub1[1:100], MEs_hub2[1:100])

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough","MYC",
                                               MEs_hubs),
                       layer="data") %>%
  subset(Annot.rough %in% c("CD8T cells")) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, MYC)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, MYC, target.genes) %>%
  filter(MYC!=0 & target.genes!=0)

### Plot the association
plot_cor=
  ggplot(expr, aes(x=MYC, y=target.genes)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
  theme_minimal() +
  # coord_cartesian(ylim=c(0,0.25)) +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="left", label.y.npc="top") +
  labs(x="Expr. of MYC", y="Avg. expr. of the hub genes", subtitle="tRNA modification\n(in CD8T)") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_MYC.vs.TRNAModification.pdf", height=2.5, width=2.5)
plot_cor
dev.off()

#####################################



### Causual analysis between expression of MYC and oxidative phosphorylation
#####################################
###

### Get the hub genes
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["yellow"]]@result # # oxidative phosphorylation
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="yellow") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1,genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough","MYC",
                                               MEs_hubs),
                       layer="data") %>%
  subset(Annot.rough %in% c("CD8T cells")) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, MYC)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, MYC, target.genes) %>%
  filter(MYC!=0 & target.genes!=0)

### Plot the association
plot_cor=
  ggplot(expr, aes(x=MYC, y=target.genes)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="center", label.y.npc="top") +
  theme_minimal() +
  labs(x="Expr. of MYC", y="Avg. expr. of the hub genes", subtitle="oxidative phosphorylation (in CD8T)") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_MYC.vs.OxidativePhosphorylation.pdf", height=2.5, width=4)
plot_cor
dev.off()

# ### Causual analysis
# if (!require("CausalImpact")) install.packages("CausalImpact")
# library(CausalImpact)
# 
# target.genes_on_MYC=expr %>%
#   mutate(MYCpos_or_neg=ifelse(MYC<=quantile(expr$MYC, 0.1), "neg", "pos")) %>%
#   mutate(target.genespos_or_neg=ifelse(target.genes<=quantile(expr$target.genes, 0.1), "neg", "pos"))
# 
# # check if MYC is the cause of target.genes
# mat_pre=target.genes_on_MYC %>% subset(MYCpos_or_neg=="pos") %>% select(MYC, age) %>% as.matrix()
# mat_post=target.genes_on_MYC %>% subset(MYCpos_or_neg=="neg") %>% select(MYC, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_MYC=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_MYC$summary
# plot_impact_MYC=plot(impact_MYC, "cumulative")
# 
# # check if target genes are the cause
# mat_pre=target.genes_on_MYC %>% subset(target.genespos_or_neg=="pos") %>% select(target.genes, age) %>% as.matrix()
# mat_post=target.genes_on_MYC %>% subset(target.genespos_or_neg=="neg") %>% select(target.genes, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_target.genes=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_target.genes$summary
# plot_impact_target.genes=plot(impact_target.genes, "cumulative")
# 
# # merge the results for plotting
# merged_df=data.table::rbindlist(list(plot_impact_MYC$data %>% mutate(intervention="AR"),
#                                      plot_impact_target.genes$data %>% mutate(intervention="oxidative phosphorylation")))
# summary_MYC=data.frame(AbsEffect=sprintf("%.2e",impact_MYC$summary$AbsEffect[1]),
#                        RelEffect=sprintf("%.2f%%",impact_MYC$summary$RelEffect[1]),
#                        AbsEffect.lower=sprintf("%.2e",impact_MYC$summary$AbsEffect.lower[1]),
#                        AbsEffect.upper=sprintf("%.2e",impact_MYC$summary$AbsEffect.upper[1]),
#                        RelEffect.lower=sprintf("%.2f%%",impact_MYC$summary$RelEffect.lower[1]),
#                        RelEffect.upper=sprintf("%.2f%%",impact_MYC$summary$RelEffect.upper[1]),
#                        pval=sprintf("%.3f",impact_MYC$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="AR")
# 
# summary_target.genes=data.frame(AbsEffect=sprintf("%.2e",impact_target.genes$summary$AbsEffect[1]),
#                                 RelEffect=sprintf("%.2f%%",impact_target.genes$summary$RelEffect[1]),
#                                 AbsEffect.lower=sprintf("%.2e",impact_target.genes$summary$AbsEffect.lower[1]),
#                                 AbsEffect.upper=sprintf("%.2e",impact_target.genes$summary$AbsEffect.upper[1]),
#                                 RelEffect.lower=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.lower[1]),
#                                 RelEffect.upper=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.upper[1]),
#                                 pval=sprintf("%.3f",impact_target.genes$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="oxidative phosphorylation")
# 
# summary_df_merged=rbind(summary_MYC, summary_target.genes)
# summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
# saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# 
# # plot
# library(gridExtra)
# library(ggplot2)
# merged_df$intervention=forcats::fct_relevel(merged_df$intervention, c("oxidative phosphorylation","AR"))
# plot_infer=
#   ggplot(merged_df, aes(x=time, y=response)) +
#   facet_wrap(~intervention, scales="free") +
#   geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
#   geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
#   theme_classic() +
#   labs(y="Cumulative effect (%)", x="metacell #") +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=11),
#         strip.text=element_text(size=10))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_AR.vs.OxidativePhosphorylation.pdf", width=4, height=3)
# plot_infer
# dev.off()

#####################################



### Causual analysis between expression of MYC and chromosome localization/mitotic checkpoint in CD8T
#####################################
###

go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["salmon"]]@result # chromosome localization/mitotic checkpoint
genes_selected=paste0(go_objs_subset[grepl("chromosome",go_objs_subset$Description),"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="salmon") %>% # chromosome localization/mitotic checkpoint
  arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough","MYC",
                                               MEs_hubs),
                       layer="data") %>%
  subset(Annot.rough %in% c("CD8T cells")) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, MYC)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, MYC, target.genes) %>%
  # mutate(MYC=log10(MYC), target.genes=log10(target.genes)) %>%
  filter(MYC!=0 & target.genes!=0)

### Plot the association
plot_cor=
  ggplot(expr, aes(x=MYC, y=target.genes)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="left", label.y.npc="top") +
  theme_minimal() +
  labs(x="Expr. of MYC", y="Avg. expr. of the hub genes", subtitle="chromosome localization/\nmitotic checkpoint (in CD8T)") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_MYC.vs.ChromosomeLocalization.pdf", height=2.5, width=2.5)
plot_cor
dev.off()

# ### Causual analysis
# if (!require("CausalImpact")) install.packages("CausalImpact")
# library(CausalImpact)
# 
# target.genes_on_MYC=expr %>%
#   mutate(MYCpos_or_neg=ifelse(MYC<=quantile(expr$MYC, 0.1), "neg", "pos")) %>%
#   mutate(target.genespos_or_neg=ifelse(target.genes<=quantile(expr$target.genes, 0.1), "neg", "pos"))
# 
# # check if MYC is the cause of target.genes
# mat_pre=target.genes_on_MYC %>% subset(MYCpos_or_neg=="pos") %>% select(MYC, age) %>% as.matrix()
# mat_post=target.genes_on_MYC %>% subset(MYCpos_or_neg=="neg") %>% select(MYC, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_MYC=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_MYC$summary
# plot_impact_MYC=plot(impact_MYC, "cumulative")
# 
# # check if target genes are the cause
# mat_pre=target.genes_on_MYC %>% subset(target.genespos_or_neg=="pos") %>% select(target.genes, age) %>% as.matrix()
# mat_post=target.genes_on_MYC %>% subset(target.genespos_or_neg=="neg") %>% select(target.genes, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_target.genes=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_target.genes$summary
# plot_impact_target.genes=plot(impact_target.genes, "cumulative")
# 
# # merge the results for plotting
# merged_df=data.table::rbindlist(list(plot_impact_MYC$data %>% mutate(intervention="MYC"),
#                                      plot_impact_target.genes$data %>% mutate(intervention="chromosome localization")))
# summary_MYC=data.frame(AbsEffect=sprintf("%.2e",impact_MYC$summary$AbsEffect[1]),
#                        RelEffect=sprintf("%.2f%%",impact_MYC$summary$RelEffect[1]),
#                        AbsEffect.lower=sprintf("%.2e",impact_MYC$summary$AbsEffect.lower[1]),
#                        AbsEffect.upper=sprintf("%.2e",impact_MYC$summary$AbsEffect.upper[1]),
#                        RelEffect.lower=sprintf("%.2f%%",impact_MYC$summary$RelEffect.lower[1]),
#                        RelEffect.upper=sprintf("%.2f%%",impact_MYC$summary$RelEffect.upper[1]),
#                        pval=sprintf("%.3f",impact_MYC$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="MYC")
# 
# summary_target.genes=data.frame(AbsEffect=sprintf("%.2e",impact_target.genes$summary$AbsEffect[1]),
#                                 RelEffect=sprintf("%.2f%%",impact_target.genes$summary$RelEffect[1]),
#                                 AbsEffect.lower=sprintf("%.2e",impact_target.genes$summary$AbsEffect.lower[1]),
#                                 AbsEffect.upper=sprintf("%.2e",impact_target.genes$summary$AbsEffect.upper[1]),
#                                 RelEffect.lower=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.lower[1]),
#                                 RelEffect.upper=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.upper[1]),
#                                 pval=sprintf("%.3f",impact_target.genes$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="chromosome localization")
# 
# summary_df_merged=rbind(summary_MYC, summary_target.genes)
# summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
# saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# 
# # plot
# library(gridExtra)
# library(ggplot2)
# merged_df$intervention=forcats::fct_relevel(merged_df$intervention, c("chromosome localization","MYC"))
# plot_infer=
#   ggplot(merged_df, aes(x=time, y=response)) +
#   facet_wrap(~intervention, scales="free") +
#   geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
#   geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
#   theme_classic() +
#   labs(y="Cumulative effect (%)", x="metacell #") +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=11),
#         strip.text=element_text(size=10))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_MYC.vs.ChromosomeLocalization.pdf", width=4, height=3)
# plot_infer
# dev.off()

#####################################



### Causual analysis between expression of MYC and defense response
#####################################
###

go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["turquoise"]]@result # defense response
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="turquoise") %>% # defense response
  arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

# MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
#   subset(module=="turquoise") %>% # defense response
#   arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
# MEs_hub2=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
#   subset(module=="turquoise") %>% # defense response
#   arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
# MEs_hubs=intersect(MEs_hub1[1:30], MEs_hub2[1:30])

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough",
                                               "NFKB2","NFKB1","RELA","RELB","REL",
                                               MEs_hubs), 
                       layer="data") %>%
  subset(Annot.rough=="CD8T cells") %>%
  mutate(NFKB=(NFKB2+NFKB1+RELA+RELB+REL)/5) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough,NFKB2,NFKB1,RELA,RELB,REL))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, NFKB)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, NFKB, target.genes) %>%
  filter(NFKB!=0 & target.genes!=0)

### Plot the association
plot_cor=
  ggplot(expr, aes(x=NFKB, y=target.genes)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="left", label.y.npc="top") +
  theme_minimal() +
  # coord_cartesian(ylim=c(0,4)) +
  labs(x="Avg. expr. of NFKB", y="Avg. expr. of the hub genes", subtitle="defense response\n(in CD8T)") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_NFKB.vs.DefenseResponse.pdf", height=2.5, width=2.5)
plot_cor
dev.off()

### Causual analysis
if (!require("CausalImpact")) install.packages("CausalImpact")
library(CausalImpact)

target.genes_on_NFKB=expr %>%
  mutate(NFKBpos_or_neg=ifelse(NFKB<=quantile(expr$NFKB, 0.1), "neg", "pos")) %>%
  mutate(target.genespos_or_neg=ifelse(target.genes<=quantile(expr$target.genes, 0.1), "neg", "pos"))

# check if NFKB is the cause of target.genes
mat_pre=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="pos") %>% select(NFKB, age) %>% as.matrix()
mat_post=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="neg") %>% select(NFKB, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))

impact_NFKB=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
impact_NFKB$summary
plot_impact_NFKB=plot(impact_NFKB, "cumulative")

# check if NFKB1 is the cause
mat_pre=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="pos") %>% select(target.genes, age) %>% as.matrix()
mat_post=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="neg") %>% select(target.genes, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))

impact_target.genes=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
impact_target.genes$summary
plot_impact_target.genes=plot(impact_target.genes, "cumulative")

# merge the results for plotting
merged_df=data.table::rbindlist(list(plot_impact_NFKB$data %>% mutate(intervention="NFKB"),
                                     plot_impact_target.genes$data %>% mutate(intervention="defense response")))
summary_NFKB=data.frame(AbsEffect=sprintf("%.2e",impact_NFKB$summary$AbsEffect[1]),
                       RelEffect=sprintf("%.2f%%",impact_NFKB$summary$RelEffect[1]),
                       AbsEffect.lower=sprintf("%.2e",impact_NFKB$summary$AbsEffect.lower[1]),
                       AbsEffect.upper=sprintf("%.2e",impact_NFKB$summary$AbsEffect.upper[1]),
                       RelEffect.lower=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.lower[1]),
                       RelEffect.upper=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.upper[1]),
                       pval=sprintf("%.3f",impact_NFKB$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="NFKB")

summary_target.genes=data.frame(AbsEffect=sprintf("%.2e",impact_target.genes$summary$AbsEffect[1]),
                        RelEffect=sprintf("%.2f%%",impact_target.genes$summary$RelEffect[1]),
                        AbsEffect.lower=sprintf("%.2e",impact_target.genes$summary$AbsEffect.lower[1]),
                        AbsEffect.upper=sprintf("%.2e",impact_target.genes$summary$AbsEffect.upper[1]),
                        RelEffect.lower=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.lower[1]),
                        RelEffect.upper=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.upper[1]),
                        pval=sprintf("%.3f",impact_target.genes$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="defense response")

summary_df_merged=rbind(summary_NFKB, summary_target.genes)
summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")

# plot
library(gridExtra)
library(ggplot2)
plot_infer=
  ggplot(merged_df, aes(x=time, y=response)) +
  facet_wrap(~intervention, scales="free") +
  geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
  theme_classic() +
  labs(y="Cumulative effect (%)", x="metacell #") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        strip.text=element_text(size=10))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_NFKB.vs.DefenseResponse.pdf", width=4, height=3)
plot_infer
dev.off()

#####################################



### Causual analysis between expression of NFKB and humoral immunity
#####################################
###

go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["blue"]]@result # humoral immunity
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="blue") %>% # humoral immunity
  arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough",
                                               "NFKB2","NFKB1","RELA","RELB","REL",
                                               MEs_hubs),
                       layer="data") %>%
  subset(Annot.rough=="B cells") %>%
  mutate(NFKB=(NFKB2+NFKB1+RELA+RELB+REL)/5) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough,NFKB2,NFKB1,RELA,RELB,REL))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, NFKB)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, NFKB, target.genes) %>%
  filter(NFKB!=0 & target.genes!=0)

### Plot the association
plot_cor=
  ggplot(expr, aes(x=NFKB, y=target.genes)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
  scale_color_manual(values=c("black","grey50")) +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="left", label.y.npc="top") +
  theme_minimal() +
  labs(x="Avg. expr. of NFKB", y="Avg. expr. of the hub genes", subtitle="humoral immunity\n(in B cells)") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_NFKB.vs.HumoralImmunity.pdf", height=2.5, width=2.5)
plot_cor
dev.off()

# ### Causual analysis
# if (!require("CausalImpact")) install.packages("CausalImpact")
# library(CausalImpact)
# 
# target.genes_on_NFKB=expr %>%
#   mutate(NFKBpos_or_neg=ifelse(NFKB<=quantile(expr$NFKB, 0.1), "neg", "pos")) %>%
#   mutate(target.genespos_or_neg=ifelse(target.genes<=quantile(expr$target.genes, 0.1), "neg", "pos"))
# 
# # check if NFKB is the cause of target.genes
# mat_pre=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="pos") %>% select(NFKB, age) %>% as.matrix()
# mat_post=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="neg") %>% select(NFKB, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_NFKB=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_NFKB$summary
# plot_impact_NFKB=plot(impact_NFKB, "cumulative")
# 
# # check if NFKB1 is the cause
# mat_pre=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="pos") %>% select(target.genes, age) %>% as.matrix()
# mat_post=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="neg") %>% select(target.genes, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_target.genes=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_target.genes$summary
# plot_impact_target.genes=plot(impact_target.genes, "cumulative")
# 
# # merge the results for plotting
# merged_df=data.table::rbindlist(list(plot_impact_NFKB$data %>% mutate(intervention="NFKB"),
#                                      plot_impact_target.genes$data %>% mutate(intervention="humoral immunity")))
# summary_NFKB=data.frame(AbsEffect=sprintf("%.2e",impact_NFKB$summary$AbsEffect[1]),
#                         RelEffect=sprintf("%.2f%%",impact_NFKB$summary$RelEffect[1]),
#                         AbsEffect.lower=sprintf("%.2e",impact_NFKB$summary$AbsEffect.lower[1]),
#                         AbsEffect.upper=sprintf("%.2e",impact_NFKB$summary$AbsEffect.upper[1]),
#                         RelEffect.lower=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.lower[1]),
#                         RelEffect.upper=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.upper[1]),
#                         pval=sprintf("%.3f",impact_NFKB$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="NFKB")
# 
# summary_target.genes=data.frame(AbsEffect=sprintf("%.2e",impact_target.genes$summary$AbsEffect[1]),
#                                 RelEffect=sprintf("%.2f%%",impact_target.genes$summary$RelEffect[1]),
#                                 AbsEffect.lower=sprintf("%.2e",impact_target.genes$summary$AbsEffect.lower[1]),
#                                 AbsEffect.upper=sprintf("%.2e",impact_target.genes$summary$AbsEffect.upper[1]),
#                                 RelEffect.lower=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.lower[1]),
#                                 RelEffect.upper=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.upper[1]),
#                                 pval=sprintf("%.3f",impact_target.genes$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="humoral immunity")
# 
# summary_df_merged=rbind(summary_NFKB, summary_target.genes)
# summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
# saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# 
# # plot
# library(gridExtra)
# library(ggplot2)
# plot_infer=
#   ggplot(merged_df, aes(x=time, y=response)) +
#   facet_wrap(~intervention, scales="free") +
#   geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
#   geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
#   theme_classic() +
#   labs(y="Cumulative effect (%)", x="metacell #") +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=11),
#         strip.text=element_text(size=10))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_NFKB.vs.HumoralImmunity.pdf", width=4, height=3)
# plot_infer
# dev.off()

#####################################



### Causual analysis between expression of NFKB and T cell differentiation
#####################################
###

go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["green"]]@result # T cell differentiation
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="green") %>% # T cell differentiation
  arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough",
                                               "NFKB2","NFKB1","RELA","RELB","REL",
                                               MEs_hubs),
                       layer="data") %>%
  subset(Annot.rough=="CD8T cells") %>%
  mutate(NFKB=(NFKB2+NFKB1+RELA+RELB+REL)/5) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough,NFKB2,NFKB1,RELA,RELB,REL))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, NFKB)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, NFKB, target.genes) %>%
  filter(NFKB!=0 & target.genes!=0)
### Plot the association
plot_cor=
  ggplot(expr, aes(x=NFKB, y=target.genes)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="left", label.y.npc="top") +
  theme_minimal() +
  labs(x="Avg. expr. of NFKB", y="Avg. expr. of the hub genes", subtitle="regulation of T cell\ndifferentiation (in CD8T)") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_NFKB.vs.TCellDifferentiation.pdf", height=2.5, width=2.5)
plot_cor
dev.off()

# ### Causual analysis
# if (!require("CausalImpact")) install.packages("CausalImpact")
# library(CausalImpact)
# 
# target.genes_on_NFKB=expr %>%
#   mutate(NFKBpos_or_neg=ifelse(NFKB<=quantile(expr$NFKB, 0.1), "neg", "pos")) %>%
#   mutate(target.genespos_or_neg=ifelse(target.genes<=quantile(expr$target.genes, 0.1), "neg", "pos"))
# 
# # check if NFKB is the cause of target.genes
# mat_pre=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="pos") %>% select(NFKB, age) %>% as.matrix()
# mat_post=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="neg") %>% select(NFKB, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_NFKB=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_NFKB$summary
# plot_impact_NFKB=plot(impact_NFKB, "cumulative")
# 
# # check if NFKB1 is the cause
# mat_pre=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="pos") %>% select(target.genes, age) %>% as.matrix()
# mat_post=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="neg") %>% select(target.genes, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_target.genes=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_target.genes$summary
# plot_impact_target.genes=plot(impact_target.genes, "cumulative")
# 
# # merge the results for plotting
# merged_df=data.table::rbindlist(list(plot_impact_NFKB$data %>% mutate(intervention="NFKB"),
#                                      plot_impact_target.genes$data %>% mutate(intervention="T cell differentiation")))
# summary_NFKB=data.frame(AbsEffect=sprintf("%.2e",impact_NFKB$summary$AbsEffect[1]),
#                         RelEffect=sprintf("%.2f%%",impact_NFKB$summary$RelEffect[1]),
#                         AbsEffect.lower=sprintf("%.2e",impact_NFKB$summary$AbsEffect.lower[1]),
#                         AbsEffect.upper=sprintf("%.2e",impact_NFKB$summary$AbsEffect.upper[1]),
#                         RelEffect.lower=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.lower[1]),
#                         RelEffect.upper=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.upper[1]),
#                         pval=sprintf("%.3f",impact_NFKB$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="NFKB")
# 
# summary_target.genes=data.frame(AbsEffect=sprintf("%.2e",impact_target.genes$summary$AbsEffect[1]),
#                                 RelEffect=sprintf("%.2f%%",impact_target.genes$summary$RelEffect[1]),
#                                 AbsEffect.lower=sprintf("%.2e",impact_target.genes$summary$AbsEffect.lower[1]),
#                                 AbsEffect.upper=sprintf("%.2e",impact_target.genes$summary$AbsEffect.upper[1]),
#                                 RelEffect.lower=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.lower[1]),
#                                 RelEffect.upper=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.upper[1]),
#                                 pval=sprintf("%.3f",impact_target.genes$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="T cell differentiation")
# 
# summary_df_merged=rbind(summary_NFKB, summary_target.genes)
# summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
# saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# 
# # plot
# library(gridExtra)
# library(ggplot2)
# merged_df$intervention=forcats::fct_relevel(merged_df$intervention, c("T cell differentiation","NFKB"))
# plot_infer=
#   ggplot(merged_df, aes(x=time, y=response)) +
#   facet_wrap(~intervention, scales="free") +
#   geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
#   geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
#   theme_classic() +
#   labs(y="Cumulative effect (%)", x="metacell #") +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=11),
#         strip.text=element_text(size=10))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_NFKB.vs.TCellDifferentiation.pdf", width=4, height=3)
# plot_infer
# dev.off()

#####################################



### Causual analysis between expression of NFKB and leukocyte-mediated cytotoxicity
#####################################
###

go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["brown"]]@result # leukocyte-mediated cytotoxicity
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="brown") %>% # leukocyte-mediated cytotoxicity
  arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough",
                                               "NFKB2","NFKB1","RELA","RELB","REL",
                                               MEs_hubs), 
                       layer="data") %>%
  subset(Annot.rough=="CD8T cells") %>%
  mutate(NFKB=(NFKB2+NFKB1+RELA+RELB+REL)/5) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough,NFKB2,NFKB1,RELA,RELB,REL))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, NFKB)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, NFKB, target.genes) %>%
  filter(NFKB!=0 & target.genes!=0)

### Plot the association
plot_cor=
  ggplot(expr, aes(x=NFKB, y=target.genes)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="left", label.y.npc="top") +
  theme_minimal() +
  labs(x="Avg. expr. of NFKB", y="Avg. expr. of the hub genes", subtitle="leukocyte-mediated\ncytotoxicity (in CD8T)") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_NFKB.vs.LeukocyteMediatedCytotoxicity.pdf", height=2.5, width=2.5)
plot_cor
dev.off()

# ### Causual analysis
# if (!require("CausalImpact")) install.packages("CausalImpact")
# library(CausalImpact)
# 
# target.genes_on_NFKB=expr %>%
#   mutate(NFKBpos_or_neg=ifelse(NFKB<=quantile(expr$NFKB, 0.1), "neg", "pos")) %>%
#   mutate(target.genespos_or_neg=ifelse(target.genes<=quantile(expr$target.genes, 0.1), "neg", "pos"))
# 
# # check if NFKB is the cause of target.genes
# mat_pre=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="pos") %>% select(NFKB, age) %>% as.matrix()
# mat_post=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="neg") %>% select(NFKB, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_NFKB=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_NFKB$summary
# plot_impact_NFKB=plot(impact_NFKB, "cumulative")
# 
# # check if NFKB1 is the cause
# mat_pre=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="pos") %>% select(target.genes, age) %>% as.matrix()
# mat_post=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="neg") %>% select(target.genes, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_target.genes=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_target.genes$summary
# plot_impact_target.genes=plot(impact_target.genes, "cumulative")
# 
# # merge the results for plotting
# merged_df=data.table::rbindlist(list(plot_impact_NFKB$data %>% mutate(intervention="NFKB"),
#                                      plot_impact_target.genes$data %>% mutate(intervention="leukocyte-mediated cytotoxicity")))
# summary_NFKB=data.frame(AbsEffect=sprintf("%.2e",impact_NFKB$summary$AbsEffect[1]),
#                         RelEffect=sprintf("%.2f%%",impact_NFKB$summary$RelEffect[1]),
#                         AbsEffect.lower=sprintf("%.2e",impact_NFKB$summary$AbsEffect.lower[1]),
#                         AbsEffect.upper=sprintf("%.2e",impact_NFKB$summary$AbsEffect.upper[1]),
#                         RelEffect.lower=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.lower[1]),
#                         RelEffect.upper=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.upper[1]),
#                         pval=sprintf("%.3f",impact_NFKB$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="NFKB")
# 
# summary_target.genes=data.frame(AbsEffect=sprintf("%.2e",impact_target.genes$summary$AbsEffect[1]),
#                                 RelEffect=sprintf("%.2f%%",impact_target.genes$summary$RelEffect[1]),
#                                 AbsEffect.lower=sprintf("%.2e",impact_target.genes$summary$AbsEffect.lower[1]),
#                                 AbsEffect.upper=sprintf("%.2e",impact_target.genes$summary$AbsEffect.upper[1]),
#                                 RelEffect.lower=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.lower[1]),
#                                 RelEffect.upper=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.upper[1]),
#                                 pval=sprintf("%.3f",impact_target.genes$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="leukocyte-mediated cytotoxicity")
# 
# summary_df_merged=rbind(summary_NFKB, summary_target.genes)
# summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
# saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# 
# # plot
# library(gridExtra)
# library(ggplot2)
# plot_infer=
#   ggplot(merged_df, aes(x=time, y=response)) +
#   facet_wrap(~intervention, scales="free") +
#   geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
#   geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
#   theme_classic() +
#   labs(y="Cumulative effect (%)", x="metacell #") +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=11),
#         strip.text=element_text(size=10))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_NFKB.vs.LeukocyteMediatedCytotoxicity.pdf", width=4, height=3)
# plot_infer
# dev.off()

#####################################



### Causual analysis between expression of MYC and regulation of hemopoiesis
#####################################
###

go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["midnightblue"]]@result # regulation of hemopoiesis
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="midnightblue") %>% # regulation of hemopoiesis
  arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough","MYC",
                                               MEs_hubs), 
                       layer="data") %>%
  subset(Annot.rough %in% c("CD8T cells")) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, MYC)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, MYC, target.genes) %>%
  filter(MYC!=0, target.genes!=0)

### Plot the association with NFKB2
plot_cor=
  ggplot(expr, aes(x=MYC, y=target.genes)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="center", label.y.npc="top") +
  theme_minimal() +
  labs(x="Expr. of MYC", y="Avg. expr. of the hub genes", subtitle="regulation of hemopoiesis (in CD8T)") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_MYC.vs.RegulationOfHemopoiesis.pdf", height=2.5, width=4)
plot_cor
dev.off()

# ### Plot the association with RELB
# expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough","RELB",
#                                                MEs_hubs), 
#                        layer="data") %>%
#   subset(Annot.rough %in% c("CD8T cells")) %>%
#   dplyr::select(-c(donor_id,sex,Annot.rough))
# target.genes=rowMeans(expr %>% dplyr::select(-c(RELB, age)), na.rm=T)
# expr$target.genes=target.genes
# expr=expr %>%
#   dplyr::select(age, RELB, target.genes) %>%
#   filter(RELB!=0, target.genes!=0)
# plot_cor=
#   ggplot(expr, aes(x=RELB, y=target.genes)) +
#   ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
#   geom_smooth(method="lm", show.legend=T, fill="grey93", linewidth=0.5, color="black") +
#   ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="center", label.y.npc="top") +
#   theme_minimal() +
#   labs(x="Expr. of RELB", y="Avg. expr. of the hub genes", subtitle="regulation of hemopoiesis") +
#   theme(plot.background=element_rect(color="transparent", fill="white"),
#         plot.subtitle=element_text(size=11),
#         axis.title=element_text(size=10),
#         axis.text=element_text(size=9),
#         legend.text=element_text(size=9),
#         legend.title=element_text(size=10),
#         strip.text=element_text(size=10)) +
#   guides(color=guide_legend(override.aes=list(fill="transparent")))

# ### Causual analysis
# if (!require("CausalImpact")) install.packages("CausalImpact")
# library(CausalImpact)
# 
# target.genes_on_NFKB=expr %>%
#   mutate(NFKBpos_or_neg=ifelse(NFKB<=quantile(expr$NFKB, 0.1), "neg", "pos")) %>%
#   mutate(target.genespos_or_neg=ifelse(target.genes<=quantile(expr$target.genes, 0.1), "neg", "pos"))
# 
# # check if NFKB is the cause of target.genes
# mat_pre=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="pos") %>% select(NFKB, age) %>% as.matrix()
# mat_post=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="neg") %>% select(NFKB, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_NFKB=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_NFKB$summary
# plot_impact_NFKB=plot(impact_NFKB, "cumulative")
# 
# # check if NFKB1 is the cause
# mat_pre=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="pos") %>% select(target.genes, age) %>% as.matrix()
# mat_post=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="neg") %>% select(target.genes, age) %>% as.matrix()
# pre.period=c(1, nrow(mat_pre))
# post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
# 
# impact_target.genes=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
# impact_target.genes$summary
# plot_impact_target.genes=plot(impact_target.genes, "cumulative")
# 
# # merge the results for plotting
# merged_df=data.table::rbindlist(list(plot_impact_NFKB$data %>% mutate(intervention="NFKB"),
#                                      plot_impact_target.genes$data %>% mutate(intervention="regulation of hemopoises")))
# summary_NFKB=data.frame(AbsEffect=sprintf("%.2e",impact_NFKB$summary$AbsEffect[1]),
#                         RelEffect=sprintf("%.2f%%",impact_NFKB$summary$RelEffect[1]),
#                         AbsEffect.lower=sprintf("%.2e",impact_NFKB$summary$AbsEffect.lower[1]),
#                         AbsEffect.upper=sprintf("%.2e",impact_NFKB$summary$AbsEffect.upper[1]),
#                         RelEffect.lower=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.lower[1]),
#                         RelEffect.upper=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.upper[1]),
#                         pval=sprintf("%.3f",impact_NFKB$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="NFKB")
# 
# summary_target.genes=data.frame(AbsEffect=sprintf("%.2e",impact_target.genes$summary$AbsEffect[1]),
#                                 RelEffect=sprintf("%.2f%%",impact_target.genes$summary$RelEffect[1]),
#                                 AbsEffect.lower=sprintf("%.2e",impact_target.genes$summary$AbsEffect.lower[1]),
#                                 AbsEffect.upper=sprintf("%.2e",impact_target.genes$summary$AbsEffect.upper[1]),
#                                 RelEffect.lower=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.lower[1]),
#                                 RelEffect.upper=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.upper[1]),
#                                 pval=sprintf("%.3f",impact_target.genes$summary$p[1])) %>%
#   mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
#          `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
#   select(Effect, `95% CI`, pval) %>%
#   mutate(Intervention="regulation of hemopoises")
# 
# summary_df_merged=rbind(summary_NFKB, summary_target.genes)
# summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
# saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
# 
# # plot
# library(gridExtra)
# library(ggplot2)
# merged_df$intervention=forcats::fct_relevel(merged_df$intervention, c("regulation of hemopoises","NFKB"))
# plot_infer=
#   ggplot(merged_df, aes(x=time, y=response)) +
#   facet_wrap(~intervention, scales="free") +
#   geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
#   geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
#   theme_classic() +
#   labs(y="Cumulative effect (%)", x="metacell #") +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=11),
#         strip.text=element_text(size=10))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_NFKB.vs.RegulationOfHemopoises.pdf", width=4, height=3)
# plot_infer
# dev.off()

#####################################



### Causual analysis between expression of NFKB and antigen processing and presentation
#####################################
###
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["black"]]@result # antigen processing & presentation
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="black") %>% # antigen processing & presentation
  arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough",
                                               "NFKB2","NFKB1","RELA","RELB","REL",
                                               MEs_hubs[!grepl("HLA-",MEs_hubs)]), # use the MHC-I related molecules in the module
                       layer="data") %>%
  subset(Annot.rough %in% c("CD8T cells")) %>%
  mutate(NFKB=(NFKB2+NFKB1+RELA+RELB+REL)/5) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough,NFKB2,NFKB1,RELA,RELB,REL))
target.genes=rowMeans(expr %>% dplyr::select(-c(NFKB, age)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, NFKB, target.genes) %>%
  filter(NFKB!=0 & target.genes!=0) %>%
  mutate(group="MHC-I\n(in CD8T)\n")
expr2=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough",
                                                "NFKB2","NFKB1","RELA","RELB","REL",
                                               MEs_hubs[grepl("HLA-",MEs_hubs)]), # use the MHC-I related molecules in the module
                       layer="data") %>%
  subset(Annot.rough %in% c("DCs","Monocytes")) %>%
  mutate(NFKB=(NFKB2+NFKB1+RELA+RELB+REL)/5) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough,NFKB2,NFKB1,RELA,RELB,REL))
target.genes=rowMeans(expr2 %>% dplyr::select(-c(NFKB, age)), na.rm=T)
expr2$target.genes=target.genes
expr2=expr2 %>%
  dplyr::select(age, NFKB, target.genes) %>%
  filter(NFKB!=0 & target.genes!=0) %>%
  mutate(group="MHC-II\n(in DCs/Monocytes)\n")
expr_total=rbind(expr, expr2)

### Plot the association
plot_cor=
  ggplot(expr_total, aes(x=NFKB, y=target.genes, color=group)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  scale_color_manual(values=c("black","grey50")) +
  geom_smooth(aes(group=group), method="lm", show.legend=T, fill="grey93", linewidth=0.5) +
  ggpubr::stat_cor(aes(group=group), method="pearson", show.legend=F, size=3.5, label.x.npc="left", label.y.npc="top") +
  theme_minimal() +
  labs(x="Avg. expr. of NFKB", y="Avg. expr. of the hub genes", subtitle="antigen processing and\npresentation") +
  theme(plot.background=element_rect(color="transparent", fill="transparent"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.position="right",
        legend.direction="vertical",
        legend.box.spacing=unit(2, "pt"),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent"), title="Submodule"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_NFKB.vs.AntigenProcessingAndPresentation.pdf", height=2.5, width=3.5)
plot_cor
dev.off()

### Causual analysis
if (!require("CausalImpact")) install.packages("CausalImpact")
library(CausalImpact)

target.genes_on_NFKB=expr %>%
  mutate(NFKBpos_or_neg=ifelse(NFKB<=quantile(expr$NFKB, 0.1), "neg", "pos")) %>%
  mutate(target.genespos_or_neg=ifelse(target.genes<=quantile(expr$target.genes, 0.1), "neg", "pos"))

# check if NFKB is the cause of target.genes
mat_pre=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="pos") %>% select(NFKB, age) %>% as.matrix()
mat_post=target.genes_on_NFKB %>% subset(NFKBpos_or_neg=="neg") %>% select(NFKB, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))

impact_NFKB=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
impact_NFKB$summary
plot_impact_NFKB=plot(impact_NFKB, "cumulative")

# check if NFKB1 is the cause
mat_pre=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="pos") %>% select(target.genes, age) %>% as.matrix()
mat_post=target.genes_on_NFKB %>% subset(target.genespos_or_neg=="neg") %>% select(target.genes, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))

impact_target.genes=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
impact_target.genes$summary
plot_impact_target.genes=plot(impact_target.genes, "cumulative")

# merge the results for plotting
merged_df=data.table::rbindlist(list(plot_impact_NFKB$data %>% mutate(intervention="NFKB"),
                                     plot_impact_target.genes$data %>% mutate(intervention="antigen processing and presentation")))
summary_NFKB=data.frame(AbsEffect=sprintf("%.2e",impact_NFKB$summary$AbsEffect[1]),
                        RelEffect=sprintf("%.2f%%",impact_NFKB$summary$RelEffect[1]),
                        AbsEffect.lower=sprintf("%.2e",impact_NFKB$summary$AbsEffect.lower[1]),
                        AbsEffect.upper=sprintf("%.2e",impact_NFKB$summary$AbsEffect.upper[1]),
                        RelEffect.lower=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.lower[1]),
                        RelEffect.upper=sprintf("%.2f%%",impact_NFKB$summary$RelEffect.upper[1]),
                        pval=sprintf("%.3f",impact_NFKB$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="NFKB")

summary_target.genes=data.frame(AbsEffect=sprintf("%.2e",impact_target.genes$summary$AbsEffect[1]),
                                RelEffect=sprintf("%.2f%%",impact_target.genes$summary$RelEffect[1]),
                                AbsEffect.lower=sprintf("%.2e",impact_target.genes$summary$AbsEffect.lower[1]),
                                AbsEffect.upper=sprintf("%.2e",impact_target.genes$summary$AbsEffect.upper[1]),
                                RelEffect.lower=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.lower[1]),
                                RelEffect.upper=sprintf("%.2f%%",impact_target.genes$summary$RelEffect.upper[1]),
                                pval=sprintf("%.3f",impact_target.genes$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="antigen processing and presentation")

summary_df_merged=rbind(summary_NFKB, summary_target.genes)
summary_df_merged_origin=readRDS("~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")
summary_df_merged_origin=rbind(summary_df_merged_origin, summary_df_merged)
saveRDS(summary_df_merged_origin, "~/Project_PBMCage/Tempt_RDS/TF.Kinase_CausualInferenceSummary.rds")

# plot
library(gridExtra)
library(ggplot2)
plot_infer=
  ggplot(merged_df, aes(x=time, y=response)) +
  facet_wrap(~intervention, scales="free") +
  geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
  theme_classic() +
  labs(y="Cumulative effect (%)", x="metacell #") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        strip.text=element_text(size=10))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_NFKB.vs.AntigenProcessingAndPresentation.pdf", width=4, height=3)
plot_infer
dev.off()

#####################################



### Causual analysis between expression of FLI1 and hemostasis
#####################################
###

# find the OtherCell-expressing genes in the module
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="red") %>% # hemostasis
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
markers_rough=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_rough.rds") %>% 
  subset(cluster=="Other cells") %>%
  dplyr::select(gene) %>% tibble::deframe()
MEs_hubs_PLT=intersect(MEs_hub1, markers_rough) # it turns out to be related to blood coagulation
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["red"]]@result
genes_selected=strsplit(paste0(go_objs_subset[grepl("coagulation",go_objs_subset$Description),"geneID"], collapse="/"), split="/")[[1]]
MEs_hubs_plt=intersect(MEs_hubs_PLT, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

# find the CD8T-expressing genes in the module
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="red") %>% # hemostasis
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
markers_rough=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_rough.rds") %>% 
  subset(cluster=="CD8T cells") %>%
  dplyr::select(gene) %>% tibble::deframe()
MEs_hubs_CD8T=intersect(MEs_hub1, markers_rough) # it turns out to be related to cell-cell adhesion
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["red"]]@result
genes_selected=strsplit(paste0(go_objs_subset[grepl("adhesion|junction|organization",go_objs_subset$Description),"geneID"], collapse="/"), split="/")[[1]]
MEs_hubs_cd8t=intersect(MEs_hubs_CD8T, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Extract the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
expr=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough",
                                               "NFKB2","NFKB1","RELA","RELB","REL",
                                               MEs_hubs_plt), 
                       layer="data") %>%
  subset(Annot.rough=="Other cells") %>%
  mutate(NFKB=(NFKB2+NFKB1+RELA+RELB+REL)/5) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough,NFKB2,NFKB1,RELA,RELB,REL))
target.genes=rowMeans(expr %>% dplyr::select(-c(age, NFKB)), na.rm=T)
expr$target.genes=target.genes
expr=expr %>%
  dplyr::select(age, NFKB, target.genes) %>%
  filter(NFKB!=0 & target.genes!=0) %>%
  mutate(group="blood coagulation\n(in platelets)\n")

expr2=Seurat::FetchData(pseudobulk_data, vars=c("donor_id","sex","age","Annot.rough",
                                                "NFKB2","NFKB1","RELA","RELB","REL",
                                                MEs_hubs_cd8t), 
                        layer="data") %>%
  subset(Annot.rough=="CD8T cells") %>%
  mutate(NFKB=(NFKB2+NFKB1+RELA+RELB+REL)/5) %>%
  dplyr::select(-c(donor_id,sex,Annot.rough,NFKB2,NFKB1,RELA,RELB,REL))
target.genes=rowMeans(expr2 %>% dplyr::select(-c(age, NFKB)), na.rm=T)
expr2$target.genes=target.genes
expr2=expr2 %>%
  dplyr::select(age, NFKB, target.genes) %>%
  filter(NFKB!=0 & target.genes!=0) %>%
  mutate(group="cell-cell junction\n(in CD8T)\n")

expr_total=rbind(expr, expr2)
expr_total$group=forcats::fct_relevel(expr_total$group, c("cell-cell junction\n(in CD8T)\n",
                                                          "blood coagulation\n(in platelets)\n"))
# plot
plot_cor=
  ggplot(expr_total, aes(x=NFKB, y=target.genes, color=group)) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  scale_color_manual(values=c("black","grey50")) +
  geom_smooth(aes(group=group), method="lm", show.legend=T, fill="grey93", linewidth=0.5) +
  ggpubr::stat_cor(method="pearson", show.legend=F, size=3.5, label.x.npc="left", label.y.npc="top") +
  theme_minimal() +
  labs(x="Avg. expr. of NFKB", y="Avg. expr. of the hub genes", subtitle="hemostasis\n") +
  theme(plot.background=element_rect(color="transparent", fill="transparent"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.position="right",
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(fill="transparent"), title="Submodule"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CorAssociation_TFs.vs.Hemostasis.pdf", height=2.5, width=4)
plot_cor
dev.off()

#####################################



### Plot the variation of the key TFs and/or kinases
#####################################
###

### Define TFs
# chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
#                 `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),NFKB=c("NFKB1","NFKB2"),KMT2A=c("KMT2A"),
#                 NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
# chosen_TFs=Reduce(c, chosen_TFs) %>% .[!duplicated(.)]
chosen_Kinases=c("JAK1","JAK3","LCK","AKT1","MAPK1","MAP3K7","MAPK8","NFKB1","NFKB2","MYC")

### Analyze the variation at the rough level in CD4 and CD8
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% subset(analysis!="all")
DEGs_Rough_selected=DEGs_Rough_subset %>%
  subset(gene %in% c(chosen_Kinases) & celltypes %in% c("CD4T cells","CD8T cells")) %>%
  mutate(gene_celltypes=paste0(gene,"_",celltypes),
         gene_analysis=paste0(gene,"_",analysis),
         gene_celltypes_analysis=paste0(gene,"_",celltypes,"_",analysis),
         celltypes_analysis=paste0(celltypes,"_",analysis))

# plot_=
  ggplot(DEGs_Rough_selected, aes(x=age, y=avg_log2FC, color=gene))+
  facet_wrap(~celltypes_analysis, ncol=4) +
  geom_point(alpha=0.25, size=0.75) +
  geom_smooth(aes(group=gene_celltypes_analysis), se=F, linewidth=0.5) +
  theme_classic() +
  # scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.direction="vertical") +
  scale_x_continuous(breaks=c(30,60,90)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))




