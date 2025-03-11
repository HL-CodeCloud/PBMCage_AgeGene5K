
### Performance on the whole unseen data (covering age from 0 to 100)
#####################################
###


### Analyze all the model predictivity on age 20-100
AgeGene5K=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_AgeGene5K_results_20TO100.txt", sep="\t")
Zhu.2023=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_Zhu.et.al.2023_results_20TO100.txt", sep="\t")
Dorum.2024=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_Dørum.2024_results_20TO100.txt", sep="\t")
Peters.2015=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_Peters.2015_results_20TO100.txt", sep="\t")

### Merge
All_DF=data.table::rbindlist(list(AgeGene5K, Zhu.2023, Dorum.2024, Peters.2015))
write.table(All_DF, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_UnseenData_20TO100_results.txt", sep="\t")

### Clean the data
All_DF=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_UnseenData_20TO100_results.txt", sep="\t")
All_DF=All_DF %>% tidyr::pivot_longer(cols=c("rsq","MAE","MAD","RMSE"), names_to="parameter", values_to="value") %>%
  subset(parameter!="rsq" & parameter!="MAD") %>%
  mutate(dataset=ifelse(grepl("Combination of ",dataset), gsub("Combination of ","",dataset), dataset))

### Use only the RNA-based model
All_DF=All_DF %>% 
  subset(model=="RNA-based") %>%
  # remove the unhealthy dataset
  subset(dataset!="GSE189050") # the SLE dataset

plot_=
  ggplot(All_DF, aes(x=dataset, y=value, fill=geneset, color=geneset)) +
  geom_bar(stat="identity", position=position_dodge(width=0.5), width=0.5) +
  facet_wrap(~parameter, scales="free_y") +
  scale_y_continuous(n.breaks=5) +
  ggsci::scale_color_d3(alpha=0.8) +
  ggsci::scale_fill_d3(alpha=0.2) +
  labs(x=NULL, y="parameter") +
  theme_light() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        strip.background=element_rect(fill="transparent", color="black"),
        strip.text=element_text(size=10, colour='black'),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10)) +
  scale_x_discrete(label=scales::label_wrap(10))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_UnseenDataPrediction_20To100.pdf", height=5, width=8)
plot(plot_)
dev.off()

#####################################



### Performance on the unseen data subset (ranging from 60 to 80 year-old)
#####################################
###

### Analyze all the model predictivity on age 60-80
AgeGene5K=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_AgeGene5K_results_60TO80.txt", sep="\t")
Zhu.2023=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_Zhu.et.al.2023_results_60TO80.txt", sep="\t")
Dorum.2024=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_Dørum.2024_results_60TO80.txt", sep="\t")
Peters.2015=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_Peters.2015_results_60TO80.txt", sep="\t")

### Merge
All_DF=data.table::rbindlist(list(AgeGene5K, Zhu.2023, Dorum.2024, Peters.2015))
write.table(All_DF, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_UnseenData_60TO80_results.txt", sep="\t")

### Clean the data
All_DF=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_UnseenData_60TO80_results.txt", sep="\t")
All_DF=All_DF %>% tidyr::pivot_longer(cols=c("rsq","MAE","MAD","RMSE"), names_to="parameter", values_to="value") %>%
  subset(parameter!="rsq" & parameter!="MAD") %>%
  mutate(dataset=ifelse(grepl("Combination of ",dataset), gsub("Combination of ","",dataset), dataset))

### Use only the RNA-based model
All_DF=All_DF %>% 
  subset(model=="RNA-based") %>%
  # remove the unhealthy dataset
  subset(dataset!="GSE189050") # the SLE dataset

plot_=
  ggplot(All_DF, aes(x=dataset, y=value, fill=geneset, color=geneset)) +
  geom_bar(stat="identity", position=position_dodge(width=0.5), width=0.5) +
  facet_wrap(~parameter, scales="free_y") +
  scale_y_continuous(n.breaks=5) +
  ggsci::scale_color_d3(alpha=0.8) +
  ggsci::scale_fill_d3(alpha=0.2) +
  labs(x=NULL, y="parameter") +
  theme_light() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        strip.background=element_rect(fill="transparent", color="black"),
        strip.text=element_text(size=10, colour='black'),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10)) +
  scale_x_discrete(label=scales::label_wrap(10))
  
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_UnseenDataPrediction_60To80.pdf", height=5, width=8)
plot(plot_)
dev.off()







