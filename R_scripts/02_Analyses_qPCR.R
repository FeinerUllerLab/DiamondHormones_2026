#This script analyses the RT-qPCR data, plots figures and performs statistical tests

rm(list=ls())

library(dplyr)
library(ggplot2)
library(emmeans)

###############  Read in dataset  ------------------------------------------------------

### READ IN DATA FOR GENE EXPRESSION EXPERIMENT --------------------------------
qPCR_est <- read.csv("C:/Users/feiner/Dropbox/13_send_to_Nathalie/qPCR_data/Estradiol/qPCR_results_estradiol.csv")
qPCR_fad <- read.csv("C:/Users/feiner/Dropbox/13_send_to_Nathalie/qPCR_data/Aromatase/qPCR_results_aromatase.csv")

### Merge both datasets
qPCR_est$Experiment <- "Estradiol"
qPCR_fad$Experiment <- "Fadrozole"
qPCR <- rbind(qPCR_est,qPCR_fad)

# Samples with an NA value are removed
qPCR <- qPCR[!(startsWith(qPCR$Content, "Unk") & is.na(qPCR$Cq)),]

# sample 2025_133_HP was loaded incorrectly - remove
qPCR <- qPCR[!qPCR$Sample %in% "2025_133_HP",] 


### This chunk is only done once to prepare the metadata for the samples that were actually used in qPCR
#data_for_export <- qPCR
#data_for_export$Sample <- gsub("_B","__Back", data_for_export$Sample)
#data_for_export$Sample <- gsub("_HP","__Hemipenes", data_for_export$Sample)
#data_for_export <- data_for_export %>% filter(!grepl("NTC", Content, ignore.case = TRUE))
#data_for_export <- data.frame(data_for_export, do.call(rbind, strsplit(as.character(data_for_export$Sample),'__')))
#data_for_export <- unique(data_for_export[,c(12)])
#write.csv(data_for_export,"C:/Users/feiner/Dropbox/13_send_to_Nathalie/Metadata/Samples_used_qPCR.csv" )


################################################ Filter by no-template controls (NTC)
# This filters the datasets by removing any datapoint with a Cq value below the minimum NTC value

NTC_min <- qPCR %>% 
  mutate(Content = gsub("NTC.*", "NTC", Content)) %>%
  filter(!is.na(Cq) & Content == "NTC") %>% 
  group_by(index_column, Target, Experiment) %>% 
  mutate(minCq = min(Cq)) %>% 
  summarise(NTC_min = mean(minCq))

# Join the minimun Cq with the original dataset
qPCR <- right_join(qPCR, NTC_min) 

# Remove all rows that contain NTC values
qPCR <- qPCR %>% 
  mutate(Content = gsub("NTC.*", "NTC", Content)) %>%
  filter(Content != "NTC")

# Remove all rows where Cq > NTC_min
qPCR <- qPCR %>% filter(Cq < NTC_min)

################################################ Filter by standard deviations (SD)
# This filters the datasets by removing any datapoint with a Cq standard deviation above 0.5
# Create dataset that summarizes data from technical replicates

qPCR_sum <- qPCR %>% group_by(Target, Sample) %>% filter(!is.na(Cq)) %>% summarise(Replicates = n(), Cq_mean = mean(Cq), Cq_std = sd(Cq))
qPCR_sum <- qPCR_sum %>% filter(Cq_std < 0.5)

################################################ Change datasets (long to wide)

### Estradiol
qPCR_sum_wide <- qPCR_sum %>% pivot_wider(names_from = Target, values_from = c("Cq_mean", "Cq_std"))
qPCR_sum_wide <-  qPCR_sum_wide[,-c(2,6:8)]

#calculate relative expression
qPCR_sum_wide$Cq_rel_ESR1 <- 2^(qPCR_sum_wide$Cq_mean_GAPDH-qPCR_sum_wide$Cq_mean_ESR1)
qPCR_sum_wide$Cq_rel_CCDC170 <- 2^(qPCR_sum_wide$Cq_mean_GAPDH-qPCR_sum_wide$Cq_mean_CCDC170)
#keep only data on relative expression
qPCR_sum_wide <- qPCR_sum_wide[,-c(2:4)]

#rename to full tissuename
qPCR_sum_wide$Sample <- gsub("_B","__Back", qPCR_sum_wide$Sample)
qPCR_sum_wide$Sample <- gsub("_HP","__Hemipenes", qPCR_sum_wide$Sample)

# create tissue type column
qPCR_sum_wide <- data.frame(qPCR_sum_wide, do.call(rbind, strsplit(as.character(qPCR_sum_wide$Sample),'__')))
qPCR_sum_wide <- qPCR_sum_wide[,c(4,5,2,3)]
colnames(qPCR_sum_wide)[c(1,2)] <- c("Sample","Tissue")

# remove samples where both ESR1 and CCDC170 failed
qPCR_sum_wide <- qPCR_sum_wide %>% filter (!is.na(Cq_rel_CCDC170) | !is.na(Cq_rel_ESR1))


################################################ Merge with imported metadata

# Read in data
Data_exp = read.csv("C:/Users/feiner/Dropbox/13_send_to_Nathalie/Final_files/Data_qPCR.csv") #84 IDs

# Merge with other data
qPCR_sum_wide_meta <- merge(qPCR_sum_wide, Data_exp, by.x="Sample", by.y="ID", all.x=T)


################################################ Generate tissue specific datasets

qPCR_B <- subset(qPCR_sum_wide_meta, Tissue=="Back")
qPCR_Hp <- subset(qPCR_sum_wide_meta, Tissue=="Hemipenes")

Data_long_B <- pivot_longer(qPCR_B, cols = c(3,4),names_to = "Gene", values_to = "Rel_Cq")
Data_long_Hp <- pivot_longer(qPCR_Hp, cols = c(3,4),names_to = "Gene", values_to = "Rel_Cq")
Data_long <- pivot_longer(qPCR_sum_wide_meta, cols = c(3,4),names_to = "Gene", values_to = "Rel_Cq")

Data_long_B <- Data_long_B %>% filter(!is.na(Rel_Cq))
Data_long_Hp <- Data_long_Hp %>% filter(!is.na(Rel_Cq))

Data_long_B$Treatment <- factor(Data_long_B$Treatment, levels = c("Control_est", "Estradiol","Control_fad","Fadrozole"))
Data_long_Hp$Treatment <- factor(Data_long_Hp$Treatment, levels = c("Control_est", "Estradiol","Control_fad","Fadrozole"))


################################################ Generate plots

pdf("C:/Users/feiner/Dropbox/13_send_to_Nathalie/Plots/qPCR_HP_ESR1.pdf", useDingbats = F)
ggplot(subset(Data_long_Hp, Gene=="Cq_rel_ESR1"), aes(x = Treatment, col = Sex, y = Rel_Cq))+
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~ Gene)+
  labs(title = paste0("Tissue: ", "Hemipenes"))+
  scale_color_manual(values=c("#cc79a7","#5db4e5"), labels = c("Female", "Male")) + 
  theme_minimal() +
  geom_jitter(aes(), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.60), size = 1.5) 
dev.off()

pdf("C:/Users/feiner/Dropbox/13_send_to_Nathalie/Plots/qPCR_HP_CCDC170.pdf", useDingbats = F)
ggplot(subset(Data_long_Hp, Gene=="Cq_rel_CCDC170"), aes(x = Treatment, col = Sex, y = Rel_Cq))+
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~ Gene)+
  labs(title = paste0("Tissue: ", "Hemipenes"))+
  scale_color_manual(values=c("#cc79a7","#5db4e5"), labels = c("Female", "Male")) + 
  theme_minimal() +
  geom_jitter(aes(), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.60), size = 1.5) 
dev.off()

pdf("C:/Users/feiner/Dropbox/13_send_to_Nathalie/Plots/qPCR_B_ESR1.pdf", useDingbats = F)
ggplot(subset(Data_long_B, Gene=="Cq_rel_ESR1"), aes(x = Treatment, col = Sex, y = Rel_Cq))+
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~ Gene)+
  labs(title = paste0("Tissue: ", "Dorsal epithelia"))+
  scale_color_manual(values=c("#cc79a7","#5db4e5"), labels = c("Female", "Male")) + 
  theme_minimal() +
  geom_jitter(aes(), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.60), size = 1.5) 
dev.off()

pdf("C:/Users/feiner/Dropbox/13_send_to_Nathalie/Plots/qPCR_B_CCDC170.pdf", useDingbats = F)
ggplot(subset(Data_long_B, Gene=="Cq_rel_CCDC170"), aes(x = Treatment, col = Sex, y = Rel_Cq))+
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~ Gene)+
  labs(title = paste0("Tissue: ", "Dorsal epithelia"))+
  scale_color_manual(values=c("#cc79a7","#5db4e5"), labels = c("Female", "Male")) + 
  theme_minimal() +
  geom_jitter(aes(), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.60), size = 1.5) 
dev.off()

################################################ Statistical analyses

#ESR1 - HP
anova_model <- aov(Rel_Cq ~ Treatment * Sex, data = subset(Data_long_Hp, Gene=="Cq_rel_ESR1" & Experiment=="Estradiol"))
summary(anova_model)
pairs(emmeans(anova_model, ~ Sex | Treatment))

anova_model <- aov(Rel_Cq ~ Treatment * Sex, data = subset(Data_long_Hp, Gene=="Cq_rel_ESR1" & Experiment=="Fadrozole"))
summary(anova_model)
pairs(emmeans(anova_model, ~ Sex | Treatment))

#CCDC170 - HP
anova_model <- aov(Rel_Cq ~ Treatment * Sex, data = subset(Data_long_Hp, Gene=="Cq_rel_CCDC170" & Experiment=="Estradiol"))
summary(anova_model)
pairs(emmeans(anova_model, ~ Sex | Treatment))

anova_model <- aov(Rel_Cq ~ Treatment * Sex, data = subset(Data_long_Hp, Gene=="Cq_rel_CCDC170" & Experiment=="Fadrozole"))
summary(anova_model)
pairs(emmeans(anova_model, ~ Sex | Treatment))

#ESR1 - B
anova_model <- aov(Rel_Cq ~ Treatment * Sex, data = subset(Data_long_B, Gene=="Cq_rel_ESR1" & Experiment=="Estradiol"))
summary(anova_model)
pairs(emmeans(anova_model, ~ Sex | Treatment))

anova_model <- aov(Rel_Cq ~ Treatment * Sex, data = subset(Data_long_B, Gene=="Cq_rel_ESR1" & Experiment=="Fadrozole"))
summary(anova_model)
pairs(emmeans(anova_model, ~ Sex | Treatment))

#CCDC170 - B
anova_model <- aov(Rel_Cq ~ Treatment * Sex, data = subset(Data_long_B, Gene=="Cq_rel_CCDC170" & Experiment=="Estradiol"))
summary(anova_model)
pairs(emmeans(anova_model, ~ Sex | Treatment))

anova_model <- aov(Rel_Cq ~ Treatment * Sex, data = subset(Data_long_B, Gene=="Cq_rel_CCDC170" & Experiment=="Fadrozole"))
summary(anova_model)
pairs(emmeans(anova_model, ~ Sex | Treatment))

# check overall expression levels
mean_Tissue_Gene <- Data_long %>%
     filter(
         Tissue %in% c("Hemipenes", "Back"),
         Gene %in% c("Cq_rel_ESR1", "Cq_rel_CCDC170")
     ) %>%
     group_by(Tissue, Gene) %>%
     summarise(
         mean_Rel_Cq = mean(Rel_Cq, na.rm = TRUE),
         .groups = 'drop'  # Ensures groups are dropped after summarising
     ) %>%
     arrange(Tissue, Gene)  # Optional: sort by Tissue, then Gene
print(mean_Tissue_Gene)
  