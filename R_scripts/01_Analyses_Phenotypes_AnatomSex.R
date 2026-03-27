#This script analyses the experimental data, plots figures and performs statistical tests of 
#the phenotypes and the anatomical sex

rm(list=ls())

library(dplyr)
library(ggplot2)
library(readxl)

###############  Read in dataset  ------------------------------------------------------

Data_all = read.csv("C:/Users/feiner/Dropbox/MS_Diamond_Hormones/Analyses/Data/Data_phen.csv") #233 IDs

### split by experiment
Data_est <- subset(Data_all, Experiment == "Estradiol")
Data_fad <- subset(Data_all, Experiment == "Fadrozole")

###############  Explore hatching success  ---------------------------------------------

### Estradiol
Data_est %>% group_by(Treatment) %>% summarise(count = n())
Data_est_dead <- as.data.frame(Data_est %>% 
  mutate(Dead = case_when(
           is.na(Dead)  ~ "Alive",
           Dead == "Y" ~ "Dead")) %>% 
  group_by(Dead, Treatment) %>% 
  summarise(count = n()) %>% 
  group_by(Treatment) %>% 
  mutate(percent = (count / sum(count))*100) %>% 
  dplyr::select(c(percent, count, Dead, Treatment)))

fisher.test(xtabs(count ~ Treatment + Dead, data = Data_est_dead))

### Fadrozole
Data_fad %>% group_by(Treatment) %>% summarise(count = n())
Data_fad_dead <- as.data.frame(Data_fad %>% 
                                 mutate(Dead = case_when(
                                   is.na(Dead)  ~ "Alive",
                                   Dead == "Y" ~ "Dead")) %>% 
                                 group_by(Dead, Treatment) %>% 
                                 summarise(count = n()) %>% 
                                 group_by(Treatment) %>% 
                                 mutate(percent = (count / sum(count))*100) %>% 
                                 dplyr::select(c(percent, count, Dead, Treatment)))

fisher.test(xtabs(count ~ Treatment + Dead, data = Data_fad_dead))

#Note: could also use chisq.test for estradiol, but for fadrozole observations are below 5, therefore Fisher's is recommended.

###############  Explore effect of treatment on anatomical sex close to hatching (Fig. 2B) -------

### Estradiol
Data_est_fil <- Data_est %>% filter(is.na(Dead))

# anatomical sex
Data_est_fil %>% 
  mutate(Treatment = factor(Treatment)) %>% 
  group_by(Treatment, Hemipenes) %>%  
  summarise(count = n()) %>% 
  group_by(Treatment)
pdf("C:/Users/feiner/Dropbox/MS_Diamond_Hormones/Analyses/Plots/SexRatio_estradiol.pdf", useDingbats = F)
Data_est_fil %>%
  group_by(Treatment, Hemipenes) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Treatment) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  ggplot(aes(x = "", y = prop, fill = Hemipenes)) +
  geom_col(width = 1, color = "white") +
  scale_fill_manual(values=c("#5db4e5","#cc79a7")) +
  coord_polar("y", start = 0) +
  facet_wrap(~ Treatment, ncol = 2, labeller = label_both) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom")
dev.off()

# genetic sex
Data_est_fil_sex <- as.data.frame(Data_est_fil %>% group_by(Sex_genetic, Treatment) %>% summarise(count = n()))
chisq.test(xtabs(count ~ Treatment + Sex_genetic, data = Data_est_fil_sex))

### Fadrozole
Data_fad_fil <- Data_fad %>% filter(is.na(Dead))

# anatomical sex
Data_fad_fil %>% 
  mutate(Treatment = factor(Treatment)) %>% 
  group_by(Treatment, Hemipenes) %>%  
  summarise(count = n()) %>% 
  group_by(Treatment)
pdf("C:/Users/feiner/Dropbox/MS_Diamond_Hormones/Analyses/Plots/SexRatio_fadrozole.pdf", useDingbats = F)
Data_fad_fil %>%
  group_by(Treatment, Hemipenes) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Treatment) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  ggplot(aes(x = "", y = prop, fill = Hemipenes)) +
  geom_col(width = 1, color = "white") +
  scale_fill_manual(values=c("#5db4e5","#cc79a7")) +
  coord_polar("y", start = 0) +
  facet_wrap(~ Treatment, ncol = 2, labeller = label_both) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom")
dev.off()

# genetic sex
Data_fad_fil_sex <- as.data.frame(Data_fad_fil %>% group_by(Sex_genetic, Treatment) %>% summarise(count = n()))
chisq.test(xtabs(count ~ Treatment + Sex_genetic, data = Data_fad_fil_sex))

###############  Explore effect of treatment on morph determination (Fig. 4) -------

### PREPARE DATA, Dc and DD genotypes grouped into a DD/Dc group
### Estradiol 
pattern_est <- Data_est_fil %>% 
  mutate(Sex_genetic = factor(Sex_genetic),
         Treatment = factor(Treatment),
         Dorsal_pattern = factor(Dorsal_pattern),
         Genotype = case_when(
           Genotype %in% c("Dc","DD") ~ "DD/Dc",
           Genotype == "cc" ~ "cc")) %>% 
  group_by(Treatment,Dorsal_pattern,Genotype,Sex_genetic) %>%  
  summarise(count = n()) 

### Fadrozole 
pattern_fad <- Data_fad_fil %>% 
  mutate(Sex_genetic = factor(Sex_genetic),
         Treatment = factor(Treatment),
         Dorsal_pattern = factor(Dorsal_pattern),
         Genotype = case_when(
           Genotype %in% c("Dc","DD") ~ "DD/Dc",
           Genotype == "cc" ~ "cc")) %>% 
  group_by(Treatment,Dorsal_pattern,Genotype,Sex_genetic) %>%  
  summarise(count = n()) 

pdf("C:/Users/feiner/Dropbox/MS_Diamond_Hormones/Analyses/Plots/PatternMorphs_est.pdf", useDingbats = F)
ggplot(pattern_est, aes(x = as.factor(Treatment), fill = Dorsal_pattern))+
  geom_bar(aes(y = count), stat = "identity", alpha = 1,color = NA, linewidth = 1.2) +
  ylab("Number of individuals") +
  facet_grid(Sex_genetic~    Genotype ) +
  scale_fill_manual(values=c("#BBD132","#514F21")) +     
  #ggtitle("Estradiol")+
  theme_minimal() 
dev.off()

pdf("C:/Users/feiner/Dropbox/MS_Diamond_Hormones/Analyses/Plots/PatternMorphs_fad.pdf", useDingbats = F)
ggplot(pattern_fad, aes(x = as.factor(Treatment), fill = Dorsal_pattern))+
  geom_bar(aes(y = count), stat = "identity", alpha = 1,color = NA, linewidth = 1.2) +
  ylab("Number of individuals") +
  facet_grid(Sex_genetic~    Genotype ) +
  scale_fill_manual(values=c("#BBD132","#514F21")) +     
  #ggtitle("Estradiol")+
  theme_minimal() 
dev.off()

#### STATS
# Assign match and mismatch based on fit between genotype and phenotype
### Estradiol 
pattern_est_match <- Data_est_fil %>% 
  mutate(Sex_genetic = factor(Sex_genetic),
         Treatment = factor(Treatment),
         Dorsal_pattern = factor(Dorsal_pattern),
         Genotype = case_when(
           Genotype %in% c("Dc","DD") ~ "DD/Dc",
           Genotype == "cc" ~ "cc"),
         Mismatch = case_when(
           Genotype == "DD/Dc" & Dorsal_pattern == "Diamond" ~ "Match",
           Genotype == "cc" & Dorsal_pattern == "Chevron" ~ "Match",
           TRUE ~ "Mismatch")) 

# Run Fisher's exact test for each sex separately
pattern_est_sum_F <- as.data.frame(pattern_est_match %>% filter(Sex_genetic == "F") %>% mutate(Mismatch = factor(Mismatch)) %>%
                                     group_by(Treatment, Mismatch) %>% summarise(count = n(), .groups = 'drop'))
fisher.test(xtabs(count ~ Treatment + Mismatch, data = pattern_est_sum_F))
pattern_est_sum_M <- as.data.frame(pattern_est_match %>% filter(Sex_genetic == "M") %>% mutate(Mismatch = factor(Mismatch)) %>%
                                     group_by(Treatment, Mismatch) %>% summarise(count = n(), .groups = 'drop'))
fisher.test(xtabs(count ~ Treatment + Mismatch, data = pattern_est_sum_M))

### Fadrozole 
pattern_fad_match <- Data_fad_fil %>% 
  mutate(Sex_genetic = factor(Sex_genetic),
         Treatment = factor(Treatment),
         Dorsal_pattern = factor(Dorsal_pattern),
         Genotype = case_when(
           Genotype %in% c("Dc","DD") ~ "DD/Dc",
           Genotype == "cc" ~ "cc"),
         Mismatch = case_when(
           Genotype == "DD/Dc" & Dorsal_pattern == "Diamond" ~ "Match",
           Genotype == "cc" & Dorsal_pattern == "Chevron" ~ "Match",
           TRUE ~ "Mismatch")) 

pattern_fad_sum_F <- as.data.frame(pattern_fad_match %>% filter(Sex_genetic == "F") %>% mutate(Mismatch = factor(Mismatch)) %>%
                                     group_by(Treatment, Mismatch) %>% summarise(count = n(), .groups = 'drop'))
fisher.test(xtabs(count ~ Treatment + Mismatch, data = pattern_fad_sum_F))
pattern_fad_sum_M <- as.data.frame(pattern_fad_match %>% filter(Sex_genetic == "M") %>% mutate(Mismatch = factor(Mismatch)) %>%
                                     group_by(Treatment, Mismatch) %>% summarise(count = n(), .groups = 'drop'))
fisher.test(xtabs(count ~ Treatment + Mismatch, data = pattern_fad_sum_M))

###############  Explore effect of treatment and temperature on sinuosity of diamond patterns (Fig. 5) -------
Data_est_fil_dia <- Data_est_fil %>% filter(Dorsal_pattern == "Diamond")
Data_fad_fil_dia <- Data_fad_fil %>% filter(Dorsal_pattern == "Diamond")
Data_fil_dia <- rbind(Data_est_fil_dia, Data_fad_fil_dia)

#calculate sinuosity scores
Data_fil_dia$Sinuosity <- ((Data_fil_dia$Length_L + Data_fil_dia$Length_R)/2)/Data_fil_dia$Midline

Male_dia <- Data_fil_dia %>% filter(Genotype %in% c("Dc","DD") & Sex_genetic == "M")
Male_dia$Group <- "Male_est"
Female_dia_Control <- Data_fil_dia %>% filter(Genotype %in% c("Dc","DD") & Sex_genetic == "F" & Treatment == "Control")
Female_dia_Control$Group <- "Female_Control"
Female_dia_est <- Data_fil_dia %>% filter(Genotype %in% c("Dc","DD") & Sex_genetic == "F" & Treatment == "Estradiol")
Female_dia_est$Group <- "Female_est"

#t.test(Male_dia$Sinuosity, Female_dia_Control$Sinuosity)
#t.test(Male_dia$Sinuosity, Female_dia_est$Sinuosity)
#t.test(Female_dia_Control$Sinuosity, Female_dia_est$Sinuosity)

Dia_all <- rbind(Male_dia, Female_dia_Control, Female_dia_est)

anova_model <- aov(Sinuosity ~ Group, data=Dia_all)
summary(anova_model)
TukeyHSD(anova_model)

pdf("C:/Users/feiner/Dropbox/MS_Diamond_Hormones/Analyses/Plots/Sinuosity_Sex.pdf", useDingbats = F)
ggplot(data=Dia_all, aes(x = Group, y = Sinuosity, col = Sex_genetic)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=c("#cc79a7","#5db4e5"), labels = c("Female", "Male")) + ylim(1,1.8) +
  theme_minimal() + geom_jitter(aes(), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.60), size = 1.5) 
dev.off()

#### ------------ Read in temperature data
# Read in data
Inc_exp <-  read_excel("C:/Users/feiner/Dropbox/MS_Diamond_Hormones/Analyses/Data/Summary_incubation_changed.xlsx")

# keep only individuals (diamonds that were scored)
Inc_exp <- Inc_exp %>% filter(!is.na(Midline))

#calculate sinuosity scores
Inc_exp$Sinuosity <- ((Inc_exp$Length_L + Inc_exp$Length_R)/2)/Inc_exp$Midline

# sort temperature into groups
Inc_exp <- Inc_exp %>%
  mutate(Temperature = case_when(
      Temperature %in% c("21", "22") ~ "Low",
      Temperature == "27" ~ "Mid",
      Temperature %in% c("32", "33") ~ "High"))

Inc_exp$Temperature <- factor(Inc_exp$Temperature, levels = c("Low", "Mid","High"))

pdf("C:/Users/feiner/Dropbox/MS_Diamond_Hormones/Analyses/Plots/Sinuosity_Temp.pdf", useDingbats = F)
ggplot(data=subset(Inc_exp, Genotype!="cc"), aes(x = Temperature, y = Sinuosity, col = Temperature)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values=c("#4984B8","#A47076","#FF5C34"), labels = c("Low", "Mid","High")) + ylim(1,1.8) +
  theme_minimal() + geom_jitter(aes(), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.60), size = 1.5) 
dev.off()

anova_model <- aov(Sinuosity ~ Temperature, data=subset(Inc_exp, Genotype!="cc"))
summary(anova_model)
TukeyHSD(anova_model)

#### ------------ Enrichment of genes with 8 or more EREs:
# total number of genes: 22156
# genes with 8 or more EREs: 879
# (numbers come from Natalia's ERE prediction)
(879+1)/(22156+1)


