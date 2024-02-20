library(dplyr)
library(readr)
library(ggpubr)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggrepel)


========================Processing of Data==================================
  
#import and add columns for conditions and run

#control 
Sample1 <- data.frame(Sample1, Run= rep("SRR3471108",nrow(Sample1)), Condition = rep("Control",nrow(Sample1)))
names(Sample1) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

Sample2 <- data.frame(Sample2, Run= rep("SRR3471109",nrow(Sample2)), Condition = rep("Control",nrow(Sample2)))
names(Sample2) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

Sample3 <- data.frame(Sample3, Run= rep("SRR3471110",nrow(Sample3)), Condition = rep("Control",nrow(Sample3)))
names(Sample3) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

#Reov_T3D_infected
Sample4 <- data.frame(Sample4, Run= rep("SRR3471111",nrow(Sample4)), Condition = rep("ReoV.Infected",nrow(Sample4)))
names(Sample4) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

Sample5 <- data.frame(Sample5, Run= rep("SRR3471112",nrow(Sample5)), Condition = rep("ReoV.Infected",nrow(Sample5)))
names(Sample5) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

Sample6 <- data.frame(Sample6, Run= rep("SRR3471113",nrow(Sample6)), Condition = rep("ReoV.Infected",nrow(Sample6)))
names(Sample6) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

#Mutant infected
Sample7 <- data.frame(Sample7, Run= rep("SRR3471114",nrow(Sample7)), Condition = rep("muReoV.Infected",nrow(Sample7)))
names(Sample7) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

Sample8 <- data.frame(Sample8, Run= rep("SRR3471115",nrow(Sample8)), Condition = rep("muReoV.Infected",nrow(Sample8)))
names(Sample8) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

Sample9 <- data.frame(Sample9, Run= rep("SRR3471116",nrow(Sample9)), Condition = rep("muReoV.Infected",nrow(Sample9)))
names(Sample9) <- c("CHROM","POS","REF","TOTAL","A","C","G","T","Run","Condition")

#Comine datasets for
#control
Control <- rbind(Sample1,Sample2,Sample3) 

#ReoV infection only
ReoV<- rbind(Sample4,Sample5,Sample6)
#Mutant ReoV only
muReoV<-rbind(Sample7,Sample8,Sample9)
#Infection
Infection <- rbind(ReoV,muReoV)
#Calculate editing frequency for each site: distinguishes between A->G and T->c

for (n in 1:nrow(Control)){
  if (Control$REF[n] == "A"){
    Control$Editing.rate[n] = Control$G[n]/Control$TOTAL[n]
  }
  else{
    Control$Editing.rate[n] = Control$C[n]/Control$TOTAL[n]
  }
} 
#ReoV
for (n in 1:nrow(ReoV)){
  if (ReoV$REF[n] == "A"){
    ReoV$Editing.rate[n] = ReoV$G[n]/ReoV$TOTAL[n]
  }
  else{
    ReoV$Editing.rate[n] = ReoV$C[n]/ReoV$TOTAL[n]
  }
} 
#muReoV
for (n in 1:nrow(muReoV)){
  if (muReoV$REF[n] == "A"){
    muReoV$Editing.rate[n] = muReoV$G[n]/muReoV$TOTAL[n]
  }
  else{
    muReoV$Editing.rate[n] = muReoV$C[n]/muReoV$TOTAL[n]
  }
} 
#Filter for significant sites
Control <- filter(Control, Editing.rate <= .99 & Editing.rate >= 0.01)
Control <- filter(Control, Editing.rate <= .49 | Editing.rate >= 0.51)
ReoV <- filter(ReoV, Editing.rate <= .99 & Editing.rate >= 0.01)
ReoV <- filter(ReoV, Editing.rate <= .49 | Editing.rate >= 0.51)
muReoV <- filter(muReoV, Editing.rate <= .99 & Editing.rate >= 0.01)
muReoV <- filter(muReoV, Editing.rate <= .49 | Editing.rate >= 0.51)

#Find mean editing rate

Controls_avg <- Control %>% group_by(POS) %>%summarise(Mean_Editing = mean(Editing.rate))


Control <- left_join(Controls_avg, 
                          Control %>% group_by(POS) %>% 
                            summarise_at(vars(-group_cols()), .funs = ~paste(unique(.), collapse ="_")) %>% 
                            ungroup()
)


ReoV_avg <- ReoV %>% group_by(POS) %>%summarise(Mean_Editing = mean(Editing.rate))


ReoV <- left_join(ReoV_avg, 
                     ReoV %>% group_by(POS) %>% 
                       summarise_at(vars(-group_cols()), .funs = ~paste(unique(.), collapse ="_")) %>% 
                       ungroup()
)

muReoV_avg <- muReoV %>% group_by(POS) %>%summarise(Mean_Editing = mean(Editing.rate))


muReoV <- left_join(muReoV_avg, 
                  muReoV %>% group_by(POS) %>% 
                    summarise_at(vars(-group_cols()), .funs = ~paste(unique(.), collapse ="_")) %>% 
                    ungroup()
)


#Whole sample staitics before  REDI annotaion
Whole_sample_statistics_non_annotated <-rbind(Control,ReoV,muReoV)
write.csv(Whole_sample_statistics_non_annotated, "Whole_sample_statistics_non_annotated.csv")

#Read REDIportal variant data
REDIportal <- read.delim("TABLE1_mm10.txt")

#Only non-dbsnp sites
REDIportal <- filter(REDIportal, dbsnp == "-")
#Add site annotations from REDIportal
names(REDIportal)[2] <- "POS"
names(REDIportal)[5] <- "Gene"
REDIportal <- REDIportal[,c(2,5,9,10,11)]
#before REDI control had 487 positions
Control <- left_join(Control,REDIportal, by='POS')
Control<- Control[complete.cases(Control$Gene), ]
#ReoV (784 postions)
ReoV <- left_join(ReoV,REDIportal, by='POS')
ReoV<- ReoV[complete.cases(ReoV$Gene), ]
#muReoV (1284)
muReoV <- left_join(muReoV,REDIportal, by='POS')
muReoV<- muReoV[complete.cases(muReoV$Gene), ]
#Annotaed sample satistics 
Whole_sample_statistics_annotated <- rbind(Control,ReoV,muReoV)
write.csv(Whole_sample_statistics_annotated, "Whole_sample_statistics_annotated.csv")

#Sites that are common in control and infection 
Common_sites_ControlvsInfection <- left_join(Control,Infection, by='POS')
Common_sites_ControlvsInfection<- Common_sites_ControlvsInfection[complete.cases(Common_sites_ControlvsInfection$Mean_Editing.y), ]
names(Common_sites_ControlvsInfection)[17] <- "Mean_Editing_Infection"

write.csv(Common_sites_ControlvsInfection, "ForCHi.csv")
#Find delta
Common_sites_ControlvsInfection <- Common_sites_ControlvsInfection %>%
  mutate( Delta= Mean_Editing_Control - Mean_Editing_Infection)


write.csv(Common_sites_ControlvsInfection, "CommonSites_ControlvsInfection.csv")
#Run T-test
Common_sites_ControlvsInfection$P_Value <- t.test(Common_sites_ControlvsInfection$Mean_Editing_Control, Common_sites_ControlvsInfection$Mean_Editing_Infection)$p.value

#SitesCommon to ReoV and muReoV
Common_sites_Infection_only <- left_join(ReoV, muReoV, by='POS')
Common_sites_Infection_only<- Common_sites_Infection_only[complete.cases(Common_sites_Infection_only$Mean_Editing.y), ]
write.csv(Common_sites_Infection_only,"CommonSites_InfectionOnly.csv")
#find delta
Common_sites_Infection_only <- Common_sites_Infection_only %>%
  mutate( Delta= Mean_Editing_ReoV - Mean_Editing_muReoV)
#Run t-TEST (Its fishy)
Common_sites_Infection_only$P_Value <- t.test(Common_sites_Infection_only$Mean_Editing.x, Common_sites_Infection_only$Mean_Editing.y)$p.value

#Find unique positions
unique_ReoV <-anti_join(ReoV,muReoV, by='POS')
Unique_muReoV <- anti_join(muReoV,ReoV,by='POS')

========================Statistics====================================================================
#ANOVA 
ANOVA_TEST$Condition<- factor(ANOVA_TEST$Condition, levels=c("Control","ReoV", "muReoV"))


baseformula <- " ~ Condition"
for (i in 2:ncol(ANOVA_TEST)) {
  formula <- paste(colnames(ANOVA_TEST)[i], baseformula, sep="")
  
  p <- summary(aov(as.formula(formula), data=ANOVA_TEST))[[1]][["Pr(>F)"]][1]
  
  print(paste(formula, ": p=", p, sep=""))
}
===============================================Plots========================================================

#Figure1
#A
TPM <- ggplot(ADAR_TPM, aes(x=Gene, y=TPM, color=Condition)) +
  geom_boxplot() +
  labs(title="Average expression of ADARs across all samples", x="Gene", y="TPM") +
  scale_color_manual(values=c("#009E73", "#E69F00", "#0072B2")) +
  theme_classic()


#B
Editing_regions <- Whole_sample_statistics_annotated %>% dplyr::count(Func.wgEncodeGencodeBasicVM16)
Editing_regions <- Editing_regions[order(Editing_regions$n,decreasing=T),]
labels = c("3' UTR","Intronic", "ncRNA_intronic", "Exonic", "Intergenic","upstream","UTR5", "ncRNA_exonic")
Editing_regions$Func.wgEncodeGencodeBasicVM16 <- labels %>% factor(levels = labels)
region_plot <- ggplot(data = Editing_regions, aes(x = "", y = n, fill = Func.wgEncodeGencodeBasicVM16)) + labs(fill="Genomic Region") + geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + theme_void() + geom_text(aes(label = n), position = position_stack(vjust = 0.5))
region_plot
#C
AEI <- ggplot(AEI, aes(x= Variant, y=Index, color= Condition)) +
  geom_boxplot() +
  scale_color_manual(values=c("#009E73", "#0072B2","#E69F00")) +
  theme_classic()

#D
MeanEditingRate <- ggplot(Whole_sample_statistics_annotated, aes(x=Condition, y=Mean_Editing, color= Condition)) +
  geom_boxplot() +
  labs( y="Mean_Editing_Rate", x="") +
  scale_color_manual(values=c("#009E73", "#0072B2","#E69F00")) +
  theme_classic() +
  theme(legend.position="none")

#Figure2

#A
A <- ggplot(Shared_Sites_Control_Infection, aes(x = Mean_Editing_Control, y = Mean_Editing_Infection, label = Gene)) +
  geom_point(color = "#E69F00") +
  geom_smooth(method = "lm", color = "#0072B2") +
  geom_text(aes(label = Gene), hjust = 1, vjust = 1, size=1.5) +
  labs(x = "Mean editing rate per site (Control)", y = "Mean editing rate per site (Infection)",
       title = "Mean editing rate for sites shared between control and infection (ReoV+muReoV)") +
  theme_classic()
#B
B <- ggplot(Shared_Sites_Control_Infection, aes(x = Gene, y = Delta)) +
  geom_col( fill = "#0072B2") +
  scale_fill_manual(values = c("#0072B2")) +
  ggtitle("Change in editing for sites shared between control and infection (ReoV+muReoV)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#C

C <- ggplot(Infection_only_shared, aes(x = Mean_Editing_ReoV, y = Mean_Editing_muReoV, label = Gene)) +
  geom_point(color ="#E69F00" ) +
  geom_smooth(method = "lm", color = "#0072B2") +
  geom_text(aes(label = Gene), hjust = 1, vjust= 1,size=1.5) +
  labs(x = "Mean editing rate per site (ReoV)", y = "Mean editing rate per site (muReoV)",
       title = "Mean editing rate for sites shared between ReoV and muReoV infection") +
  theme_classic()
#D
D <- ggplot(Infection_only_shared, aes(x = Gene, y = Delta)) +
  geom_col( fill = "#0072B2") +
  scale_fill_manual(values = c("#0072B2")) +
  ggtitle("Change in editing for sites shared between ReoV and muReoV infection") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#Combine all 4 plots for Figure2

ggarrange(A, B, C,D,
          labels = c("A", "B", "C","D"), ncol = 2, nrow = 2)

#Figure3
pathwayplot <- ggplot(KEGG_Pathways, aes(x = Percentage_of_gene_list, y = `KEGG Pathway`, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Control" = "#009E73", "ReoV-Infected" = "#0072B2", "muReoV-Infected" = "#E69F00")) +
  labs(x = "Percentage of gene list", y = "KEGG Pathway", title = "KEGG pathway representation of genes containing edited sites") +
  theme_minimal() +
  theme(axis.text.y = element_text(size=6))

#Supplemental_Figures
#Supplemental_Figures 1 & 2 were made using InteractiVenn web tool by Heberle et al., 2015.



