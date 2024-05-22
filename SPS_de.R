#This is a file for running differential expression of PEA data for the sepsis study

# Library
library(limma)
library(docstring)
library(ggrepel)
library(patchwork)

# directory 
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

# Main directory 
main_dir <- getwd()

#### Set directory ----
source('SPS_prep.R', local = TRUE)

### Download data ---- 
# KF metadata 290124
KS_meta <- read.csv("Export/KS_matched_meta.csv") 

#  metadata
df_DA_meta <- read.csv("Data/Export_meta_v2.0.2.csv") 

df_meta_pnc<-read.csv('Export/DE_meta_pnc.csv')

Olink_meta <-read.csv('Export/Olink_meta.csv')

# Get NPX Data
df_NPX_unfiltered<-read.csv('Data/Export_expr_v2.0.1.csv')

#### GET info from BINN and featseek ----

#replicated protein are "IL6"   "TNF"   "CXCL8"
replicated_olink<-c("OID20101","OID20074","OID20153","OID21237","OID21276","OID21430","OID20473","OID20563","OID20631","OID20848","OID20911","OID20997")


oid_feat_seek<-read.csv('Export/oid_feat_seek')%>%
  pull(selected)

pnc_uniprot <- read.csv('~/Documents/github/BINN_elin/data/pnc_biomarkers_new.csv', sep = ',')%>%
  colnames()

oid_pnc_binn<-Olink_meta%>%
  filter(UniProt %in% pnc_uniprot)%>%
  pull(OlinkID)

ecoli_uniprot <- read.csv('~/Documents/github/BINN_elin/data/ecoli_biomarkers_new.csv', sep = ',')%>%
  colnames()

oid_ecoli_binn<-Olink_meta%>%
  filter(UniProt %in% ecoli_uniprot)%>%
  pull(OlinkID)

sps_uniprot <- read.csv('~/Documents/github/BINN_elin/data/sps_biomarkers.csv', sep = ',')%>%
  colnames()

oid_sps_binn<-Olink_meta%>%
  filter(UniProt %in% sps_uniprot)%>%
  pull(OlinkID)

### Data wrangling ----

### NPX clean up

DAids_SPS<-
  KS_meta%>%
  filter(Till.djup.analys==1)%>%
  pull(DAid)


duplicated_proteins<- c("P01375","P05231","P10145")# note certain proteins are duplicates


sps_data<-
  df_NPX_unfiltered%>%
  filter(DAid %in% DAids_SPS)%>%
  column_to_rownames('DAid')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('OlinkID')%>%
  left_join(Olink_meta, by='OlinkID')%>%
  select(UniProt, starts_with('DA'))%>%
  filter(!UniProt %in% duplicated_proteins)%>%
  column_to_rownames('UniProt')%>%
  t()%>%
  as.data.frame()


df_NPX<-
  df_DA_meta%>% # Takes the DA meta file
  filter(Cohort=='SPS')%>% #Filter the file for Sepsis
  left_join(df_NPX_unfiltered, by = 'DAid')%>% #Joins the meta file with the expression file by DAid
  select(-c(SampleID:Origin, Cohort, Class), -all_of(replicated_olink)) #Removes some non useful columns

# Refactor metadata
sps_meta <- KS_meta%>% #Takes Raw meta file
  filter(Till.djup.analys == 1) %>% #Filter df to only include patients who meets criteria for study
  select(DAid, 
         HDA_diagnosis = HDA.diagnosis...1..Pnc...2..Mycop...3..E.coli.pye...4..S.aur.bacter...5..Strep., 
         Sepsis = Sepsis..1.ja..0.nej, 
         Alive90 = Alive.at.day.90..0..No...1..Yes.,
         Sex = Gender..1..female...2..male., 
         Age = Age.at.sampling,
         Bacterimia_pnc = If..Pnc..S..pne.in.blood.culture...0..No...1..Yes.,
         Bacterimia_ecoli = If.E.coli.pyel..was.blood.cult.pos...0..No...1..Yes.,
         Co_infection = Micro.finding..co.infection......5.d.from.presentation.in.addition.to.HDA.agent..0..No...1..Yes.,
         PI = Cohort.PI,
         BMI = BMI) %>% #Selects the parameters and renames the columns
  mutate(HDA_diagnosis = factor(HDA_diagnosis, levels = 1:5, labels = c("Pnc", "Mycop", "E.coli", "S.aur.bacter", "Strep")), 
         Sepsis = factor(Sepsis, levels = c(0,1), labels = c(FALSE, TRUE)), 
         Alive90 = factor(Alive90, levels = c(0,1), labels = c(FALSE, TRUE)),
         Sex = factor(Sex, levels =c(1,2), labels = c('Female', 'Male')),
         Age = round(Age),
         Bacterimia_pnc = na_if(Bacterimia_pnc, 0),
         Bacterimia_ecoli = na_if(Bacterimia_ecoli, 0),
         Co_infection = ifelse(Co_infection == 1, TRUE, FALSE))%>% #Changes the construct of the useful columns so the structure is neater
  unite(Bacterimia, Bacterimia_pnc:Bacterimia_ecoli, na.rm=TRUE)%>% #unites the bacterimia column
  mutate(Bacterimia = factor(Bacterimia, levels = '1', labels = TRUE))%>% #turns the bacterimia column into factor
  relocate(Alive90, .after = last_col()) #Relocates 'Alive90' to be the last column in the df


## Find DAid for NPX expression

full_df <- left_join(sps_meta, df_NPX, by ='DAid') #Create a complete df with all meta and and all NPX values in wide format

### E coli df --
e_coli_df<-
  sps_meta%>%
  filter(HDA_diagnosis == 'E.coli')%>%
  left_join(df_NPX, by = 'DAid')

### Pncc df

df_pnc <- df_meta_pnc%>%
  relocate(Alive90, .after = last_col())%>%
  inner_join(., df_NPX, by = 'DAid')



#### Method function ----
method_func <- function(df, method = NA, disease = NA) {
  if (is.na(method) && is.na(disease)) {
    df_func <- df
  } else if (method == 'featseek' && is.na(disease)) {
    df_func <- df %>%
      select(DAid:Alive90, all_of(oid_feat_seek))
  } else if (method == 'BINN' && is.na(disease)) {
    df_func <- df %>%
      select(DAid:Alive90, all_of(oid_sps_binn))
  } else if (is.na(method) && disease == 'pnc') {
    df_func <- df %>%
      filter(HDA_diagnosis == 'Pnc')
  } else if (is.na(method) && disease == 'ecoli') {
    df_func <- df %>%
      filter(HDA_diagnosis == 'E.coli')
  } else if (method == 'featseek' && disease == 'pnc') {
    df_func <- df %>%
      select(DAid:Alive90, all_of(oid_feat_seek)) %>%
      filter(HDA_diagnosis == 'Pnc')
  } else if (method == 'BINN' && disease == 'pnc') {
    df_func <- df %>%
      select(DAid:Alive90, all_of(oid_pnc_binn)) %>%
      filter(HDA_diagnosis == 'Pnc')
  } else if (method == 'BINN' && disease == 'ecoli') {
    df_func <- df %>%
      select(DAid:Alive90, all_of(oid_ecoli_binn)) %>%
      filter(HDA_diagnosis == 'E.coli')}
  else if (method == 'featseek' && disease == 'ecoli') {
    df_func <- df %>%
      select(DAid:Alive90, all_of(oid_feat_seek)) %>%
      filter(HDA_diagnosis == 'E.coli')}
  else {
    df_func <- df
  }
  
  return(df_func)
}



#### Running DE ----

de_func<- function(df, method = NA, disease = NA, column_name){
  
  #'Differential Expression Function
  #'
  #'This function generates a differantinal expression matrix.
  #'The output from this function is a tibble consisting of 'DAid,OlinkID,logFC,CI.L,CI.R,AveExpr,t, P.Value, adj.P.Val, B,UniProt,Assay,Panel,sig'.
  #'The tibble can be used to generate volcanoplot or boxplot for the most sigifican proteins
  #'
  #'@param df: The dataframe input. Could be 'e_coli_df', 'df_pnc' or 'full_df'
  #'@param column_name: The case-control variable. The column name should be a character string identical to a column name in the df. It should also be an boleanean.
  #'
  #'@examples: 
  #'de_func(df_pnc, 'Severe.pneumonia')
  #'de_func(e_coli_df, 'Sepsis')
  #'de_func(full_df, 'Alive90')
  
  
  df<-method_func(df, method = method, disease = disease)
  
  # Design a model
  design <- model.matrix(~0 + as.factor(df[[column_name]]))
  colnames(design) <- c("control", 'case')
  
  # Make contrast
  contrast <- makeContrasts(Diff = case - control, levels = design)
  
  # Fit linear model to each protein assay
  dat_fit <- 
    df%>% 
    select(!HDA_diagnosis:Alive90)%>% 
    column_to_rownames("DAid") %>% 
    t()

  fit <- lmFit(dat_fit, design = design,  method = "robust", maxit = 10000)
  
  # Apply contrast
  contrast_fit <- contrasts.fit(fit, contrast)
  
  # Apply empirical Bayes smoothing to the SE
  ebays_fit <- eBayes(contrast_fit)
  
  # Extract DE results
  DE_results <-
    topTable(ebays_fit,
             n = nrow(ebays_fit$p.value), 
             adjust.method = "fdr", 
             confint = TRUE)
  
  DE_res <- 
    DE_results %>% 
    as_tibble(rownames = "OlinkID")%>%
    left_join(Olink_meta, by = 'OlinkID')%>%
    mutate(sig = case_when(adj.P.Val < 0.05 & logFC < 0 ~ "significant down",
                           adj.P.Val < 0.05 & logFC > 0 ~ "significant up", 
                           T ~ "not significant"))
  
  return(DE_res)
  
}

### Generate plot function ---

plot_volc<- function(df, method = NA, disease = NA, column_name) {
  
  #' Volcano plot function
  #' 
  #' This function has a volcano plot as a output. 
  #' It takes a dataframe and a column name containing an boleanean variable as a input.
  #' 
  #'@param df: The dataframe input. Could be 'e_coli_df', 'df_pnc' or 'full_df'
  #'@param column_name: column_name: The case-control variable. The column name should be a character string identical to a column name in the df. It should also be an boleanean.
  #'
  #'@examples: 
  #'de_func(df_pnc, 'Severe.pneumonia')
  #'de_func(e_coli_df, 'Sepsis')
  #'de_func(full_df, 'Alive90')
  #'
  
  
  plot_df<-de_func(df, method = method, disease = disease, column_name)
  
  proteins<-nrow(plot_df)
  samples<-nrow(method_func(df, method = method, disease = disease))
  
  
  plot_df%>%
    ggplot(aes(x=logFC, y=-log10(adj.P.Val))) +
    geom_point( aes(color = sig),
                alpha = 0.7,
                show.legend = F) +
    geom_text_repel(data=plot_df%>% filter(sig != 'not significant'),
                    aes(label = Assay),
                    size=3,
                    min.segment.length =0.1,
                    max.overlaps = 7)+
    labs(caption  = paste("Proteins included:", proteins, '\nSamples:', samples))+
    plot_theme
}

# Funciton for boxplot

plot_box<- function(df, method = NA, disease = NA, column_name, most_sig = 12){
  
  #' Box plot function
  #' 
  #' This function has a box plot as a output. 
  #' It takes a dataframe and a column name containing an boleanean variable as a input.
  #' It also takes a most_sig input to decide number of plots to generate.
  #' Note that if there are multiple ad.p.val of same value all will be kept.
  #' 
  #'@param df: The dataframe input. Could be 'e_coli_df', 'df_pnc' or 'full_df'
  #'@param column_name: The case-control variable. The column name should be a character string identical to a column name in the df. It should also be an boleanean.
  #'@param most_sig: a integer input for number of most significant proteins to be ploted. Defult is 10.
  #'
  #'@examples: 
  #'de_func(df_pnc, 'Severe.pneumonia', most_sig = 5)
  #'de_func(e_coli_df, 'Sepsis')
  #'de_func(full_df, 'Alive90', most_sig = 15)
  #'
  
  
  #select the n most significanttly expressed protien by ad.p.val
  #note that if there are multiple ad.p.val of same value all will be kept.
  
  #box_df<- de_func(df, method = method, disease = disease, column_name)
  
  sig<-
    de_func(df, method = method, disease = disease, column_name)%>%
    filter(sig != 'not significant')%>%
    arrange(adj.P.Val)%>%
    slice_min(adj.P.Val, n= most_sig)%>%
    pull(Assay)
  
  df%>%
    pivot_longer(cols= -(DAid:Alive90),
                 names_to = 'OlinkID',
                 values_to = 'NPX')%>%
    left_join(Olink_meta, by = 'OlinkID')%>% #Note use of Olink Meta file
    filter(Assay %in% sig)%>%
    mutate(Assay = factor(Assay, levels = sig[0:most_sig]))%>%
    ggplot(aes(x= !!sym(column_name), y=NPX))+
    geom_point(size=0.4)+
    geom_boxplot(aes(fill= !!sym(column_name)), 
                 alpha= 0.5)+
    facet_wrap(~Assay, scales = "free", ncol=4)+
    plot_theme
  
}


### Funciton for density plots ---
plot_dens<-function(df, method = NA, disease = NA, column_name, most_sig = 12){
  
  #' Density plot function
  #' 
  #' This function has a density plot as a output. 
  #' It takes a dataframe and a column name containing an boleanean variable as a input.
  #' It also takes a most_sig input to decide number of plots to generate.
  #' Note that if there are multiple ad.p.val of same value all will be kept.
  #' 
  #'@param df: The dataframe input. Could be 'e_coli_df', 'df_pnc' or 'full_df'
  #'@param column_name: The case-control variable. The column name should be a character string identical to a column name in the df. It should also be an boleanean.
  #'@param most_sig: a integer input for number of most significant proteins to be ploted. Defult is 10.
  #'
  #'@examples: 
  #'de_func(df_pnc, 'Severe.pneumonia', most_sig = 5)
  #'de_func(e_coli_df, 'Sepsis')
  #'de_func(full_df, 'Alive90', most_sig = 15)
  #'
  
  
  sig<-
    de_func(df, method = method, disease = disease, column_name)%>%
    arrange(adj.P.Val)%>%
    slice_min(adj.P.Val, n= most_sig)%>%
    pull(Assay)
  
  df%>%
    pivot_longer(cols= -(DAid:Alive90),
                 names_to = 'OlinkID',
                 values_to = 'NPX')%>%
    left_join(Olink_meta, by = 'OlinkID')%>% #Note use of Olink Meta file
    filter(Assay %in% sig)%>%
    mutate(Assay = factor(Assay, levels = sig[0:most_sig]))%>%
    ggplot(aes(x=NPX, fill=!!sym(column_name)))+
    geom_density(alpha=0.5)+
    facet_wrap(~Assay, scales = "free", ncol=4)+
    plot_theme
  
}


#### GSEA ----

### Plots ---
volc_ecoli<- plot_volc(e_coli_df, column_name = 'Sepsis')+ggtitle('Volcano plot E coli')
volc_pnc <- plot_volc(df_pnc, column_name = 'Severe.pneumonia')+ggtitle('Volcano plot \n Streptococcus Pneumoniae')

volc_ecoli+volc_pnc+volc_both+plot_annotation( tag_levels = 'A')


box_ecoli<-plot_box(e_coli_df, column_name = 'Sepsis')
box_pnc <- plot_box(df_pnc, column_name = 'Severe.pneumonia')

volc_ecoli+box_ecoli+
  plot_annotation( tag_levels = 'A')
  
volc_pnc+box_pnc+
  plot_annotation( tag_levels = 'A')


volc_ecoli<- plot_volc(e_coli_df, method = 'featseek', column_name = 'Sepsis')
box_ecoli<-plot_box(e_coli_df, method = 'featseek', column_name = 'Sepsis')

(umap_ecoli_fs + (volc_ecoli/box_ecoli))+ plot_annotation( tag_levels = 'A')

volc_pnc <- plot_volc(df_pnc, method = 'featseek', column_name = 'Severe.pneumonia')
box_pnc <- plot_box(df_pnc, method = 'featseek', column_name = 'Severe.pneumonia')

(umap_pnc_fs+(volc_pnc/box_pnc))+
  plot_annotation( tag_levels = 'A')

ecoli<-de_func(e_coli_df, column_name = 'Sepsis')%>%
  mutate(sig = case_when(adj.P.Val < 0.05 ~ TRUE))%>%
  filter(sig == TRUE)%>%
  select(assay =Assay, adj_p_val_ecoli = adj.P.Val)

ecoli_feat<-de_func(e_coli_df, method = 'featseek', column_name = 'Sepsis')%>%
  mutate(sig = case_when(adj.P.Val < 0.05 ~ TRUE))%>%
  filter(sig == TRUE)%>%
  select(assay_feat=Assay, adj_p_val_ecoli = adj.P.Val)

intersect(ecoli$assay, ecoli_feat$assay_feat)%>%
  length()

pnc<-de_func(df_pnc, column_name = 'Severe.pneumonia')%>%
  filter(adj.P.Val < 0.05)

pnc_feat<-de_func(df_pnc, method = 'featseek', column_name = 'Severe.pneumonia')%>%
  filter(adj.P.Val < 0.05)%>%
  pull(OlinkID)

pnc_feat%>%
  filter(OlinkID %in% pnc$OlinkID)%>%
  nrow()

length(pnc_feat$assay_feat)

unique(pnc$assay %in% pnc_feat$assay_feat)

pnc<-de_func(df_pnc, column_name = 'Severe.pneumonia')%>%
  mutate(sig = case_when(adj.P.Val < 0.05 ~ TRUE))%>%
  filter(sig == TRUE)%>%
  select(Assay, adj_p_val_pnc = adj.P.Val)

pnc<-
  de_func(df_pnc, column_name = 'Severe.pneumonia')%>%
  mutate(adj_p_val_pnc = adj.P.Val)

ecoli<-
  de_func(e_coli_df, column_name = 'Sepsis')%>%
  mutate(adj_p_val_ecoli = adj.P.Val)

diff_df<-
  left_join(ecoli,pnc, by = 'Assay')%>%
  mutate(significance = case_when(adj_p_val_ecoli < 0.05 & adj_p_val_pnc < 0.05 ~ "both",
                                  adj_p_val_ecoli > 0.05 & adj_p_val_pnc > 0.05 ~ "none",
                                  adj_p_val_ecoli > 0.05 & adj_p_val_pnc < 0.05 ~ "only Streptococcus \n Pneumoniae",
                                  adj_p_val_ecoli < 0.05 & adj_p_val_pnc > 0.05 ~ "only E coli")) 

volc_both<-
  diff_df%>% 
  ggplot(aes(x=-log10(adj_p_val_ecoli), y=-log10(adj_p_val_pnc), color = significance))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")+
  geom_point()+
  geom_text_repel(aes(label=Assay), show.legend = F)+
  labs(title = 'Volcano plot both')+
  plot_theme

oid_pnc_binn_box<-c("PSMA1","NTF3","NEFL","IL6R","ARNT","ANGPTL3","NRP1","ROBO2","NOS1","IL2")
oid_ecoli_binn_box<-c('RASA1',
                      'TBL1X',
                      'CLPS',
                      'PSMA1',
                      'FGF2',
                      'IL6R',
                      'PPY',
                      'TNFSF10',
                      'FABP6',
                      'TDGF1')

df_pnc%>%
  pivot_longer(cols= -(DAid:Alive90),
               names_to = 'OlinkID',
               values_to = 'NPX')%>%
  left_join(Olink_meta, by = 'OlinkID')%>%
  filter(Assay %in% oid_pnc_binn_box)%>%
  mutate(Assay = factor(Assay, levels = oid_ecoli_binn_box[0:10]))%>%
  ggplot(aes(x= Severe.pneumonia, y=NPX))+
  geom_point(size=0.4)+
  geom_boxplot(aes(fill= Severe.pneumonia), 
               alpha= 0.5)+
  facet_wrap(~Assay, scales = "free")+
  plot_theme

a<-e_coli_df%>%
  pivot_longer(cols= -(DAid:Alive90),
               names_to = 'OlinkID',
               values_to = 'NPX')%>%
  left_join(Olink_meta, by = 'OlinkID')%>%
  filter(Assay %in% oid_ecoli_binn_box)%>%
  mutate(Assay = factor(Assay, levels = oid_ecoli_binn_box[0:10]))%>%
  ggplot(aes(x= Sepsis, y=NPX))+
  geom_point(size=0.4)+
  geom_boxplot(aes(fill= Sepsis), 
               alpha= 0.5)+
  facet_wrap(~Assay, scales = "free")+
  plot_theme


plot_df<-de_func(e_coli_df, column_name = 'Sepsis')

proteins<-nrow(plot_df)
samples<-nrow(method_func(e_coli_df))


b<-plot_df%>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point( aes(color = sig),
              alpha = 0.7) +
  geom_text_repel(data=plot_df%>% 
                    filter(Assay %in% oid_ecoli_binn_box),
                  aes(label = Assay),
                  size=3,
                  box.padding = 0.8,
                  segment.size = 0.5,
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  arrow = arrow(length = unit(0.010, "npc")),
                  nudge_x = .15,
                  nudge_y = .5)+
  labs(caption  = paste("Proteins included:", proteins, '\nSamples:', samples))+
  plot_theme

(b+a)+plot_annotation( tag_levels = 'A')

pnc_de<- de_func(df_pnc, column_name = 'Severe.pneumonia')%>%
  mutate(hda = 'pnc')

ecoli_de<- de_func(e_coli_df, column_name = 'Sepsis')%>%
  mutate(hda = 'ecoli')

volc_fs_ecoli<-plot_volc(e_coli_df, method = 'featseek', column_name = 'Sepsis')
box_fs_ecoli<-plot_box(e_coli_df, method = 'featseek', column_name = 'Sepsis')

volc_fs_ecoli<-plot_volc(df_pnc, method = 'featseek', column_name = 'Severe.pneumonia')
box_fs_ecoli<-plot_box(df_pnc, method = 'featseek', column_name = 'Severe.pneumonia')


##### Volcplot of BINN prot

# Plot for E coli
plot_df<-de_func(e_coli_df, column_name = 'Sepsis')

proteins<-nrow(plot_df)
samples<-nrow(method_func(e_coli_df))


plot_df%>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point( aes(color = sig),
              alpha = 0.7,
              show.legend = F) +
  geom_text_repel(data=plot_df%>% filter(Assay %in% c('MED18','COL6A3','SDC4',
                                                    'HMOX1',
                                                    'CSF3',
                                                    'TBL1X',
                                                    'TNFRSF13B',
                                                    'SOST',
                                                    'PDGFRB',
                                                    'FGF19')),
                  aes(label = Assay),
                  size=4,
                  min.segment.length =0.1,
                  max.overlaps = 7)+
  labs( title = 'E coli BINN proteins',
    caption  = paste("Proteins included:", proteins, '\nSamples:', samples))+
  plot_theme


#plot for pnc

plot_df<-de_func(df_pnc, column_name = 'Severe.pneumonia')

proteins<-nrow(plot_df)
samples<-nrow(method_func(df_pnc))


plot_df%>%
  ggplot(aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point( aes(color = sig),
              alpha = 0.7,
              show.legend = F) +
  geom_text_repel(data=plot_df%>% filter(Assay %in% c('SDC4',
                                                      'CSF3',
                                                      'SOST',
                                                      'FGF19',
                                                      'HMOX1',
                                                      'PDGFRB',
                                                      'CXCL16',
                                                      'DLL1',
                                                      'NOTCH1',
                                                      'COL6A3')),
                  aes(label = Assay),
                  size=4,
                  min.segment.length =0.1,
                  max.overlaps = 7)+
  labs( title = 'Pneumococcal Pneumonia BINN proteins',
    caption  = paste("Proteins included:", proteins, '\nSamples:', samples))+
  plot_theme
