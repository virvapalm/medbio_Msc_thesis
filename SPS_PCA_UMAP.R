### Setup ----
library(tidymodels)
library(docstring)
library(umap)
library(embed)
library(patchwork)

#### Set directory ----
## directory 
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

## Main directory 
main_dir <- getwd()

source('SPS_prep.R', local = TRUE)

### Download data ---- 
# KF metadata 290124
KS_meta <- read.csv("Export/KS_matched_meta.csv") 

#  metadata
df_DA_meta <- read.csv("Data/Export_meta_v2.0.2.csv")

clin_meta<- read.csv('Data/final_metadata_all_batches.csv')

Olink_meta <-read.csv('Export/Olink_meta.csv')

df_meta_pnc<-read.csv('Export/DE_meta_pnc.csv')



# Get NPX Data
df_NPX_unfiltered<-read.csv('Data/Export_expr_v2.0.1.csv')

#replicated protein are "IL6"   "TNF"   "CXCL8"
replicated_olink<-c("OID20101","OID20074","OID20153","OID21237","OID21276","OID21430","OID20473","OID20563","OID20631","OID20848","OID20911","OID20997")

# BINN and featseek OID

oid_feat_seek<-read.csv('Export/oid_feat_seek')%>%
  pull(selected)

pnc_uniprot <- read.csv('~/Documents/github/BINN_elin/data/pnc_biomarkers_new.csv', sep = ',')%>%
  colnames()

oid_pnc_binn<-Olink_meta%>%
  filter(Assay %in% pnc_uniprot)%>%
  pull(OlinkID)

ecoli_uniprot <- read.csv('~/Documents/github/BINN_elin/data/ecoli_biomarkers_new.csv', sep = ',')%>%
  colnames()

oid_ecoli_binn<-Olink_meta%>%
  filter(Assay %in% ecoli_uniprot)%>%
  pull(OlinkID)

sps_uniprot <- read.csv('~/Documents/github/BINN_elin/data/sps_biomarkers.csv', sep = ',')%>%
  colnames()

oid_sps_binn<-Olink_meta%>%
  filter(UniProt %in% sps_uniprot)%>%
  pull(OlinkID)


### Data wrangling ----

### NPX clean up

df_NPX<-
  df_DA_meta%>% # Takes the DA meta file
  filter(Cohort=='SPS')%>% #Filter the file for Sepsis
  left_join(df_NPX_unfiltered, by = 'DAid')%>% #Joins the meta file with the expression file by DAid
  select(-c(SampleID:Origin, Cohort, Class), -all_of(replicated_olink)) #Removes some non useful columns

sps_meta <- KS_meta%>%
  filter(Till.djup.analys == 1) %>%
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
         BMI = BMI) %>%
  mutate(HDA_diagnosis = factor(HDA_diagnosis, levels = 1:5, labels = c("Pnc", "Mycop", "E.coli", "S.aur.bacter", "Strep")), 
         Sepsis = factor(Sepsis, levels = c(0,1), labels = c(FALSE, TRUE)), 
         Alive90 = factor(Alive90, levels = c(0,1), labels = c(FALSE, TRUE)),
         Sex = factor(Sex, levels =c(1,2), labels = c('Female', 'Male')),
         Age = round(Age),
         Bacterimia_pnc = na_if(Bacterimia_pnc, 0),
         Bacterimia_ecoli = na_if(Bacterimia_ecoli, 0),
         Co_infection = ifelse(Co_infection == 1, TRUE, FALSE))%>%
  unite(Bacterimia, Bacterimia_pnc:Bacterimia_ecoli, na.rm=TRUE)%>%
  mutate(Bacterimia = factor(Bacterimia, levels = '1', labels = TRUE))%>%
  relocate(Alive90, .after = last_col())


# df_NPX_sps<-
#   df_DA_meta%>%
#   filter(Cohort=='SPS')%>%
#   left_join(df_NPX_unfiltered, by = 'DAid')%>%
#   #select(DAid, all_of(oid_feat_seek))%>%
#   select(DAid, starts_with('OID'))%>%
#   left_join(clin_meta, by = 'DAid')%>%
#   select(DAid, Cohort, Class, Disease, starts_with('OID'))%>%
#   right_join(sps_meta, by = 'DAid')%>%
#   relocate(HDA_diagnosis:Alive90, .after = Disease)%>%
#   mutate(Cohort = replace_na(Cohort, 'SPS'),
#          Class = replace_na(Class, 'Infection'),
#          Disease = ifelse(is.na(Disease),
#                           case_when(HDA_diagnosis == 'Pnc' ~ 'Pneumococcal Pneumonia',
#                                     HDA_diagnosis == 'E.coli' ~ 'E.coli pyelonephritis',
#                                     TRUE ~ Disease), Disease),
#          Sepsis = paste(Sepsis, HDA_diagnosis, sep = '_'))%>%
#   select(DAid, HDA_diagnosis,Sepsis, Cohort, Class, Disease, Sex, Age, starts_with('OID'))
# 
# df_NPX_infl<-
#   df_DA_meta%>%
#   filter(Cohort=='INFL')%>%
#   left_join(df_NPX_unfiltered, by = 'DAid')%>%
#   #select(DAid, all_of(oid_feat_seek))%>%
#   select(DAid, starts_with('OID'))%>%
#   left_join(clin_meta, by = 'DAid')%>%
#   filter(Class=='Healthy')%>%
#   mutate(Sex = factor(Gender, levels =c(0,1), labels = c('Female', 'Male')))%>%
#   select(-c(Vial.barcode:Updated, Diagnose, Subcategory:DAid_ext, Gender, BMI))%>%
#   relocate(Sex,.after = Cohort)%>%
#   relocate(Age:Class, .after=DAid)%>%
#   mutate(HDA_diagnosis = 'Healthy',
#          Sepsis = 'Healthy')%>%
#   select(DAid, HDA_diagnosis, Sepsis, Cohort, Class, Disease, Sex, Age, starts_with('OID'))


# Refactor metadata


## Find DAid for NPX expression

#complete_df<-rbind(df_NPX_sps, df_NPX_infl)

full_df <- left_join(sps_meta, df_NPX, by ='DAid')


### E coli df --
e_coli_df<-
  sps_meta%>%
  filter(HDA_diagnosis == 'E.coli')%>%
  left_join(df_NPX, by = 'DAid')

ecoli_test<-
  e_coli_df%>%
  select(DAid:Alive90, all_of(ecoli_feat))

### Pncc df

df_pnc <- df_meta_pnc%>%
  relocate(Alive90, .after = last_col())%>%
  inner_join(., df_NPX, by = 'DAid')

df_pnc_test<-df_pnc%>%
  select(DAid:Alive90, all_of(pnc_feat))

### Creating func ---

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



pca_prep_fn<-function(df, thold = 1, method = NA, disease = NA){
  
  #' PCA prep function
  #' 
  #' This function has a large recipie output.
  #' It takes a dataframe and a potential threshold input.The threshold defult is 0.75. The threshold can be increased/decreased to incresease/decrease number of included variables.
  #' 
  #'@param df: The dataframe input. Could be 'e_coli_df', 'df_pnc' or 'full_df'
  #'@param thold: a integer input for threshold value. should be between 0 and 1. A value close to 1 include more proteins.  
  #'
  #'@examples: 
  #'pca_prep_fn(df_pnc)
  #'pca_prep_fn(e_coli_df, thold = 1)
  #'pca_prep_fn(full_df, thold = 0.2)
  #'
  
  df_funk<- method_func(df, method = method, disease = disease)
  
  set.seed(123)
  
  pca_prep<-
    recipe(~., data = df_funk) %>%
    update_role(DAid:Alive90,
                new_role = "id") %>%
    step_corr(all_numeric(), threshold = thold) %>%
    step_impute_knn(all_numeric())%>%
    step_pca(all_predictors())%>%
    prep()
    
}

pca_vis<- function(df, thold=1, method = NA, disease = NA, color_var){
  
  #' PCA visualization funciton
  #' 
  #' This function has a dotplot output.
  #' It takes a dataframe and a column name.
  #' It also takes a thold input. The defult thold is 0.75. 
  #' The threshold can be increased/decreased to incresease/decrease number of included variables.
  #' 
  #' The function takes the df and generates as large recipie using the pca_prep_fn. The recipie is then used for the vizualisation.
  #' For more info on the pca_prep_fn, see pca_prep_fn. 
  #' 
  #' 
  #'@param df: The dataframe input. Could be 'e_coli_df', 'df_pnc' or 'full_df'
  #'@param color_var: a varible to which dye the dots by. Should be a character string of dateframe column name.
  #'@param thold: a integer input for threshold value. should be between 0 and 1. A value close to 1 include more proteins.  
  #'
  #'@examples: 
  #'pca_vis(df_pnc, 'Sepsis')
  #'pca_vis(e_coli_df, thold = 1, 'Alive90')
  #'pca_vis(full_df, thold = 0.2, 'PI')
  #'
  
  pca_prep<-pca_prep_fn(df, method = method, disease = disease) 
  
  #removals <- length(pca_prep$steps[[1]]$removals)
  #keep <- df%>%select(!DAid:Alive90)%>%ncol() - removals
  
  proteins<-
    pca_prep$var_info%>%
    filter(role=='predictor')%>%
    nrow()
  
  samples<-pca_prep$template%>%nrow()
  
  pca_prep%>%
    bake(new_data = NULL) %>%
    as_tibble()%>%
    ggplot(aes(x = PC1, y = PC2, color = !!sym(color_var))) + # !!sym() for using a string as a variable name
    geom_point(alpha = 0.7)+
    labs(caption  = paste("Proteins included:", proteins, '\nSamples', samples))+
    plot_theme
}

scree_fn<-function(df, thold=1, method = NA, disease = NA){
  
  #' Scree funciton
  #' 
  #' This function has a dotplot output of the vriance explained by the 9 first components in PCA.
  #' It takes a dataframe and a potential threshold input.The threshold defult is 0.75. The threshold can be increased/decreased to incresease/decrease number of included variables.
  #' 
  #' The function takes the df and generates as large recipie using the pca_prep_fn. The recipie is then used for the vizualisation.
  #' For more info on the pca_prep_fn, see pca_prep_fn. 
  #' 
  #' 
  #'@param df: The dataframe input. Could be 'e_coli_df', 'df_pnc' or 'full_df'
  #'@param thold: a integer input for threshold value. should be between 0 and 1. A value close to 1 include more proteins.  
  #'
  #'@examples: 
  #'scree_fn(df_pnc)
  #'scree_fn(e_coli_df, thold = 1)
  #'scree_fn(full_df, thold = 0.2)
  #'
  
  
  
  pca_prep<-pca_prep_fn(df, thold) 
  
  pca_summary_stat_grep<-summary(pca_prep$steps[[3]]$res)
  
  pca_summary_stat_grep_importance<-pca_summary_stat_grep$importance
    
  pca_summary_stat<-
    pca_summary_stat_grep_importance%>%
    as.data.frame()%>%
    select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9)%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column(., var= 'PC')%>%
    rename(variance = `Proportion of Variance`,
           std = `Standard deviation`,
           culmative_proportion = `Cumulative Proportion`)
    
    pca_summary_stat %>%
      ggplot(aes(x= PC, y = variance))+
      geom_point( size= 2)+
      scale_y_continuous(expand = c(0,1)) +
      labs(y ="",
           title ="Proportion of variance")+
      plot_theme
    
}

umap_fn<-function(df, thold=1, method = NA, disease = NA, color_var, neighbours = 5){
  #' UMAP funciton
  #' 
  #' This function has a dotplot output.
  #' It takes a dataframe and a column name.
  #' It also takes a thold input. The defult thold is 0.75. 
  #' The threshold can be increased/decreased to incresease/decrease number of included variables.
  #' 
  #' The function takes the df and generates as large recipie. The recipie is then used for the visualization. 
  #' 
  #' 
  #'@param df: The dataframe input. Could be 'e_coli_df', 'df_pnc' or 'full_df'
  #'@param color_var: a varible to which dye the dots by. Should be a character string of dateframe column name.
  #'@param thold: a integer input for threshold value. should be between 0 and 1. A value close to 1 include more proteins.  
  #'
  #'@examples: 
  #'pca_vis(df_pnc, 'Sepsis')
  #'pca_vis(e_coli_df, thold = 1, 'Alive90')
  #'pca_vis(full_df, thold = 0.2, 'PI')
  #'
  
  
  data<-method_func(df, method = method, disease = disease)
  
  set.seed(123)
  
  umap_prep <-
    recipe(~., data = data) %>%
    update_role(DAid:Alive90,
                new_role = "id") %>%
    step_corr(all_numeric(), threshold = 1) %>%
    step_impute_knn(all_numeric())%>%
    step_umap(all_predictors(),
              neighbors = neighbours)%>%
    prep()
  
  #removals <- length(umap_prep$steps[[1]]$removals)
  #keep <- df%>%select(!DAid:Alive90)%>%ncol() - removals
  
  proteins<-
    umap_prep$var_info%>%
    filter(role=='predictor')%>%
    nrow()
  
  samples<-umap_prep$template%>%nrow()
  
  umap_prep%>%
    bake(new_data = NULL) %>%
    as_tibble() %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = !!sym(color_var))) + # !!sym() for using a string as a variable name
    geom_point(alpha =0.7)+
    labs(caption  = paste("Proteins included:", proteins, '\nSamples:', samples))+
    plot_theme
  
}


### PCA following Feature selection --- 

e_coli_pca_all<-pca_vis(e_coli_df, color_var = 'Sepsis') 
pnc_pca_all<-pca_vis(df_pnc, color_var = 'Severe.pneumonia')

(e_coli_pca_all+pnc_pca_all)/(e_coli_umap_all+ pnc_umap_all)+
  plot_annotation( tag_levels = 'A')

e_coli_pca_featseek<-pca_vis(e_coli_df,method = 'featseek', color_var = 'Sepsis')+ggtitle('E.Coli')
pnc_pca_featseek<-pca_vis(df_pnc, method = 'featseek', color_var = 'Severe.pneumonia') + ggtitle('Pnc')

e_coli_pca_featseek + pnc_pca_featseek +plot_annotation( tag_levels = 'A')

e_coli_pca_binn<-pca_vis(e_coli_df,method = 'BINN', color_var = 'Sepsis')+ggtitle('E.Coli')
pnc_pca_binn<-pca_vis(df_pnc, method = 'BINN', color_var = 'Severe.pneumonia')+ + ggtitle('Pnc')

e_coli_umap_all<-umap_fn(e_coli_df, color_var = 'Sepsis') 
pnc_umap_all<-umap_fn(df_pnc, color_var = 'Severe.pneumonia')

e_coli_umap_featseek<-umap_fn(e_coli_df,method = 'featseek', color_var = 'Sepsis')
pnc_umap_featseek<-umap_fn(df_pnc, method = 'featseek', color_var = 'Severe.pneumonia')

e_coli_umap_binn<-umap_fn(e_coli_df,method = 'BINN', color_var = 'Sepsis')
pnc_umap_binn<-umap_fn(df_pnc, method = 'BINN', color_var = 'Severe.pneumonia')

(e_coli_umap_featseek+pnc_umap_featseek)+plot_annotation( tag_levels = 'A')

pnc_pca_diff<- pca_vis(df_pnc_test, color_var = 'Severe.pneumonia')
umap_pnc_diff<-umap_fn(df_pnc_test, color_var = 'Severe.pneumonia')

ecoli_pca_diff<- pca_vis(ecoli_test, color_var = 'Sepsis')
ecoli_umap_diff<-umap_fn(ecoli_test, color_var = 'Sepsis')

(e_coli_pca_all+ecoli_pca_diff)+ plot_annotation( tag_levels = 'A')
(pnc_pca_all+pnc_pca_diff) + plot_annotation( tag_levels = 'A')

(e_coli_umap_all + ecoli_umap_diff)/(pnc_umap_all + umap_pnc_diff) + plot_annotation( tag_levels = 'A')

(e_coli_umap_binn+pnc_umap_binn)+ plot_annotation( tag_levels = 'A')

umap_ecoli_fs<-umap_fn(ecoli_test, color_var = 'Sepsis')
umap_pnc_fs<-umap_fn(df_pnc_test, color_var = 'Severe.pneumonia' )

umap_ecoli_binn<- umap_fn(e_coli_df, method = 'BINN', color_var = 'Sepsis')



volc_fs_ecoli+(umap_ecoli_fs/box_fs_ecoli)+
  plot_annotation(tag_levels = 'A')

volc_fs_ecoli + (umap_pnc_fs/box_fs_ecoli)+
  plot_annotation(tag_levels = 'A')
