library(eulerr)
library(ggplotify)
de_results <- readRDS("diff_df.rds")

pnc <- 
  de_results |> 
  filter(adj_p_val_pnc < 0.05) |> 
  pull(Assay)

ecoli <- 
  de_results |> 
  filter(adj_p_val_ecoli < 0.05) |> 
  pull(Assay)


y <- list("Pnc" = pnc, 
          "Ecoli" = ecoli)
plot(euler(y, shape = "ellipse"), quantities = TRUE, fills = c("#F2CDC4","#C2E6DF")) |> as.ggplot()

ggsave(savepath("overlap_ecoli_pnc.png"), h = 3, w = 3)



### Diff prot


ecoli_all<-de_func(e_coli_df, column_name = 'Sepsis')%>%
  filter(adj.P.Val<0.05)%>%
  pull(OlinkID)

ecoli_fs<-de_func(e_coli_df, method = 'featseek', column_name = 'Sepsis')%>%
  filter(adj.P.Val<0.05)%>%
  filter(OlinkID %in% oid_feat_seek)%>%
  pull(Assay)

ecoli_binn<-Olink_meta%>%
  filter(Assay %in% c('MED18','COL6A3','SDC4',
                      'HMOX1',
                      'CSF3',
                      'TBL1X',
                      'TNFRSF13B',
                      'SOST',
                      'PDGFRB',
                      'FGF19'))%>%
  pull(OlinkID)

ecoli<-list(
  'All' = ecoli_all,
  'FS' = oid_feat_seek,
  'BINN' = ecoli_binn
)

plot(euler(ecoli_euler, shape = "circle"), quantities = TRUE) |> as.ggplot()+
  labs( title = 'Protein overlap between cohort and method')+
  theme(plot.title = element_text(hjust = 0.5))
 

ecoli_euler<-list(
  'DE Ecoli' = ecoli_all,
  'DE Pnc'= pnc_all,
  'uFS' = oid_feat_seek
)

pnc_all<-de_func(df_pnc, column_name = 'Severe.pneumonia')%>%
  filter(adj.P.Val<0.05)%>%
  pull(OlinkID)

pnc_euler<-list(
  'DE' = pnc_all,
  'uFS' = oid_feat_seek,
  'BINN' = oid_pnc_binn
)

b<-plot(euler(pnc_euler, shape = "circle"), quantities = TRUE) |> as.ggplot()+
  labs( title = 'Pneumococcal pneumonia cohort')+
  theme(plot.title = element_text(hjust = 0.5))

Olink_meta%>%
  filter(OlinkID %in% oid_feat_seek)%>%
  filter(OlinkID %in% pnc_all)%>%
  filter(!OlinkID %in% ecoli_all)%>%
  pull(Assay)

Olink_meta%>%
  filter(OlinkID %in% oid_feat_seek)%>%
  filter(OlinkID %in% pnc_all)%>%
  filter(!OlinkID %in% ecoli_all)%>%
  #filter(Assay %in% binn_pnc)
  pull(Assay)


a+b+plot_annotation(tag_levels = 'A')


binn_ecoli<-c('MED18','COL6A3','SDC4',
'HMOX1',
'CSF3',
'TBL1X',
'TNFRSF13B',
'SOST',
'PDGFRB',
'FGF19')

binn_pnc<-c('SDC4',
            'CSF3',
            'SOST',
            'FGF19',
            'HMOX1',
            'PDGFRB',
            'CXCL16',
            'DLL1',
            'NOTCH1',
            'COL6A3')

new_list<-list(
  'binn ecoli'= binn_ecoli,
  'binn pnc' = binn_pnc
)

plot(euler(new_list, shape = "circle"), quantities = TRUE) |> as.ggplot()+
  labs( title = 'Pneumococcal pneumonia cohort')+
  theme(plot.title = element_text(hjust = 0.5))

Olink_meta%>%
  filter(Assay %in% binn_ecoli)%>%
  filter(!Assay %in% binn_pnc)
