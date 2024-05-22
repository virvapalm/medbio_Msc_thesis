
library(DOSE)
library(clusterProfiler)
library(AnnotationDbi)
library(ReactomePA)


reactome_db <- readRDS("Export/reactome.rds")
olink_meta <- read_csv("Export/Olink_meta.csv")


#Settings
minGSSize <- 10
maxGSSize <- Inf
pvalueCutoff <- 1
qvalueCutoff <- 1

universe <-
  olink_meta %>% 
  pull(UniProt) %>% 
  unique()

de_ecoli <- de_func(e_coli_df, column_name = 'Sepsis')

de_pnc <- de_func(df_pnc, column_name = 'Severe.pneumonia')



# GSEA GO
original_gene_list <- de_pnc$logFC #de_results contains results from differential expression analysis

# We name the vector with the UniProt IDs
names(original_gene_list) <- de_pnc$UniProt

# We sort the named logFC vector based on FC
gene_list = sort(original_gene_list, decreasing = TRUE)

gse <- gseGO(geneList = gene_list, 
             ont = "BP", 
             keyType = "UNIPROT", 
             #  nPerm = 10000, 
             #  minGSSize = 3, 
             #  maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "BH")

gse%>%
  as_tibble()%>%
  filter(p.adjust < 0.05)%>%
  view()

# KEGG ORA

ecoli <- de_ecoli |> 
  filter(adj.P.Val < 0.05) |> 
  pull(UniProt)

pnc <- 
  de_pnc |> 
  filter(adj.P.Val < 0.05) |> 
  pull(UniProt)

kegg_res <- 
  enrichKEGG(gene = pnc, 
             universe = universe, 
             keyType = "uniprot", 
             organism = "hsa", 
             pAdjustMethod = "BH")

dotplot(kegg_res, showCategory=10) 


kegg_res@result |> 
  as_tibble() |> 
  dplyr::select(category, subcategory, Description, p.adjust)%>%
  view()


dotplot(kegg_res, showCategory=10, split=".sign") + facet_grid(.~.sign)



# Reactome ORA
database <- 
  reactome_db %>% 
  filter(Species == "Homo sapiens") %>% 
  dplyr::select(term = Event_Name, gene = UniProt) %>% 
  distinct() |> 
  filter(gene %in% universe)

partition <- 
  de_ecoli |> 
  filter(adj.P.Val < 0.05) |> 
  #left_join(olink_meta, by = "Assay") |> 
  distinct(gene = UniProt)

res <- 
  enricher(partition$gene, 
           universe = universe,
           TERM2GENE = database, 
           minGSSize = minGSSize,
           maxGSSize = maxGSSize,
           pvalueCutoff = pvalueCutoff,
           qvalueCutoff = qvalueCutoff) |> 
  as_tibble()

res |> 
  arrange(-p.adjust)%>%
  dplyr::select(ID, Description, p.adjust)%>%
  view()


### Reactome enrichment analysis
universe_RA<-bitr(universe, fromType = "UNIPROT", 
                  toType = "ENTREZID", 
                  OrgDb = 'org.Hs.eg.db')

partition_RA <- 
  de_pnc |> 
  filter(adj.P.Val < 0.05) |> 
  #left_join(olink_meta, by = "Assay") |> 
  distinct(gene = UniProt)%>%
  pull(gene)

partion_RA<-bitr(partition_RA, fromType = "UNIPROT", 
                 toType = "ENTREZID", 
                 OrgDb = 'org.Hs.eg.db')

res <- 
  enrichPathway(partion_RA$ENTREZID, 
           universe = universe_RA$ENTREZID,
           minGSSize = minGSSize,
           maxGSSize = maxGSSize,
           pvalueCutoff = pvalueCutoff,
           qvalueCutoff = qvalueCutoff) |> 
  as_tibble()

res |> 
  arrange(-p.adjust)%>%
  head(n=10)

### Reactome gse

# GSEA GO
names<-bitr(de_ecoli$UniProt, fromType = "UNIPROT", 
            toType = "ENTREZID", 
            OrgDb = org.Hs.eg.db,
            drop = TRUE)

#de_results contains results from differential expression analysis
original_gene_list <- de_ecoli$logFC

# We name the vector with the UniProt IDs

names(original_gene_list)<-names$ENTREZID

# We sort the named logFC vector based on FC
gene_list = sort(original_gene_list, decreasing = TRUE)

pathway <- gsePathway(geneList = gene_list,
                      keyType = "UNIPROT", 
                      #  nPerm = 10000, 
                      #  minGSSize = 3, 
                      #  maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = TRUE, 
                      OrgDb = "org.Hs.eg.db", 
                      pAdjustMethod = "BH")%>%
  as_tibble()

ggplot(plot, aes(x = reorder(Description, -Count), y = Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 10 Pathways", x = "Pathway", y = "Gene Count") +
  theme_minimal()


reactome_db %>% 
  filter(Species == "Homo sapiens") %>%
  head()
s