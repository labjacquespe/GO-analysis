# Install necessary packages if not already installed
if (!requireNamespace("clusterProfiler", quietly = TRUE)) install.packages("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

# Load the required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(grid)
library(ggrepel)

#set wd
setwd("/Volumes/homedir$/Ongoing projects/MF GDM Gen3G validation study/gene expression data from Frede/R code")

# Step 1: Load the CSV files for males and females
male_data <- read_csv("male_strata.csv")
female_data <- read_csv("female_strata.csv")

# Step 2: Filter genes with logFC > 0.2 or logFC < -0.2, and P.Value < 0.05
filtered_genes_male <- male_data %>%
  filter((logFC > 0.2 | logFC < -0.2) & P.Value < 0.05)

filtered_genes_female <- female_data %>%
  filter((logFC > 0.2 | logFC < -0.2) & P.Value < 0.05)

# Step 3: Remove extra numbers from Ensembl IDs (removes versioning)
filtered_genes_male$gene_id <- gsub("\\..*", "", filtered_genes_male$Name)
filtered_genes_female$gene_id <- gsub("\\..*", "", filtered_genes_female$Name)

# Step 4: Map Ensembl IDs to Entrez IDs for male and female datasets
mapped_genes_male <- bitr(filtered_genes_male$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
mapped_genes_female <- bitr(filtered_genes_female$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Step 5: Perform GO enrichment analysis
go_enrich_male <- enrichGO(
  gene          = mapped_genes_male$ENTREZID,  # Use Entrez IDs
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",  # "BP", "MF", or "CC" for individual ontology
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

go_enrich_female <- enrichGO(
  gene          = mapped_genes_female$ENTREZID,  # Use Entrez IDs
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",  # "BP", "MF", or "CC" for individual ontology
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)


# Step 6: Extract significant GO terms and associated genes

# For males
significant_male_go <- go_enrich_male@result %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(ID, Description, geneID, p.adjust)

# For females
significant_female_go <- go_enrich_female@result %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(ID, Description, geneID, p.adjust)

# Step 7: Expand the 'geneID' column into individual genes (separated by "/")
significant_male_go <- significant_male_go %>%
  mutate(genes = strsplit(geneID, "/")) %>%
  unnest(genes)

significant_female_go <- significant_female_go %>%
  mutate(genes = strsplit(geneID, "/")) %>%
  unnest(genes)

# Step 8: Convert Entrez IDs to Gene Symbols (gene names)

# For males: Map Entrez IDs to gene names
gene_names_male <- bitr(significant_male_go$genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# For females: Map Entrez IDs to gene names
gene_names_female <- bitr(significant_female_go$genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Merge the gene symbols back with the GO results
significant_male_go <- significant_male_go %>%
  left_join(gene_names_male, by = c("genes" = "ENTREZID")) %>%
  select(ID, Description, SYMBOL, p.adjust) %>%
  rename(Gene_Name = SYMBOL)

significant_female_go <- significant_female_go %>%
  left_join(gene_names_female, by = c("genes" = "ENTREZID")) %>%
  select(ID, Description, SYMBOL, p.adjust) %>%
  rename(Gene_Name = SYMBOL)

# Step 9: Combine the gene names for each GO term into a single field (genes separated by "/")

# For males: Group by GO term and collapse gene names
significant_male_go <- significant_male_go %>%
  group_by(ID, Description, p.adjust) %>%
  summarise(Gene_Names = paste(Gene_Name, collapse = "/")) %>%
  ungroup()

# For females: Group by GO term and collapse gene names
significant_female_go <- significant_female_go %>%
  group_by(ID, Description, p.adjust) %>%
  summarise(Gene_Names = paste(Gene_Name, collapse = "/")) %>%
  ungroup()

# Step 10: Save the results as CSV files

# For males
write_csv(significant_male_go, "significant_go_terms_male_logFC2_p05.csv")

# For females
write_csv(significant_female_go, "significant_go_terms_female_logFC2_p05.csv")








#### merge into one dotplot
# Step 1: Extract results and add a new column to indicate the group (male/female)
male_results <- go_enrich_male@result %>%
  mutate(Group = "Male", Ontology = ONTOLOGY)  # Add "Male" group label and Ontology

female_results <- go_enrich_female@result %>%
  mutate(Group = "Female", Ontology = ONTOLOGY)  # Add "Female" group label and Ontology

# Step 2: Combine the two results based on the common GO terms (you can merge if you want only shared terms)
combined_results <- bind_rows(male_results, female_results)

# Step 3: Filter for significant GO terms (optional, e.g., adjust p-value < 0.05)
combined_results <- combined_results %>%
  filter(p.adjust < 0.05)  # Filter by adjusted p-value for significance

# Step 4: Rename the ONTOLOGY values to full names
combined_results <- combined_results %>%
  mutate(Ontology = case_when(
    ONTOLOGY == "BP" ~ "Biological Process",
    ONTOLOGY == "CC" ~ "Cellular Component",
    ONTOLOGY == "MF" ~ "Molecular Function",
    TRUE ~ ONTOLOGY  # Leave as is for any other cases (if present)
  ))

# Step 4: Remove rows with zero or missing `Count` (these rows will not have dots) per Ontology facet
combined_results <- combined_results %>%
  group_by(Ontology) %>%
  filter(!is.na(Count) & Count > 0) %>%
  ungroup()

combined_results <- combined_results %>%
  mutate(Ontology = case_when(
    ONTOLOGY == "BP" ~ "Biological Process",
    ONTOLOGY == "MF" ~ "Molecular Function",
    TRUE ~ NA_character_  # Remove Cellular Component (CC) or other cases
  )) %>%
  filter(!is.na(Ontology))  # Only keep BP and MF


# Step 5: Create a dotplot with vertical faceting by Ontology (BP, MF, CC) and adjust facet placement - Figure 2C
dotplot_combined <- ggplot(combined_results, aes(x = Group, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +  # Color by p-value
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),  # Adjust text size for y-axis
    plot.margin = unit(c(1.2, 1, 1, 2), "cm"),  # Increase margin for y-axis labels
    strip.placement = "outside",  # Place facet labels outside
    strip.text.y.right = element_text(angle = 0),  # Keep facet labels horizontal on the right
    panel.spacing = unit(0.01, "lines"),  # Reduce spacing between panels
  ) +
  labs(
    x = "Fetal sex", 
    y = "GO Terms",
    color = "Adjusted p-value", 
    size = "Count"
  ) +
  facet_wrap(Ontology ~ .,ncol=1,scale="free_y",switch="y",drop=FALSE)  # Use full names for facet

# Display the plot
print(dotplot_combined)


##### plot DEG in heatmap - Figure 2B ####### 

# Step 1: Find the common genes that are differentially expressed in both males and females
common_genes <- intersect(filtered_genes_male$Description, filtered_genes_female$Description)

# Step 2: Create a dataset for the heatmap using only the common differentially expressed genes
# Merge male and female datasets based on the "Description" column
male_common <- male_data %>%
  filter(Description %in% common_genes) %>%
  select(Description, logFC) %>%
  rename(male = logFC)

female_common <- female_data %>%
  filter(Description %in% common_genes) %>%
  select(Description, logFC) %>%
  rename(female = logFC)

# Merge the male and female data
heatmap_data <- merge(male_common, female_common, by = "Description", all = TRUE)

# Step 5: Add a column to flag genes with the same directionality (logFC_male * logFC_female > 0)
heatmap_data <- heatmap_data %>%
  mutate(directionality = ifelse(male * female > 0, 1, 0))

# Step 6: Sort the data:
# - First by directionality (1: same, 0: opposite)
# - Then by logFC_male for genes with opposite directionality
heatmap_data <- heatmap_data %>%
  arrange(desc(directionality), 
          ifelse(directionality == 1, male, NA),
          ifelse(directionality == 0, male, NA))

# Step 7: Convert the data frame to a matrix for heatmap (remove the "Description" and "directionality" columns)
heatmap_matrix <- as.matrix(heatmap_data[, c("male", "female")])

# Assign gene names as row names
rownames(heatmap_matrix) <- heatmap_data$Description

# Step 8: Define the color scale with white representing logFC = 0
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

# Set the color breaks to center white at 0
breaks <- seq(min(-1.5), 
              max(1.5, na.rm = TRUE), 
              length.out = 51)  # 50 breaks ensure white is centered at 0

# Step 9: Create the heatmap using pheatmap with rescaled colors
p<-pheatmap(heatmap_matrix,
         color = color_palette,  # Custom color palette
         breaks = breaks,  # Custom breaks to center white at 0
         cluster_rows = FALSE,  # Disable clustering for custom sorting
         cluster_cols = FALSE,  # Do not cluster columns (males vs females)
         show_rownames = TRUE,  # Show gene names
         show_colnames = TRUE,  # Show male/female labels
         scale = "none",  # Do not scale the logFC values
         fontsize_row = 12,  # Adjust font size for gene names
         fontsize_col = 12, # Adjust font size for column names
         legend = TRUE)  


