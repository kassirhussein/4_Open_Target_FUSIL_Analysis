# For all the Open Target Data set 

# To see Gene targets Paralogue Presence for gene target between withdrawn and approved

library(biomaRt)  # Load biomaRt for querying Ensembl gene annotations

# Extract unique human gene symbols from OpenTargets dataset
human_gene_list <- opentarget_data %>%
  dplyr::select(approvedSymbol) %>%
  unique()

# Flatten the data frame to a simple character vector of gene names
human_gene_list2 <- human_gene_list$approvedSymbol

# Connect to the Ensembl Human dataset
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query paralog information for each gene using biomaRt
human_gene_paralogues <- getBM(attributes = c("ensembl_gene_id",                     # Ensembl gene ID
                                              "external_gene_name",                 # HGNC gene symbol
                                              "hsapiens_paralog_ensembl_gene",      # Paralog Ensembl ID
                                              "hsapiens_paralog_associated_gene_name", # Paralog gene name
                                              "hsapiens_paralog_perc_id",           # % identity gene → paralog
                                              "hsapiens_paralog_perc_id_r1"),       # % identity paralog → gene
                               filters = "external_gene_name",                      # Filter by gene symbol
                               values = human_gene_list2,
                               mart = human)

# Rename for readability: keeping only the symbol for consistency
human_gene_paralogues <- human_gene_paralogues %>%
  rename(gene_symbol = external_gene_name )

# Identify genes that have no paralogues or very low similarity (< 30%)
na_paralogue_count <- human_gene_paralogues %>%
  group_by(gene_symbol) %>%
  filter(all(is.na(hsapiens_paralog_perc_id) | hsapiens_paralog_perc_id < 30 ))

# Identify genes with meaningful paralogues (similarity > 30%)
Genes_with_paralogue <- human_gene_paralogues %>%
  group_by(gene_symbol) %>%
  na.omit() %>%
  filter(hsapiens_paralog_perc_id > 30)


# --------------------- Withdrawn Genes Analysis ---------------------

# Select unique withdrawn gene targets
withdrawn_genes <- withdrawn_drugs %>%
  dplyr::select(approvedSymbol) %>%
  unique()

# Create reference lists for classification
with_paralogue <- unique(Genes_with_paralogue$gene_symbol)
without_paralogue <- unique(na_paralogue_count$gene_symbol)

# Classify each withdrawn gene by paralogue status: with, without, both, or not found
withdrawn_genes$paralogue_status <- ifelse(
  withdrawn_genes$approvedSymbol %in% with_paralogue & withdrawn_genes$approvedSymbol %in% without_paralogue, 
  "both",
  ifelse(
    withdrawn_genes$approvedSymbol %in% with_paralogue, 
    "with_paralogues",
    ifelse(
      withdrawn_genes$approvedSymbol %in% without_paralogue,
      "without_paralogues",
      "not_found"
    )
  )
)

# Count how many genes fall into each paralogue status
status_counts <- table(withdrawn_genes$paralogue_status)

# Convert counts to percentages
status_percentages <- prop.table(status_counts) * 100

# Convert to data frame for plotting
status_df <- as.data.frame(status_percentages)
colnames(status_df) <- c("paralogue_status", "percentage")

# Visualize withdrawn genes by paralogue status
ggplot(status_df, aes(x = paralogue_status, y = percentage, fill = paralogue_status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5) +
  labs(title = "Percentage of Genes by Paralogous Status ( Withdrawn )",
       x = "Paralogue Status",
       y = "Percentage") +
  theme_minimal()


# --------------------- Approved Genes Analysis ---------------------

# Select unique approved gene targets
approved_genes <- opentarget_approved %>%
  dplyr::select(approvedSymbol) %>%
  unique()

# Reuse lists of genes with and without paralogues
with_paralogue <- unique(Genes_with_paralogue$gene_symbol)
without_paralogue <- unique(na_paralogue_count$gene_symbol)

# Classify each approved gene by paralogue status
approved_genes$paralogue_status <- ifelse(
  approved_genes$approvedSymbol %in% with_paralogue & approved_genes$approvedSymbol %in% without_paralogue, 
  "both",
  ifelse(
    approved_genes$approvedSymbol %in% with_paralogue, 
    "with_paralogues",
    ifelse(
      approved_genes$approvedSymbol %in% without_paralogue,
      "without_paralogues",
      "not_found"
    )
  )
)

# Count approved genes per paralogue status
status_counts <- table(approved_genes$paralogue_status)

# Convert counts to percentages
status_percentages <- prop.table(status_counts) * 100

# Prepare data for visualization
status_df <- as.data.frame(status_percentages)
colnames(status_df) <- c("paralogue_status", "percentage")

# Plot approved genes by paralogue status
ggplot(status_df, aes(x = paralogue_status, y = percentage, fill = paralogue_status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5) +
  labs(title = "Percentage of Genes by Paralogous Status ( Approved )",
       x = "Paralogue Status",
       y = "Percentage") +
  theme_minimal()



#############################################################
# --------------------- With FUSIL DATA ---------------------
#############################################################

# For all the Open Target Data set 

# To see Gene targets Paralogue Presence for gene target between withdrawn and approved

library(biomaRt)  # Load biomaRt for querying Ensembl gene annotations

# Extract unique human gene symbols from OpenTargets dataset
human_gene_list <- opentarget_data %>%
  dplyr::select(approvedSymbol) %>%
  unique()

# Flatten the data frame to a simple character vector of gene names
human_gene_list2 <- human_gene_list$approvedSymbol

# Connect to the Ensembl Human dataset
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query paralog information for each gene using biomaRt
human_gene_paralogues <- getBM(attributes = c("ensembl_gene_id",                     # Ensembl gene ID
                                              "external_gene_name",                 # HGNC gene symbol
                                              "hsapiens_paralog_ensembl_gene",      # Paralog Ensembl ID
                                              "hsapiens_paralog_associated_gene_name", # Paralog gene name
                                              "hsapiens_paralog_perc_id",           # % identity gene → paralog
                                              "hsapiens_paralog_perc_id_r1"),       # % identity paralog → gene
                               filters = "external_gene_name",                      # Filter by gene symbol
                               values = human_gene_list2,
                               mart = human)

# Rename for readability: keeping only the symbol for consistency
human_gene_paralogues <- human_gene_paralogues %>%
  rename(gene_symbol = external_gene_name )

# Identify genes that have no paralogues or very low similarity (< 30%)
na_paralogue_count <- human_gene_paralogues %>%
  group_by(gene_symbol) %>%
  filter(all(is.na(hsapiens_paralog_perc_id) | hsapiens_paralog_perc_id < 30 ))

# Identify genes with meaningful paralogues (similarity > 30%)
Genes_with_paralogue <- human_gene_paralogues %>%
  group_by(gene_symbol) %>%
  na.omit() %>%
  filter(hsapiens_paralog_perc_id > 30)



# Integrate FUSIL DATA

opentarget_data_fusil <- opentarget_data %>%
  left_join(fusil, by = c ( "approvedSymbol" = "gene_symbol")) %>%
  na.omit()



# --------------------- Withdrawn Genes Analysis ---------------------

# Select unique withdrawn gene targets
withdrawn_genes <- opentarget_data_fusil %>%
  filter(status == "Withdrawn") %>%
  dplyr::select(approvedSymbol, fusil) %>%
  unique()

# Create reference lists for classification
with_paralogue <- unique(Genes_with_paralogue$gene_symbol)
without_paralogue <- unique(na_paralogue_count$gene_symbol)

# Classify each withdrawn gene by paralogue status: with, without, both, or not found
withdrawn_genes$paralogue_status <- ifelse(
  withdrawn_genes$approvedSymbol %in% with_paralogue & withdrawn_genes$approvedSymbol %in% without_paralogue, 
  "both",
  ifelse(
    withdrawn_genes$approvedSymbol %in% with_paralogue, 
    "with_paralogues",
    ifelse(
      withdrawn_genes$approvedSymbol %in% without_paralogue,
      "without_paralogues",
      "not_found"
    )
  )
)


withdrawn_genes2 <- withdrawn_genes %>%
  count(fusil, paralogue_status) %>%
  group_by(fusil) %>%
  mutate(percentage = (n/sum(n))*100)


# Visualize withdrawn genes by paralogue status
ggplot(withdrawn_genes2, aes(x = fusil, y = percentage, fill = paralogue_status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5) +
  labs(title = "Percentage of Genes by Paralogous Status ( Withdrawn )",
       x = "Paralogue Status",
       y = "Percentage") +
  theme_minimal()



# --------------------- Approved Genes Analysis ---------------------

# Select unique withdrawn gene targets
approved_genes <- opentarget_data_fusil %>%
  filter(phase == "4" | status == "Completed") %>%
  dplyr::select(approvedSymbol, fusil) %>%
  unique()

# Create reference lists for classification
with_paralogue <- unique(Genes_with_paralogue$gene_symbol)
without_paralogue <- unique(na_paralogue_count$gene_symbol)

# Classify each withdrawn gene by paralogue status: with, without, both, or not found
approved_genes$paralogue_status <- ifelse(
  approved_genes$approvedSymbol %in% with_paralogue & approved_genes$approvedSymbol %in% without_paralogue, 
  "both",
  ifelse(
    approved_genes$approvedSymbol %in% with_paralogue, 
    "with_paralogues",
    ifelse(
      approved_genes$approvedSymbol %in% without_paralogue,
      "without_paralogues",
      "not_found"
    )
  )
)


approved_genes2 <- approved_genes %>%
  count(fusil, paralogue_status) %>%
  group_by(fusil) %>%
  mutate(percentage = (n/sum(n))*100)


# Visualize withdrawn genes by paralogue status
ggplot(approved_genes2, aes(x = fusil, y = percentage, fill = paralogue_status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5) +
  labs(title = "Percentage of Genes by Paralogous Status ( Approved )",
       x = "Paralogue Status",
       y = "Percentage") +
  theme_minimal()
















