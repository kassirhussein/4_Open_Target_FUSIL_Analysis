# Title: "A Systems-Level View of Drug Targeting: From Molecular Action to Disease Classes"


# Load tidyverse for data manipulation and visualization
library(tidyverse)
library(scales)

# ------------------------ DATA IMPORT ------------------------

# Import human gene paralogues data
human_gene_paralogues <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/paralogues/human_gene_paralogues.csv", 
                                  col_types = cols(...1 = col_skip()))

# Import FUSIL gene essentiality data
fusil <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/fusil.csv")

# Import Open Targets drug data
opentarget_data <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/Open Target drug data/Opentarget_data.csv", 
                            col_types = cols(...1 = col_skip()))

# Fill missing values in 'status' column with 'Unknown status'
opentarget_data <- opentarget_data %>%
  replace_na(list(status = 'Unknown status'))

# Import DrugBank drug-target interaction data
DrugBank_data <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/Drug bank data/DrugBank_data.csv", 
                          col_types = cols(...1 = col_skip()))

DrugBank_data <- DrugBank_data[!is.na(DrugBank_data$gene_name), ]

# --- Import a List of Protein Coding GENES ---

pc_genes <- read.delim("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/locus_groups/protein-coding_gene.txt")


# ------------------------ CLUSTERING ACTION TYPES ------------------------

# Define a mapping from specific action types to broader clusters
action_cluster_map <- list(
  "Inhibitors/Blockers" = c("INHIBITOR", "BLOCKER", "ANTISENSE INHIBITOR", 
                            "CROSS-LINKING AGENT", "DISRUPTING AGENT"),
  "Agonists/Activators" = c("AGONIST", "ACTIVATOR", "PARTIAL AGONIST", 
                            "INVERSE AGONIST", "OPENER"),
  "Modulators" = c("MODULATOR", "POSITIVE MODULATOR", "POSITIVE ALLOSTERIC MODULATOR",
                   "NEGATIVE ALLOSTERIC MODULATOR", "STABILISER"),
  "Antagonists" = c("ANTAGONIST", "ALLOSTERIC ANTAGONIST"),
  "Agents" = c("BINDING AGENT", "RELEASING AGENT", "SUBSTRATE"),
  "Biological Constructs & Others" = c("EXOGENOUS GENE", "EXOGENOUS PROTEIN",
                                       "HYDROLYTIC ENZYME", "PROTEOLYTIC ENZYME",
                                       "RNAI INHIBITOR", "OTHER")
)

# Function to map each action type to its cluster
map_to_cluster <- function(action) {
  for (cluster in names(action_cluster_map)) {
    if (action %in% action_cluster_map[[cluster]]) {
      return(cluster)
    }
  }
  return("Unclassified")
}

# Apply the clustering function to the Open Targets dataset
opentarget_data <- opentarget_data %>%
  mutate(actionCluster = sapply(actionType, map_to_cluster))

# ------------------------ FILTERING APPROVED DRUGS ------------------------

# Filter DrugBank for non-withdrawn drugs (approved)
drugbank_approved <- DrugBank_data %>%
  filter(!grepl("withdrawn", groups, ignore.case = TRUE))
  

# Filter Open Targets for drugs with 'Completed' status (assumed approved)
opentarget_approved <- opentarget_data %>%
  filter(phase == "4" | status == "Completed")

# ------------------------ DATASET COMPARISON ------------------------

# Count number of unique drugs in Open Targets and DrugBank
count_drugs_opentarget <- n_distinct(opentarget_approved$prefName)
count_drugs_drubank <- n_distinct(drugbank_approved$name)

# Prepare data for bar plot comparison
dataset_names <- c("Open_Targets", "Drug_Bank")
unique_drug_counts <- c(3563, 2146)  # hardcoded for plotting

plot_data <- data.frame(
  Dataset = dataset_names,
  Unique_Drugs = unique_drug_counts
)

# Bar plot: Unique drug counts in Open Targets vs DrugBank
ggplot(plot_data, aes(x = Dataset, y = Unique_Drugs)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Unique Drugs in Each Dataset",
    x = "Dataset",
    y = "Number of Unique Drugs"
  ) +
  theme_minimal()

# ------------------------ DEGREE DISTRIBUTION: OPEN TARGETS ------------------------

# Drugs per target (how many drugs target each gene)
nb_drugs_per_target_opentarget <- opentarget_approved %>%
  dplyr::select(approvedSymbol, prefName) %>%
  unique() %>%
  count(approvedSymbol)

# Targets per drug (how many genes each drug targets)
nb_target_per_drug_opentarget <- opentarget_approved %>%
  dplyr::select(approvedSymbol, prefName) %>%
  unique() %>%
  count(prefName)

# Rename columns for clarity
colnames(nb_drugs_per_target_opentarget) <- c("name", "count")
colnames(nb_target_per_drug_opentarget) <- c("name", "count")

# Frequency distribution for each degree (number of interactions)
target_degree_freq <- nb_drugs_per_target_opentarget %>%
  count(count, name = "frequency") %>%
  mutate(category = "Drugs per Target", degree = count) %>%
  dplyr::select(degree, frequency, category)

drug_degree_freq <- nb_target_per_drug_opentarget %>%
  count(count, name = "frequency") %>%
  mutate(category = "Targets per Drug", degree = count) %>%
  dplyr::select(degree, frequency, category)

# Combine both distributions
degree_distribution <- bind_rows(target_degree_freq, drug_degree_freq)

# Plot: Linear scale
ggplot(degree_distribution, aes(x = degree, y = frequency, color = category)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Degree Distribution in Open Targets Bipartite Network",
    x = "Degree (Number of Interactions)",
    y = "Frequency",
    color = "Category"
  ) +
  theme_minimal()

# Plot: Log-log scale
ggplot(degree_distribution, aes(x = degree, y = frequency, color = category)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Degree Distribution in Open Targets Bipartite Network (Log-Log Scale)",
    x = "Degree (log scale)",
    y = "Frequency (log scale)",
    color = "Category"
  ) +
  theme_minimal()

# ------------------------ DEGREE DISTRIBUTION: DRUGBANK ------------------------

# Drugs per target (DrugBank)
nb_drugs_per_target_drugbank <- drugbank_approved %>%
  dplyr::select(name, gene_name) %>%
  unique() %>%
  na.omit() %>%
  count(gene_name)

# Targets per drug (DrugBank)
nb_target_per_drug_drugbank <- drugbank_approved %>%
  dplyr::select(name, gene_name) %>%
  unique() %>%
  na.omit() %>%
  count(name)

# Rename for consistency
colnames(nb_drugs_per_target_drugbank) <- c("name", "count")
colnames(nb_target_per_drug_drugbank) <- c("name", "count")

# Degree distributions
target_degree_freq <- nb_drugs_per_target_drugbank %>%
  count(count, name = "frequency") %>%
  mutate(category = "Drugs per Target", degree = count) %>%
  dplyr::select(degree, frequency, category)

drug_degree_freq <- nb_target_per_drug_drugbank %>%
  count(count, name = "frequency") %>%
  mutate(category = "Targets per Drug", degree = count) %>%
  dplyr::select(degree, frequency, category)

# Combine
degree_distribution <- bind_rows(target_degree_freq, drug_degree_freq)

# Plot: Linear scale
ggplot(degree_distribution, aes(x = degree, y = frequency, color = category)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Degree Distribution in DrugBank Bipartite Network",
    x = "Degree (Number of Interactions)",
    y = "Frequency",
    color = "Category"
  ) +
  theme_minimal()

# Plot: Log-log scale
ggplot(degree_distribution, aes(x = degree, y = frequency, color = category)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Degree Distribution in DrugBank Bipartite Network (Log-Log Scale)",
    x = "Degree (log scale)",
    y = "Frequency (log scale)",
    color = "Category"
  ) +
  theme_minimal()




# Observing how many of the approved drugs target are Protein Coding Genes --------

# Total approved drug targets with OT

# Extract unique gene targets
unique_targets <- unique(opentarget_approved$approvedSymbol)

# Ensure gene symbol column is named consistently
protein_gene_list <- unique(pc_genes$symbol)

# Check which targets are protein-coding
is_protein_coding <- unique_targets %in% protein_gene_list

# Count results
total_targets <- length(unique_targets)
protein_coding_count <- sum(is_protein_coding)
non_protein_coding_count <- total_targets - protein_coding_count

# Display summary
cat("Total unique gene targets:", total_targets, "\n")
cat("Protein-coding gene targets:", protein_coding_count, "\n")
cat("Non-protein-coding gene targets:", non_protein_coding_count, "\n")

non_protein_genes <- unique_targets[!is_protein_coding]
head(non_protein_genes)  # View some of them

# Create a summary dataframe
summary_pc_ot <- data.frame(
  Category = c("Protein-Coding OT", "Non-Protein-Coding OT"),
  Count = c(protein_coding_count, non_protein_coding_count)
) %>%
  mutate(percentage = (Count/sum(Count)*100))

# Extract the pewrcenatge values

protein_coding_percentage <- summary_pc_ot[1,3]
non_protein_coding_percentage <- summary_pc_ot[2,3]



# Plot
ggplot(summary_pc_ot, aes(x = Category, y = percentage, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_minimal() +
  labs(title = "Gene Targets of Approved Drugs",
       y = "Percentage of Unique Gene Targets",
       x = "") +
  scale_fill_manual(values = c("#4CAF50", "#F44336")) +
  theme(legend.position = "none")



# Total approved drug targets with DrugBank

# Extract unique gene targets
unique_targets_db <- unique(drugbank_approved$gene_name)

# Ensure gene symbol column is named consistently
protein_gene_list <- unique(pc_genes$symbol)

# Check which targets are protein-coding
is_protein_coding_db <- unique_targets_db %in% protein_gene_list

# Count results
total_targets_db <- length(unique_targets_db)
protein_coding_count_db <- sum(is_protein_coding_db)
non_protein_coding_count_db <- total_targets_db - protein_coding_count_db

# Display summary
cat("Total unique gene targets DB:", total_targets_db, "\n")
cat("Protein-coding gene targets DB:", protein_coding_count_db, "\n")
cat("Non-protein-coding gene targets DB:", non_protein_coding_count_db, "\n")

non_protein_genes_db <- unique_targets_db[!is_protein_coding_db]
head(non_protein_genes_db)  # View some of them

# Create a summary dataframe
summary_pc_db <- data.frame(
  Category = c("Protein-Coding DB", "Non-Protein-Coding DB"),
  Count = c(protein_coding_count_db, non_protein_coding_count_db)
) %>%
  mutate(percentage = (Count/sum(Count)*100))

# Extract the pewrcenatge values

protein_coding_percentage_db <- summary_pc_db[1,3]
non_protein_coding_percentage_db <- summary_pc_db[2,3]


# Plot
ggplot(summary_pc_db, aes(x = Category, y = percentage, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_minimal() +
  labs(title = "Gene Targets of Approved Drugs",
       y = "Percentage of Unique Gene Targets",
       x = "") +
  scale_fill_manual(values = c("#4CAF50", "#F44336")) +
  theme(legend.position = "none")


#

# Compare the 2 Datasets

# Example data (replace with your actual numbers)

# Create summary dataframe
compare_df <- data.frame(
  Dataset = rep(c("Open Targets", "Drug Bank"), each = 2),
  Category = rep(c("Protein-Coding", "Non-Protein-Coding"), times = 2),
  Count = c(protein_coding_count, non_protein_coding_count,
            protein_coding_count_db, non_protein_coding_count_db),
  Percentage = c (protein_coding_percentage, non_protein_coding_percentage,
                  protein_coding_percentage_db, non_protein_coding_percentage_db)
  
)

ggplot(compare_df, aes(x = Dataset, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_minimal() +
  labs(title = "Comparison of Gene Targets between OpenTargets and DrugBank",
       y = "Percentage of Unique Gene Targets",
       x = "") +
  scale_fill_manual(values = c("#4CAF50", "#F44336"))





# Observing how many of the approved drugs target are fusil Genes --------

# Total approved drug targets with OT

# Extract unique gene targets
unique_targets <- unique(opentarget_approved$approvedSymbol)

# Ensure gene symbol column is named consistently
fusil_gene_list <- unique(fusil$gene_symbol)

# Check which targets are protein-coding
is_fusil_gene <- unique_targets %in% fusil_gene_list

# Count results
total_targets <- length(unique_targets)
fusil_gene_count <- sum(is_fusil_gene)
non_fusil_count <- total_targets - fusil_gene_count

# Display summary
cat("Total unique gene targets:", total_targets, "\n")
cat("Fusil gene targets:", fusil_gene_count, "\n")
cat("Non-fusil gene targets:", non_fusil_count, "\n")

non_fusil_genes <- unique_targets[!is_fusil_gene]
head(non_fusil_genes)  # View some of them

# Create a summary dataframe
summary_fusil_ot <- data.frame(
  Category = c("Fusil OT", "Non-Fusil OT"),
  Count = c(fusil_gene_count, non_fusil_count)
) %>%
  mutate(percentage = (Count/sum(Count)*100))



# Extract the pewrcenatge values

fusil_gene_percentage <- summary_fusil_ot[1,3]
non_fusil_percentage <- summary_fusil_ot[2,3]


# Plot
ggplot(summary_fusil_ot, aes(x = Category, y = percentage, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_minimal() +
  labs(title = "Gene Targets of Approved Drugs",
       y = "Percentage of Unique Gene Targets",
       x = "") +
  scale_fill_manual(values = c("#4CAF50", "#F44336")) +
  theme(legend.position = "none")



# Total approved drug targets with DrugBank

# Extract unique gene targets
unique_targets_db <- unique(drugbank_approved$gene_name)

# Ensure gene symbol column is named consistently
fusil_gene_list <- unique(fusil$gene_symbol)

# Check which targets are protein-coding
is_fusil_gene_db <- unique_targets_db %in% fusil_gene_list

# Count results
total_targets <- length(unique_targets_db)
fusil_gene_count_db <- sum(is_fusil_gene_db)
non_fusil_count_db <- total_targets - fusil_gene_count_db

# Display summary
cat("Total unique gene targets:", total_targets, "\n")
cat("Fusil gene targets:", fusil_gene_count_db, "\n")
cat("Non-fusil gene targets:", non_fusil_count_db, "\n")

non_fusil_genes_db <- unique_targets_db[!is_fusil_gene_db]
head(non_fusil_genes_db)  # View some of them

# Create a summary dataframe
summary_fusil_db <- data.frame(
  Category = c("Fusil DB", "Non-Fusil DB"),
  Count = c(fusil_gene_count_db, non_fusil_count_db)
) %>%
  mutate(percentage = (Count/sum(Count)*100))



# Extract the pewrcenatge values

fusil_gene_percentage_db <- summary_fusil_db[1,3]
non_fusil_percentage_db <- summary_fusil_db[2,3]


# Plot
ggplot(summary_fusil_db, aes(x = Category, y = percentage, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  theme_minimal() +
  labs(title = "Gene Targets of Approved Drugs",
       y = "Percentage of Unique Gene Targets",
       x = "") +
  scale_fill_manual(values = c("#4CAF50", "#F44336")) +
  theme(legend.position = "none")



# Compare the 2 Datasets

# Example data (replace with your actual numbers)

# Create summary dataframe
compare_df <- data.frame(
  Dataset = rep(c("Open Targets", "Drug Bank"), each = 2),
  Category = rep(c("In The FUSIL", "Not In FUSIL"), times = 2),
  Count = c(fusil_gene_count, non_fusil_count,
            fusil_gene_count_db, non_fusil_count_db),
  Percentage = c (fusil_gene_percentage, non_fusil_percentage,
                  fusil_gene_percentage_db, non_fusil_percentage_db)
)

ggplot(compare_df, aes(x = Dataset, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_minimal() +
  labs(title = "Comparison of Gene Targets between OpenTargets and DrugBank",
       y = "Percentage of Unique Gene Targets",
       x = "") +
  scale_fill_manual(values = c("#4CAF50", "#F44336"))




# ------------------------ GENERAL OPEN TARGET DATA OVERVIEW ------------------------

# Count total number of unique drugs and unique indications in Open Targets
length(unique(opentarget_data$prefName)) # There are 4626 different drugs
length(unique(opentarget_data$label))    # There are 2683 different indications

# ------------------------ DRUG STATUS DISTRIBUTION ------------------------

# Count the number of unique drug-status pairs and calculate percentage per status
count_drugs_status <- opentarget_data %>%
  dplyr::select(prefName, status) %>%
  distinct() %>%
  count(status) %>%
  mutate(percentage = n * 100 / sum(n))

# Lollipop plot to visualize number of drugs by status
ggplot(count_drugs_status, aes(x = reorder(status, n), y = n)) +
  geom_segment(aes(xend = status, y = 0, yend = n), color = "gray") +
  geom_point(size = 5, color = "tomato") +
  coord_flip() +
  labs(title = "Number of Drugs by Status",
       x = "Status", y = "Count") +
  theme_minimal()

# ------------------------ INTEGRATION WITH FUSIL DATA ------------------------

# Join Open Targets with FUSIL data by gene symbol and remove NA values
drug_fusil <- opentarget_data %>%
  left_join(fusil, by = c("approvedSymbol" = "gene_symbol")) %>%
  na.omit()

# Number of unique mechanisms of action for genes with FUSIL data
length(unique(drug_fusil$mechanismOfAction))

# Count unique drug–FUSIL category combinations
fusil_count <- drug_fusil %>%
  dplyr::select(prefName, fusil) %>%
  unique() %>%
  count(fusil)

# ------------------------ WITHDRAWN DRUGS & FUSIL RELATIONSHIP ------------------------

# Filter data for withdrawn drugs
withdrawn_drugs <- drug_fusil %>%
  filter(status == c("Withdrawn"))

# Count FUSIL categories among withdrawn drugs
fusil_count_withdrawn <- withdrawn_drugs %>%
  dplyr::select(prefName, fusil) %>%
  unique() %>%
  count(fusil)

# Compute percentage of withdrawn drugs per FUSIL category
percentage_fusil_withdrawn <- fusil_count_withdrawn %>%
  left_join(fusil_count, by = "fusil") %>%
  mutate(percentage = n.x / n.y * 100)

# Set FUSIL category order
percentage_fusil_withdrawn$fusil <- factor(percentage_fusil_withdrawn$fusil, 
                                           levels = c("CL", "DL", "SV", "VP", "VnP"))

# Plot percentage of withdrawn drugs by FUSIL category
ggplot(percentage_fusil_withdrawn, aes(x = fusil, y = percentage, fill = fusil)) +
  geom_bar(stat = "identity") +
  labs(
    title = "FUSIL Category Distribution Among Withdrawn Drug Targets",
    x = "FUSIL Category",
    y = "Percentage of Targets"
  ) +
  theme_minimal()

# ------------------------ COMPLETED DRUGS & FUSIL RELATIONSHIP ------------------------

# Filter data for completed drugs (assumed approved)
completed_drugs <- drug_fusil %>%
  filter(status == c("Completed"))

# Count FUSIL categories among completed drugs
fusil_count_completed <- completed_drugs %>%
  dplyr::select(prefName, fusil) %>%
  unique() %>%
  count(fusil)

# Compute percentage of completed drugs per FUSIL category
percentage_fusil_completed <- fusil_count_completed %>%
  left_join(fusil_count, by = "fusil") %>%
  mutate(percentage = n.x / n.y * 100)

# Set FUSIL category order
percentage_fusil_completed$fusil <- factor(percentage_fusil_completed$fusil, 
                                           levels = c("CL", "DL", "SV", "VP", "VnP"))

# Plot percentage of completed drugs by FUSIL category
ggplot(percentage_fusil_completed, aes(x = fusil, y = percentage, fill = fusil)) +
  geom_bar(stat = "identity") +
  labs(
    title = "FUSIL Category Distribution Among Completed Drug Targets",
    x = "FUSIL Category",
    y = "Number of Targets"
  ) +
  theme_minimal()

# ------------------------ DRUG COUNT PER GENE FOR FUSIL BINS ------------------------

# Count how many completed drugs target each gene
gene_drugs_count <- opentarget_data %>%
  filter(phase == "4" | status == "Completed") %>%
  dplyr::select(approvedSymbol, prefName) %>%
  unique() %>%
  group_by(approvedSymbol) %>%
  tally() %>%
  rename("gene_target_drug_count" = "n")

# Join with FUSIL annotation for those genes
gene_drugs_fusil <- gene_drugs_count %>%
  left_join(fusil, by = c("approvedSymbol" = "gene_symbol")) %>%
  na.omit()

# Order FUSIL categories
gene_drugs_fusil$fusil <- factor(gene_drugs_fusil$fusil, 
                                 levels = c("CL", "DL", "SV", "VP", "VnP"))

# Density plot: distribution of drug counts per gene by FUSIL category
ggplot(gene_drugs_fusil, aes(x = gene_target_drug_count, fill = fusil)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Drug Counts per Gene by FUSIL",
       x = "Drug Count",
       y = "Density") +
  theme_minimal()

# ------------------------ PROPORTIONS OF DRUG PER GENE BY FUSIL ------------------------

# Summarize average drug count per gene in each FUSIL category
gene_drugs_fusil_count <- gene_drugs_count %>%
  left_join(fusil, by = c("approvedSymbol" = "gene_symbol")) %>%
  na.omit() %>%
  group_by(fusil) %>%
  summarise(
    gene_count = n(),                      # Total genes per FUSIL category
    total_drug_count = sum(gene_target_drug_count)  # Total drugs across those genes
  ) %>%
  mutate(proportion = total_drug_count / gene_count)

# Set FUSIL category order
gene_drugs_fusil_count$fusil <- factor(gene_drugs_fusil_count$fusil, 
                                       levels = c("CL", "DL", "SV", "VP", "VnP"))

# Plot: average number of drugs per gene in each FUSIL category
ggplot(gene_drugs_fusil_count, aes(x = fusil, y = proportion, fill = fusil)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Proportions of Drug per Gene per FUSIL Category",
       x = "FUSIL Category",
       y = "Proportion",
       fill = "FUSIL bin") +
  theme_minimal()


# ------------------------ DRUG & INDICATION COUNTS PER GENE BY FUSIL BIN ------------------------

# Filter Open Targets data to include only 'Completed' (approved) drugs
opentarget_data_approved <- opentarget_data %>%
  filter(phase == "4" | status == "Completed")

# Count the number of distinct drugs per gene
gene_drugs_count <- opentarget_data %>%
  filter(phase == "4" | status == "Completed") %>%
  dplyr::select(approvedSymbol, prefName) %>%
  unique() %>%
  group_by(approvedSymbol) %>%
  tally() %>%
  rename("gene_target_drug_count" = "n")

# Count the number of distinct indications per gene
gene_drugs_count_indication <- opentarget_data %>%
  filter(phase == "4" | status == "Completed") %>%
  dplyr::select(approvedSymbol, label) %>%
  unique() %>%
  group_by(approvedSymbol) %>%
  tally() %>%
  rename("indication_count" = "n")

# Combine drug and indication counts per gene
drug_and_indication_count_per_gene <- gene_drugs_count %>%
  full_join(gene_drugs_count_indication, by = "approvedSymbol")

# Join with FUSIL data for each gene and remove rows with missing values
gene_drugs_fusil_df <- drug_and_indication_count_per_gene %>%
  left_join(fusil, by = c("approvedSymbol" = "gene_symbol")) %>%
  na.omit()

# Define order for FUSIL categories
gene_drugs_fusil_df$fusil <- factor(gene_drugs_fusil_df$fusil, 
                                    levels = c("CL", "DL", "SV", "VP", "VnP"))

# Reshape to long format for plotting both drug and indication counts
gene_drugs_fusil_df_long <- gene_drugs_fusil_df %>%
  pivot_longer(cols = c(gene_target_drug_count, indication_count),
               names_to = "count_type",
               values_to = "count")

# Summarize total counts for each FUSIL category and count type
summary_df <- gene_drugs_fusil_df_long %>%
  group_by(fusil, count_type) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  mutate(count_type = recode(count_type,
                             "gene_target_drug_count" = "Drugs per Gene Target",
                             "indication_count" = "Indications per Gene Target"))


summary_df <- summary_df %>%
  group_by(count_type) %>%
  mutate(percentage = (total_count/ sum(total_count))*100)

# Line plot showing total drug and indication counts by FUSIL category

ggplot(summary_df, aes(x = fusil, y = total_count, color = count_type, group = count_type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_text(aes(label = total_count), vjust = -0.5, size = 3, color = "black") +
  scale_y_continuous(labels = function(x) sprintf("%.1f%%", x)) +
  labs(
    title = "Drug and Indication Counts per Target(Gene) per FUSIL Category",
    x = "FUSIL", y = "Count",
    color = "Count Type"
  ) +
  theme_minimal()

# ------------------------ DRUG–INDICATION RATIO PER FUSIL CATEGORY ------------------------

# Calculate average number of drugs per indication per gene for each FUSIL category
# High ratios suggest multiple drugs used for the same indication; lower implies broader use
gene_drugs_fusil_proportion <- gene_drugs_fusil_df %>% 
  mutate(drug_indication_ratio = gene_target_drug_count / indication_count) %>%
  group_by(fusil) %>%
  summarise(
    gene_count = n(),
    avg_ratio = mean(drug_indication_ratio)
  )

# Set FUSIL category order
gene_drugs_fusil_proportion$fusil <- factor(gene_drugs_fusil_proportion$fusil, 
                                            levels = c("CL", "DL", "SV", "VP", "VnP"))

# Bar plot: Average drug-to-indication ratio per gene across FUSIL categories
ggplot(gene_drugs_fusil_proportion, aes(x = fusil, y = avg_ratio, fill = fusil)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Average Indication per drug target gene for each by FUSIL Category",
       x = "FUSIL Category",
       y = "Drugs per Indication (avg per gene)",
       fill = "FUSIL bin") +
  theme_minimal()

# ------------------------ INDICATION-WEIGHTED DRUG IMPACT ------------------------

# Multiply number of drugs and indications to estimate total therapeutic impact
# High values indicate highly targeted genes with wide clinical relevance
gene_drugs_fusil_proportion2 <- gene_drugs_fusil_df %>%
  mutate(total_impact = gene_target_drug_count * indication_count) %>%
  group_by(fusil) %>%
  summarise(
    gene_count = n(),
    Average_total_impact = mean(total_impact)
  )

# Reorder FUSIL levels
gene_drugs_fusil_proportion2$fusil <- factor(gene_drugs_fusil_proportion2$fusil, 
                                             levels = c("CL", "DL", "SV", "VP", "VnP"))

# Plot: Mean of drug × indication product per gene by FUSIL category
ggplot(gene_drugs_fusil_proportion2, aes(x = fusil, y = Average_total_impact, fill = fusil)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Indication-Weighted Drug Targeting per Gene by FUSIL Category",
       x = "FUSIL Category",
       y = "Average Drug × Indication per Gene",
       fill = "FUSIL bin") +
  theme_minimal()


# ------------------------ IMPORT TOP-LEVEL EFO TERMS ------------------------

# Read EFO ontology relationships file (term ↔ ancestor mappings)
efo_ancestor_terms <- read_delim("C:/Users/HP-ssd/Desktop/Short term project2/Analysis script/ontologies scripts/ontologies data wrangling/data_efo/efo_ancestor_terms.txt", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

# ------------------------ CLEAN EFO DATA ------------------------

# Remove 'efo:' prefix and replace underscore with colon in ID fields
efo_ancestor_terms$efo_term_id <- gsub("^efo:", "", efo_ancestor_terms$efo_term_id)
efo_ancestor_terms$efo_ancestor_id <- gsub("^efo:", "", efo_ancestor_terms$efo_ancestor_id)

efo_ancestor_terms$efo_term_id <- sub("_", ":", efo_ancestor_terms$efo_term_id)
efo_ancestor_terms$efo_ancestor_id <- sub("_", ":", efo_ancestor_terms$efo_ancestor_id)



# ------------------------ JOIN EFO TERMS TO OPEN TARGET DATA ------------------------

# Join curated disease labels from EFO with approved drugs in Open Targets
efo_ancestor_terms_joined <- opentarget_approved %>%
  full_join(efo_ancestor_terms, by = c("label" = "efo_term_description")) %>%
  na.omit()


# ------------------------ Seeing drugs distribution per disease class ------------------------

drug_classes_to_remove <- c(
  "urogenital neoplasm", "disease related to transplantation", 
  "connective tissue disease", "poisoning", "syndromic disease", 
  "obstetric disorder", "post-infectious disorder", "chromosomal disorder", "disorder of visual system")


count_disease <- efo_ancestor_terms_joined %>%
  filter(!(efo_ancestor_description %in% drug_classes_to_remove)) %>%
  dplyr::select(prefName, efo_ancestor_description) %>%
  unique() %>%
  count(efo_ancestor_description) %>%
  mutate(sum = sum(n)) %>%
  mutate(percentage = n/sum *100)

#write.csv(count_disease, "C:/Users/HP-ssd/Desktop/count_dx.csv")

ggplot(count_disease, aes(x = reorder(efo_ancestor_description, percentage), y = percentage)) +
  geom_bar(stat = "identity", fill = "tomato") +
  coord_flip() +
  scale_y_log10() +  # log scale for percentage
  labs(
    title = "Drugs Main Disease Category Indicaton",
    x = "Diagnosis Category",
    y = "Log10(Percentage)"
  ) +
  theme_minimal()


# ------------------------ MAP DRUG-TARGET PAIRS TO FUSIL ------------------------

# Extract drug-gene pairs from joined dataset
drug_gene <- efo_ancestor_terms_joined %>%
  dplyr::select(approvedSymbol, prefName) %>%
  unique()

# Join drug-gene pairs with FUSIL data to annotate each target with FUSIL bin
drug_fusil <- fusil %>%
  full_join(drug_gene, by = c("gene_symbol" = "approvedSymbol")) %>%
  na.omit() %>%
  dplyr::select(3, 4, 5)  # dplyr::selects fusil, gene, and drug

# Count number of genes in each FUSIL bin
count_fusil_genes <- drug_fusil %>%
  count(fusil)

# ------------------------ FOR LOOP: FUSIL DISTRIBUTION PER DISEASE CATEGORY ------------------------

# Get unique EFO ancestor terms (disease classes)
unique_drug_indication_class <- unique(efo_ancestor_terms_joined$efo_ancestor_description)

# Create list to hold results by disease class
fusil_by_class <- list()

# Loop through each disease class
for (cls in unique_drug_indication_class) {
  
  # Filter drug-gene pairs for current disease class
  drug_gene <- efo_ancestor_terms_joined %>%
    filter(efo_ancestor_description == cls) %>%
    dplyr::select(approvedSymbol, prefName) %>%
    distinct()
  
  # Join with FUSIL and count occurrences per FUSIL bin
  tmp <- fusil %>%
    full_join(drug_gene, by = c("gene_symbol" = "approvedSymbol")) %>%
    na.omit() %>%
    dplyr::select(3, 4, 5) %>%
    count(fusil) %>%
    left_join(count_fusil_genes, by = "fusil") %>%
    mutate(
      percentage = (n.x / n.y) * 100,  # Proportion of this disease class in each bin
      drug_class = cls
    )
  
  # Set FUSIL order for consistent plotting
  tmp$fusil <- factor(tmp$fusil, levels = c("CL", "DL", "SV", "VP", "VnP"))
  
  # Save to list for combining later
  fusil_by_class[[cls]] <- tmp
  
  # Plot bar chart per class (optional: will auto-display in R environment)
  p <- ggplot(tmp, aes(x = fusil, y = percentage, fill = fusil)) +
    geom_bar(stat = "identity") +
    labs(
      title = paste(cls, "Drugs Distribution among the FUSIL Categories"),
      x = "FUSIL bin",
      y = "Count of Drugs / Targets"
    )
  print(p)
}

# ------------------------ VISUALIZE FUSIL DISTRIBUTION ACROSS DISEASES ------------------------

# Combine all disease class summaries into one data frame
combined_fusil <- bind_rows(fusil_by_class)

# Optional: Filter out some drug classes to reduce clutter in plots
drug_classes_to_keep <- c("urogenital neoplasm", "premature aging syndrome", "chromosomal disorder", 
                          "obstetric disorder", "upper digestive tract disorder", 
                          "disease related to transplantation", "idiopathic disease", 
                          "otorhinolaryngologic disease", "acute disease", "poisoning", 
                          "radiation-induced disorder", "ulcer disease", "perinatal disease", 
                          "obstetric disorder", "post-infectious disorder")

# Reorder and filter
combined_fusil <- combined_fusil %>%
  group_by(fusil) %>%
  mutate(drug_class = fct_reorder(drug_class, percentage)) %>%
  filter(!(drug_class %in% drug_classes_to_keep))

# ------------------------ PLOT 1: FACET BY FUSIL BIN ------------------------

ggplot(combined_fusil, aes(x = drug_class, y = percentage, fill = drug_class)) +
  geom_bar(stat = "identity", nrow = 1) +
  facet_wrap(~ fusil) +
  coord_flip() +
  labs(
    title = "Distribution of Disease Classes Within Each FUSIL Bin",
    x = "Disease Class",
    y = "Percentage of Drug-Target Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  geom_text(aes(label = round(percentage, 1)), hjust = -0.2, size = 3)

# ------------------------ PLOT 2: FACET BY DISEASE CLASS ------------------------

ggplot(combined_fusil, aes(x = fusil, y = percentage, fill = drug_class)) +
  geom_bar(stat = "identity", nrow = 1) +
  facet_wrap(~ drug_class) +
  coord_flip() +
  labs(
    title = "Distribution of Disease Classes Within Each FUSIL Bin",
    x = "FUSIL bin",
    y = "Percentage of Drug-Target Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  geom_text(aes(label = round(percentage, 1)), hjust = -0.2, size = 3)
