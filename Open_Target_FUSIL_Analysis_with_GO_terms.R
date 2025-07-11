# Title: "A Systems-Level View of Drug Targeting: From Molecular Action to Disease Classes"


# Load tidyverse for data manipulation and visualization
library(tidyverse)

# ------------------------ DATA IMPORT ------------------------

# Import MHRA side effect dataset and remove rows with NA
mhra_se_data <- read_csv("C:/Users/HP-ssd/Desktop/Short term project2/mhra/mhra_se_data.csv", 
                         col_types = cols(...1 = col_skip()))
mhra_se_data <- mhra_se_data %>%
  na.omit()

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
  filter(status == "Completed")



library(ontologyIndex)

# Stratifyin Drug Target with GO SLIM -----------------

 # Step 1 : Get GO terms associtaed with Tragets/Genes using Biomart

open_target_gene_list <- unique(opentarget_data$approvedSymbol)

library(biomaRt)

ensembl = useEnsembl(biomart = "genes")

datasets = listDatasets(ensembl)


ensembl_human = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)


attributes = listAttributes(mart = ensembl_human)


opentarget_gene_go_term <- getBM(attributes = c("external_gene_name",
                                                "go_id",
                                                "name_1006",
                                                "definition_1006"),
                                 filters = "external_gene_name",
                                 values = open_target_gene_list,
                                 mart = ensembl_human)


opentarget_gene_go_term <- opentarget_gene_go_term[!(opentarget_gene_go_term$go_id == "" | opentarget_gene_go_term$name_1006 == "" | opentarget_gene_go_term$definition_1006 == ""), ]

# Step 2: Get GO Ancestor Term and map them Genes GO terms using Go_Term_OBO_to_TXT.R function

opentarget_gene_go_term_joined <- opentarget_gene_go_term %>%
  left_join(go_ancestor_arranged, by =  c ("go_id" = "go_term_id")) %>%
  na.omit() %>%
  dplyr::select(1,2,5,6,7)

# step 3 : Join with Open Target Data


opentarget_approved_goterm <- opentarget_approved %>%
  left_join(opentarget_gene_go_term_joined, by = c("approvedSymbol" = "external_gene_name")) %>%
  unique()




# ------------------------ MAP DRUG-TARGET PAIRS TO FUSIL ------------------------

# Extract drug-gene pairs from joined dataset
drug_gene <- opentarget_approved_goterm %>%
  dplyr::select(approvedSymbol, prefName) %>%
  unique()

# Join drug-gene pairs with FUSIL data to annotate each target with FUSIL bin
drug_fusil <- fusil %>%
  full_join(drug_gene, by = c("gene_symbol" = "approvedSymbol")) %>%
  na.omit() %>%
  dplyr::select(3, 4, 5)  # selects fusil, gene, and drug

# Count number of genes in each FUSIL bin
count_fusil_genes <- drug_fusil %>%
  count(fusil)

# ------------------------ FOR LOOP: FUSIL DISTRIBUTION PER GO CATEGORY ------------------------

# Get unique EFO ancestor terms (disease classes)
unique_go_indication_class <- unique(opentarget_approved_goterm$go_ancestor_description)

# Create list to hold results by disease class
fusil_by_class <- list()

# Loop through each disease class
for (cls in unique_go_indication_class) {
  
  # Filter drug-gene pairs for current disease class
  drug_gene <- opentarget_approved_goterm %>%
    filter(go_ancestor_description == cls) %>%
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
      y = "Count of Drugs / Targets (Percentages)"
    )
  print(p)
}

# ------------------------ VISUALIZE FUSIL DISTRIBUTION ACROSS DISEASES ------------------------

# Combine all disease class summaries into one data frame
combined_fusil2 <- bind_rows(fusil_by_class)


# Optional: Filter out some drug classes to reduce clutter in plots
drug_classes_to_remove <- c("vesicle targeting", "nutrient storage","cellular senescence" , "cell aggregation" ,
                            "maintenance of location in cell", "intercellular transport", "intermediate filament-based process",
                            "membrane docking" ,"cellular detoxification","protein folding","establishment or maintenance of cell polarity" ,
                            "exocytic process", "microtubule organizing center localization", "cell division","cellular process involved in reproduction in multicellular organism",
                            "cell-cell fusion","syncytium formation" )

# Reorder and filter
combined_fusil2 <- combined_fusil2 %>%
  group_by(fusil) %>%
  mutate(drug_class = fct_reorder(drug_class, percentage)) %>%
  filter(!(drug_class %in% drug_classes_to_remove))

# ------------------------ PLOT 1: FACET BY FUSIL BIN ------------------------

ggplot(combined_fusil2, aes(x = drug_class, y = percentage, fill = drug_class)) +
  geom_bar(stat = "identity", nrow = 1) +
  facet_wrap(~ fusil) +
  coord_flip() +
  labs(
    title = "Distribution of Gene/Target GO TERM Classes Within Each FUSIL Bin",
    x = "GO Term Class",
    y = "Percentage of Drug-Target Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  geom_text(aes(label = round(percentage, 1)), hjust = -0.2, size = 3)

# ------------------------ PLOT 2: FACET BY DISEASE CLASS ------------------------

ggplot(combined_fusil2, aes(x = fusil, y = percentage, fill = drug_class)) +
  geom_bar(stat = "identity", nrow = 1) +
  facet_wrap(~ drug_class) +
  coord_flip() +
  labs(
    title = "Distribution of Gene/Target GO TERM Classes Within Each FUSIL Bin",
    x = "FUSIL bin",
    y = "Percentage of Drug-Target Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  geom_text(aes(label = round(percentage, 1)), hjust = -0.2, size = 3)



