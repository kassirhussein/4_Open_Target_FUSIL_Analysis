# -------------------------------------------------------------------------
# Project      : Ontologies
# Script       : GO_obo_to_txt.R
# Description  : Imports the GO Ontology (OBO format) and extracts key 
#                relationships (parents, children, ancestors) into plain text.
# Author       : Hussein
# Contributor  : Dr. Pilar Cacheiro
# Credit       : Code structure and functions adapted from 
#                https://github.com/pilarcacheiro/Ontologies
# -------------------------------------------------------------------------







# install / load libraries ------------------------------------------------

if (!require("ontologyIndex")) install.packages("ontologyIndex")
library("ontologyIndex")

if (!require("dplyr")) install.packages("dplyr")
library("dplyr")

if (!require("tidyr")) install.packages("tidyr")
library("tidyr")

if (!require("Hmisc")) install.packages("Hmisc")
library("Hmisc")

# Load the MONDO Data ---------------------------------

go <- get_ontology("https://purl.obolibrary.org/obo/go/go-basic.obo")

#mondo id and description

go_description <- data.frame(
  go_term = as.character(row.names(as.data.frame(go[[2]]))),
  go_description = as.character(as.data.frame(go[[2]])[, 1]), 
  stringsAsFactors = F)


## parents, children and all the ancestors for a given go term as a list

parents <- go[[3]]
children <-  go[[4]]
ancestors <- go[[5]]



# set of functions to convert these lists to data frames ------------------


## get parental nodes 

get_parent_nodes <- function(parents) {
  
  go_term <- list()
  go_parents <-list()
  
  for(i in 1:length(names(parents))){
    
    go_term [[i]] <- names(parents)[i]
    go_parents[[i]] <- paste(parents[[i]], collapse=",")
  }
  
  go_parents <- data.frame(
    go_term = do.call(rbind,go_term), 
    go_parents = do.call(rbind, go_parents), stringsAsFactors = F) %>%
    mutate(go_parents = strsplit(as.character(go_parents), ",")) %>%
    unnest(go_parents) %>%
    filter(go_term != go_parents)
  
  return(go_parents)
  
}


## get children nodes 

get_children_nodes <- function(children) {
  
  go_term <- list()
  go_children <- list()
  
  for(i in 1:length(names(children))){
    
    go_term [[i]] <- names(children)[i]
    go_children[[i]] <- paste(children[[i]], collapse=",")
  }
  
  go_children <-  data.frame(
    go_term = do.call(rbind, go_term),
    go_children = do.call(rbind, go_children), stringsAsFactors = F) %>%
    mutate(go_children = strsplit(as.character(go_children), ",")) %>%
    unnest(go_children) %>%
    filter(go_term != go_children)
  
  return(go_children)
}


## get all the ancestors/ top levels

get_ancestor_nodes <- function(ancestors) {
  
  go_term <- list()
  go_ancestors <- list()
  
  for(i in 1:length(names(ancestors))) {
    
    go_term [[i]] <- names(ancestors)[i]
    go_ancestors[[i]] <- paste(ancestors[[i]], collapse=",")
  }
  
  go_ancestors <- data.frame(
    go_term = do.call(rbind, go_term),
    go_ancestors = do.call(rbind, go_ancestors),stringsAsFactors = F) %>%
    mutate(go_ancestors = strsplit(as.character(go_ancestors), ",")) %>%
    unnest(go_ancestors) %>%
    filter(go_term != go_ancestors )
  
  return(go_ancestors)
}

## get top levels of the ontology (physiological systems)

go_toplevels <- get_parent_nodes(parents) %>%
  filter(go_parents == "GO:0008150") %>%
  left_join(go_description,by = "go_term") %>%
  dplyr::select(go_term, go_description)

go_toplevels_phenotypic_abnormalities_only <- get_parent_nodes(parents) %>%
  filter(go_parents == "GO:0009987") %>%
  left_join(go_description,by = "go_term") %>%
  dplyr::select(go_term,go_description)


# apply functions ---------------------------------------------------------


go_parental_nodes <- get_parent_nodes(parents)

go_children_nodes <- get_children_nodes(children)

go_ancestor_nodes <- get_ancestor_nodes(ancestors)


# export files ------------------------------------------------------------


go_dir <- "./data_go/"


files_to_export <- list(go_description, 
                        go_parental_nodes,
                        go_children_nodes,
                        go_ancestor_nodes,
                        go_toplevels,
                        go_toplevels_phenotypic_abnormalities_only)

names(files_to_export) <- Cs(go_description,
                             go_parental_nodes,
                             go_children_nodes,
                             go_ancestor_nodes,
                             go_toplevels,
                             go_toplevels_phenotypic_abnormalities_only)


for (i in 1:length(files_to_export)){
  
  write.table(files_to_export[[i]],
              paste0(go_dir, names(files_to_export)[i], ".txt"), 
              quote = F, sep = "\t", row.names = FALSE)
  
}




# Create the directory if it doesn't exist
if (!dir.exists(go_dir)) {
  dir.create(go_dir, recursive = TRUE)
}

# Now write the files
for (i in 1:length(files_to_export)) {
  write.table(files_to_export[[i]],
              paste0(go_dir, names(files_to_export)[i], ".txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE)
}






# --- Code above adapted from https://github.com/pilarcacheiro/Ontologies ---
# Original logic by Dr. Pialr Cacheiro






# MONDO_TOP_LEVEL_DATA_WRANGLING ----------------------


library(tidyverse)

go_description <- read.delim("./data_go/go_description.txt")
go_parental_nodes <- read.delim("./data_go/go_parental_nodes.txt")
go_children_nodes <- read.delim("./data_go/go_children_nodes.txt")
go_ancestor_nodes <- read.delim("./data_go/go_ancestor_nodes.txt")
go_toplevels <- read.delim("./data_go/go_toplevels.txt")
go_toplevels_phenotypic_abnormalities_only <- read.delim("./data_go/go_toplevels_phenotypic_abnormalities_only.txt")


go_ancestor_filtered <- go_ancestor_nodes %>%
  filter(go_ancestors %in% go_toplevels_phenotypic_abnormalities_only$go_term) %>%
  left_join(go_description, by = "go_term") %>%
  left_join(go_toplevels_phenotypic_abnormalities_only, by = c("go_ancestors" = "go_term"))


# Rearrange and rename columns

go_ancestor_arranged <- go_ancestor_filtered %>%
  rename("go_term_id" = "go_term") %>%
  rename("go_ancestor_id" = "go_ancestors") %>%
  rename("go_term_description" = "go_description.x") %>%
  rename("go_ancestor_description" = "go_description.y")


# Export table

go_dir <- "./data_go/"

write.table(go_ancestor_arranged,
            paste0(go_dir, "go_ancestor_terms.txt"),  
            quote = FALSE, sep = "\t", row.names = FALSE)
