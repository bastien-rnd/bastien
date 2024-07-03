
# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(ggsci)

# DATA --------------------------------------------------------------------

## IMPORTATION
metaphlan <- read.csv("metaphlan-4.1.0_vJun23_CHOCOPhlAnSGB_202403_merged_abundances.tsv", sep = "\t")


## MISE EN FORME
metaphlan <- as_tibble(metaphlan)

species_all <- metaphlan %>% dplyr::filter(str_detect(clade_name, "s__") & !str_detect(clade_name, "t__"))

species_ird <- species_all %>% dplyr::select(contains("clade_name") | contains("_IRD_"))

species_ird_filter <- species_ird %>% filter(apply(dplyr::select(., -clade_name), 1, function(row) any(row > 0.10)))

#Créer une table taxonomique
taxonomy <- species_ird_filter %>% 
  separate(clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|", fill = "right") %>%
  dplyr::select(1:7)

#Créer une fonction pour enlever s__ du nom d'espèce
clean_species_name <- function(species) {
  sub("^s__", "", species)
}

#Créer une colonne clean à partir de la fonction
taxonomy <- taxonomy %>%
  mutate(Name = clean_species_name(Species))

#Réorganiser les colonnes
taxonomy <- taxonomy %>% dplyr::select(Name, Kingdom, Phylum, Class, Order, Family, Genus, Species, everything())

species_ird_filter <- bind_cols(taxonomy %>% dplyr::select(Name), species_ird_filter)
species_ird_filter <- species_ird_filter %>% dplyr::select(-2)

metadata <- readxl::read_excel("metadata_metagenomes.xlsx")
metadata <- as_tibble(metadata)

#Phyloseq fonctionne avec 3 inputs
# Count
# Taxonomy (rownames = ligne count)
# Metadata (rownames = col count)

# Ajoute les rownames pour n'avoir que des données d'intérêts à traiter
species_ird_filter <- species_ird_filter %>%
  tibble::column_to_rownames("Name") 

taxonomy <- taxonomy %>%
  tibble::column_to_rownames("Name") 

metadata <- metadata %>%
  tibble::column_to_rownames("sample") 


species_mat <- as.matrix(species_ird_filter)
taxa_mat <- as.matrix(taxonomy)


#Création de l'objet phyloseq
OTU = otu_table(species_mat, taxa_are_rows = TRUE)
TAX = tax_table(taxa_mat)
samples = sample_data(metadata)

carbom <- phyloseq(OTU, TAX, samples)

plot_bar(carbom, fill = "Genus") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  facet_wrap(~Genus)

plot_bar(carbom, fill = "Genus") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +


plot_bar(carbom, x="Genus", fill = "Genus", facet_grid = matrix~process) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

plot_heatmap(carbom, method = "NMDS", distance = "bray")




