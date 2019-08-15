### Get GO terms for a subset of genes found to be imprinted in repro and sterile BB

setwd("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Final_imprinting_model/GO_analysis")
library(readr)

paternal_just_sterile_all_data <- read_csv("paternal_just_sterile_all_data.csv")
paternal_just_repro_all_data <- read_csv("paternal_just_repro_all_data.csv")
maternal_just_repro_all_data <- read_csv("maternal_just_repro_all_data.csv")
maternal_just_strile_all_data <- read_csv("maternal_just_strile_all_data.csv")
greater90maternal_both_alldata <- read_delim("greater90maternal_both_alldata.txt", 
                                             "\t", escape_double = FALSE, trim_ws = TRUE)



Bumble_bee_ensemble_GO_terms <- read_delim("~/Dropbox/PhD Folder/Virtual Lab Book/Leuven_Samples/Transcription/Differential_expression/Bumble_bee_ensemble_GO_terms.txt", 
                                           "\t", escape_double = FALSE, col_names = FALSE, 
                                           trim_ws = TRUE)
colnames(Bumble_bee_ensemble_GO_terms) <- c("Symbol", "go_id")

paternal_just_sterile_all_data_GO <- merge(paternal_just_sterile_all_data, Bumble_bee_ensemble_GO_terms,
                                           by = "Symbol")

paternal_just_repro_all_data_GO <- merge(paternal_just_repro_all_data, Bumble_bee_ensemble_GO_terms,
                                           by = "Symbol")

maternal_just_repro_all_data_GO <- merge(maternal_just_repro_all_data, Bumble_bee_ensemble_GO_terms,
                                           by = "Symbol")

maternal_just_strile_all_data_GO <- merge(maternal_just_strile_all_data, Bumble_bee_ensemble_GO_terms,
                                         by = "Symbol")
maternal_top <- merge(Bumble_bee_ensemble_GO_terms, greater90maternal_both_alldata,
                      by = "Symbol")


# File from Alun's Trinotate annotation of the BB genome:
#go_ids <- read_delim("go_ids.xls", "\t", escape_double = FALSE, trim_ws = TRUE)
#library(tidyr)
#go_ids_sorted <- gather(go_ids)
#bio_processes <- go_ids_sorted[grep("biological", go_ids_sorted$value),]
#write.csv(bio_processes, file="all_gos_with_description.csv") # In excel used the 'text to columns' to make new fields

all_gos_with_description <- read_csv("all_gos_with_description.csv")


mat_top_final <- merge(maternal_top, all_gos_with_description, by = "go_id")
mat_top_final1 <- unique(mat_top_final)
write.csv(mat_top_final1, file= "top0.9_maternal_both_with_go.csv")


mat_repro_final <- merge(maternal_just_repro_all_data_GO, all_gos_with_description, by = "go_id")
mat_repro_final1 <- unique(mat_repro_final)
write.csv(mat_repro_final1, file = "maternal_just_repro_with_go.csv")

pat_repro_final <- merge(paternal_just_repro_all_data_GO, all_gos_with_description, by = "go_id")
pat_repro_final1 <- unique(mat_repro_final)
write.csv(pat_repro_final1, file = "paternal_just_repro_with_go.csv")

mat_sterile_final <- merge(maternal_just_strile_all_data_GO, all_gos_with_description, by = "go_id")
mat_sterile_final1 <- unique(mat_sterile_final)
write.csv(mat_sterile_final1, file = "maternal_just_sterile_with_go.csv")

pat_sterile_final <- merge(paternal_just_sterile_all_data_GO, all_gos_with_description, by = "go_id")
pat_sterile_final1 <- unique(pat_sterile_final)
write.csv(pat_sterile_final1, file = "paternal_just_sterile_with_go.csv")