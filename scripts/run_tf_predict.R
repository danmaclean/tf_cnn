#!/usr/bin/env Rscript

library(keras)
source("scripts/utilities.R")



pairs <- pair_dataframe(         tf_file = "lib/Ath_TF_list", 
                                 interaction_file = "lib/AtRegNet",
                                 non_tf_file = "lib/all_agis.txt"
                                 )
pair_data <- make_pair_data("lib/normalised_data.csv", 
                            pairs, 
                            probe_info_file = "lib/affy_ATH1_array_elements-2010-12-20.txt")


affy_probes <- data.frame(
                    agi_TF = pairs$TF,
                    agi_Target = pairs$Target,
                    affy_TF = agi_to_affy_mapping("lib/affy_ATH1_array_elements-2010-12-20.txt")[[pairs$TF]], 
                    affy_Target = agi_to_affy_mapping("lib/affy_ATH1_array_elements-2010-12-20.txt")[[pairs$Target]]
)

affy_probes<- affy_probes[complete.cases(affy_probes), ]

model <- load_model_hdf5("data/convnet_model.hdf5")
result_probs <- predict_proba(model, pair_data)
result_class <- predict_classes(model, pair_data)

predictions <- data.frame(
    class_prediction = result_class,
    class_likelihood = result_probs
                     )
output <- dplyr::bind_cols(affy_probes, predictions)

filename <- "data/predictions.csv"

readr::write_csv(output, filename, col_names = FALSE)