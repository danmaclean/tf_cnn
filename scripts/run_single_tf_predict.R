#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(keras)
suppressMessages(library(log4r))

dir.create("tmp", showWarnings = FALSE)
logger <- create.logger()
logfile(logger) <- file.path("tmp/skipped.log")
level(logger) <- "INFO"

source("scripts/utilities.R")
tf_agi <- get_tfs( file = "lib/Ath_TF_list" )
tf_agi <- tf_agi[as.integer(args[1])]

pairs <- specific_pair_dataframe(tf_agi, 
                                 interaction_file = "lib/AtRegNet",
                                 non_tf_file = "lib/all_agis.txt"
                                 )
pair_data <- make_pair_data("lib/normalised_data.csv", 
                            pairs, 
                            probe_info_file = "lib/affy_ATH1_array_elements-2010-12-20.txt")

if ( is.null(pair_data) ){
  info(logger, tf_agi) 
  stop("empty pair data. Skipping predict.")
}

affy_probes <- data.frame(
                    agi_TF = pairs$TF,
                    agi_Target = pairs$Target,
                    affy_TF = agi_to_affy_mapping("lib/affy_ATH1_array_elements-2010-12-20.txt")[[pairs$TF]], 
                    affy_Target = agi_to_affy_mapping("lib/affy_ATH1_array_elements-2010-12-20.txt")[[pairs$Target]]
)

affy_probes<- affy_probes[complete.cases(affy_probes), ]



model <- suppressMessages(load_model_hdf5("data/convnet_model.hdf5"))
suppressMessages(result_probs <- predict_proba(model, pair_data))
suppressMessages(result_class <- predict_classes(model, pair_data))

predictions <- data.frame(
    class_prediction = result_class,
    class_likelihood = result_probs
                     )
output <- dplyr::bind_cols(affy_probes, predictions)

filename <- paste0("tmp/predictions_", tf_agi, ".csv")

readr::write_csv(output, filename, col_names = FALSE)