library(magrittr)

# Makes array of pairs of gene expression profiles from a text file
make_pairs <- function(f, positive_pairs){

  #turn expression file to matrix
  dat <- as.matrix(readr::read_tsv(f))
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]

  #filter out rows matching control probes
  good_probes <- grep("^\\d+.*", rownames(dat), perl = TRUE)
  good_probes <- rownames(dat)[good_probes]
  dat <- dat[good_probes, ]

  #select the pairs of probes we would like to be the negative control set
  neg_tags <- c()
  while ( length(neg_tags) < length(positive_pairs$tag)) {
      new_tag <- paste0( sample(rownames(dat),1 ),
                         "-",
                         sample(rownames(dat), 1)
                        )
      if (! new_tag %in% positive_pairs$tag){
        neg_tags <- c(neg_tags, new_tag)
      }
  }

  negative_pairs <- data.frame(tag = neg_tags) %>%
    tidyr::separate(tag, into = c("TF", "Target"), sep = "-", remove = F) %>%
    dplyr::mutate(class = FALSE)

  all_pairs <- dplyr::bind_rows(positive_pairs, negative_pairs)

  ##set up the tensor to hold the example pair expression data
  x <- array(NA, dim = c(nrow(all_pairs), 2, ncol(dat) ))

  #build the tensor
  for (i in 1:nrow(all_pairs)){
    x[i,,] <- as.numeric(rbind(dat[all_pairs$TF[i], ],dat[all_pairs$Target[i], ]))
  }


  list(
    x = x,
    pair_info = all_pairs,
    y = all_pairs$class
  )

}

#extracts true target pairs from AtRegNet http://agris-knowledgebase.org/AtTFDB/
#and maps to tair probeset ids https://www.arabidopsis.org/download_files/Microarrays/Affymetrix/affy_ATH1_array_elements-2010-12-20.txt
get_positive <- function(f, locus_to_probe){

  interactions <- readr::read_tsv(f) %>%
   dplyr::select(TFLocus, TargetLocus, Confirmation, TargetType) %>%
     dplyr::filter(Confirmation %in% c("confirmed", "Confirmed")) %>%
     dplyr::filter(TargetType %in% c("direct","Direct", "Directly")) %>%
    dplyr::mutate(TF = toupper(TFLocus)) %>%
    dplyr::mutate(Target = toupper(TargetLocus)) %>%
    dplyr::select(TF, Target)

 mapping <-  readr::read_tsv(locus_to_probe) %>%
   dplyr::select(array_element_name, locus)

 dplyr::left_join(interactions, mapping, by =c("TF" = "locus") ) %>%
   dplyr::left_join(mapping, by = c("Target" = "locus" ) ) %>%
   dplyr::select(TF = array_element_name.x, Target = array_element_name.y) %>%
   dplyr::filter( ! is.na( TF), ! is.na(Target) ) %>%
   dplyr::mutate(tag = paste0(TF, "-", Target)) %>%
   dplyr::mutate(class = TRUE)  %>%
   dplyr::filter(! stringr::str_detect(TF, "AFFX"), ! stringr::str_detect(Target, "AFFX"))

}

make_data_set <- function(expression_data, regulatory_info, locus_to_probe_map ){


  positive_pairs <- get_positive(regulatory_info, locus_to_probe_map)
  d <- make_pairs(expression_data, positive_pairs)
  saveRDS(d, file = "data/arab_TF.RDS")
}

make_data_set("lib/normalised_data.csv", "lib/AtRegNet", "lib/affy_ATH1_array_elements-2010-12-20.txt")
