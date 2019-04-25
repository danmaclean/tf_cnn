predict_interactions <- function(model, data){
  
  
}

library(magrittr)
library(hashmap)
# Makes 4d array of pairs of gene expression profiles from a text file of pair names for use in convnet model
# expression_file = expression file, with row_names in first column
# pairs = data frame with columns TF, Target describing pairs of interactions,
# scale = whether to scale expression values
make_pair_data <- function(expression_file, pairs, scale=TRUE, convert_agi = TRUE, probe_info_file = "lib/affy_ATH1_array_elements-2010-12-20.txt"){

  if (convert_agi){
    pairs <- data.frame(TF = agi_to_affy_mapping(probe_info_file)[[pairs$TF]], 
                        Target = agi_to_affy_mapping(probe_info_file)[[pairs$Target]],
                        stringsAsFactors = FALSE
                        )
    pairs <- pairs[complete.cases(pairs), ]
    #if dataframe is empty (TF has no affy or no targets have affy)
    if (nrow(pairs) == 0){
      return(NULL)
    }
  }

  #turn expression file to matrix
  dat <- suppressWarnings(as.matrix(readr::read_tsv(expression_file, col_types = readr::cols() )))
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  #filter out rows matching control probes
  good_probes <- grep("^\\d+.*", rownames(dat), perl = TRUE)
  good_probes <- rownames(dat)[good_probes]
  dat <- dat[good_probes, ]
  
  ##set up the tensor to hold the example pair expression data
  x <- array(0, dim = c(nrow(pairs), 2, ncol(dat), 1 ))
  
  #build the tensor
  for (i in 1:nrow(pairs)){ #sample
    if(scale){
      scale_a <- scale(as.numeric(dat[pairs$TF[i], ]))
      scale_b <- scale(as.numeric(dat[pairs$Target[i], ]))

      for(ex in 1:ncol(dat) ){ #experiment column
        x[i,1,ex,1] <-scale_a[ex]  #    small_cels[1,i]
        x[i,2,ex,1] <-scale_b[ex] #  small_cels[2,i]
      }
      
      #cat(dat[all_pairs$TF[i], ])
      #cat(dat[all_pairs$Target[i], ])
      
    }
    else{
      for(ex in 1:ncol(dat) ){ #experiment column
        x[i,ex,ex,1] <-dat[ pairs$TF[i], ex ]  #    small_cels[1,i]
        x[1,ex,ex,2] <-dat[ pairs$Target[i], ex ] #  small_cels[2,i]
      }
      
    }
  }
  
  return(x)
  
}

#generates vector of A.thal TFs from http://planttfdb.cbi.pku.edu.cn/download.php
get_tfs <- function(file = "lib/Ath_TF_list"){
  tfs <- suppressWarnings(readr::read_tsv(file, col_types = readr::cols() ))
  return(unique(tfs$Gene_ID))
}

## generates vector of non-tfs, all unique protein coding agis 
get_non_tfs <- function(file = "lib/all_agis.txt"){
    d <- suppressWarnings(readr::read_tsv(file, col_types = readr::cols(),  col_names = FALSE ))
    return(unique(d$X1))
}

#gets all known confirmed interactions from AtRegNet
get_known_interactions <- function(file = "lib/AtRegNet"){
  suppressWarnings(readr::read_tsv(file, col_types = readr::cols() )) %>%
    dplyr::select(TFLocus, TargetLocus, Confirmation, TargetType) %>%
    dplyr::filter(Confirmation %in% c("confirmed", "Confirmed")) %>%
    dplyr::filter(TargetType %in% c("direct","Direct", "Directly")) %>%
    dplyr::mutate(TF = toupper(TFLocus)) %>%
    dplyr::mutate(Target = toupper(TargetLocus)) %>%
    dplyr::select(TF, Target)
  
}

#makes all possible pairs from known TFs in HQ interactions and non TFs
pair_dataframe <- function(interaction_file = "lib/AtRegNet", tf_file = "lib/Ath_TF_list", non_tf_file = "lib/all_agis.txt"){
  tfs <- get_tfs(file = tf_file)
  non_tfs <- get_non_tfs(file = non_tf_file)
  known_pairs <- get_known_interactions(interaction_file) %>% dplyr::mutate(tag = paste0(TF, "-", Target))
  
  expand.grid(TF = tfs, Target = non_tfs, stringsAsFactors = FALSE) %>%
    dplyr::arrange(TF) %>%
    dplyr::mutate( tag = paste0(TF, "-", Target)) %>%
    dplyr::filter( ! tag %in% known_pairs$tag ) %>%
    dplyr::select(TF, Target)
  
}


#makes all possible pairs from known TFs in HQ interactions and non TFs
specific_pair_dataframe <- function(tf_agi, interaction_file = "lib/AtRegNet", non_tf_file = "lib/all_agis.txt"){

  non_tfs <- get_non_tfs(file = non_tf_file)
  known_pairs <- get_known_interactions(interaction_file) %>% dplyr::mutate(tag = paste0(TF, "-", Target))
  
  expand.grid(TF = tf_agi, Target = non_tfs, stringsAsFactors = FALSE) %>%
    dplyr::arrange(TF) %>%
    dplyr::mutate( tag = paste0(TF, "-", Target)) %>%
    dplyr::filter( ! tag %in% known_pairs$tag ) %>%
    dplyr::select(TF, Target)
  
}

select_probes <- function(file){
  suppressWarnings(readr::read_tsv(file, col_types = readr::cols())) %>%
    dplyr::filter(is_control == "no", ! chromosome %in% c("M", "C")) %>%
    dplyr::select(array_element_name, locus)
}

#map array element to locus and vice versa
affy_to_agi_mapping <- function(file = "lib/affy_ATH1_array_elements-2010-12-20.txt"){
  d <- select_probes(file)
  return(hashmap(d$array_element_name, d$locus))
}
agi_to_affy_mapping <- function(file = "lib/affy_ATH1_array_elements-2010-12-20.txt"){
  d <- select_probes(file)
  return(hashmap(d$locus,d$array_element_name))
}

#converts a dataframe of predictions to list with members edges -> an edgelisted datafram
# and nodes -> a unique list of nodes
to_edgelist <- function(df, prob_cutoff = 0.95){
  
  df <- df %>% dplyr::filter(X6 >= prob_cutoff)
  
  tfs <- df %>% dplyr::distinct(X1) %>% dplyr::rename(agi = X1)
  targs <- df %>% dplyr::distinct(X2) %>% dplyr::rename(agi = X2)
  nodes <- dplyr::full_join(tfs, targs, by = "agi" )  %>%
    tibble::rowid_to_column("id")
  
  edges <- df %>% 
    dplyr::left_join(nodes, by = c("X1" = "agi")) %>% 
    dplyr::rename( tf = id) %>% 
    dplyr::left_join(nodes, by =c("X2" = "agi")) %>% 
    dplyr::rename(target = id) %>%
    dplyr::rename(tf_agi = X1, target_agi = X2, tf_affy = X3, target_affy = X4, keras_class = X5, keras_prob = X6)
  
  return(
    list(nodes = nodes,
         edges = edges
    )
  )
}


#given a pair of expression values, will return the class activation values for that pair as a 116 element vector
get_class_activation_values <- function(expression_pair, model, layer = "separable_conv2d_2"){
  
  last_conv_layer <- model %>% get_layer(layer)
  real_interaction_output <- model$output[,1]
  grads <- k_gradients(real_interaction_output, last_conv_layer$output)[[1]]
  pooled_grads <- k_mean(grads, axis = c(1,2,3))
  
  iterate <- k_function(list(model$input),
                        list(pooled_grads, last_conv_layer$output[1,,,])
  )
  
  c(pooled_grads_value, conv_layer_output_value) %<-% iterate(list(expression_pair))
  
  for (i in 1:8){
    conv_layer_output_value[,,i] <- conv_layer_output_value[,,i] * pooled_grads_value[[i]]
  }
  
  activations <- apply(conv_layer_output_value, c(1,2), mean)
  return(activations)
}


#normalises the activation vector
normalise_activations <- function(acts){
  acts <- scale( acts[1,] ) #scale around 0
  acts <- pmax(acts, 0) #anything < 0 to 0
  return(acts / max(acts) ) #normalise 0..1
}

#finds the start and stop of the biggest peak 
find_biggest_peak_indices <- function(norm_acts, spar = 0.4){
  spsm <- stats::smooth.spline(norm_acts, spar =spar)
  left_of_peak <- right_of_peak <- NA
  valleys <- quantmod::findValleys(spsm$y)
  
  
  if (length(valleys) == 0){ # got no valleys
    return(list(
        start = NA, 
        end = NA,
        smooth_curve = spsm
      )
    )
  }
  mx_peak <- which.max(spsm$y)
  booleans <- valleys <= mx_peak
  if (sum( valleys > mx_peak) == 0){ #no valley ends greater than max peak then peak goes off end and right peak value == 116
    right_of_peak = 116
  } else {
    right_of_peak <- valleys[! booleans][1]
  }
  ##left of peak - no valley starts before max peak then left peak is 
  if (sum( valleys < mx_peak) == 0){
    left_of_peak = 1
  } else {
    left_of_peak <- valleys[booleans][length(valleys[booleans])]
  }
  
  if (left_of_peak == right_of_peak){ #found length one peak - not a peak
    return(list(
      start = NA, 
      end = NA,
      smooth_curve = spsm
    )
    )
  }

  #return(booleans)
  return(list(
             start = left_of_peak * 2, #account for half width of final conv layer 
             end = right_of_peak * 2,
             smooth_curve = spsm
             )
    )
}

#extracts an expression sub profile for a tf target pair. left = start of 
get_expression_subprofile <- function(exprs, tf, target,left, right, sub_only = TRUE){
  if (sub_only){
    warning(paste(tf, target, left, right))
    return( exprs[c(tf,target), left:right] )
  }
  else {
    e <- exprs[c(tf,target),]
    e[, 1:(left - 1)] <- 1
    e[, (right + 1):ncol(e)] <- 1
    return(e)
  }
}

load_expression_file_as_matrix <- function(file = "lib/normalised_data.csv"){
  dat <- suppressWarnings(as.matrix(readr::read_tsv(file, col_types = readr::cols() )))
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  ndat <- matrix(as.numeric(dat), ncol = ncol(dat))
  colnames(ndat) <- colnames(dat)
  rownames(ndat) <- rownames(dat)
  return(ndat)
}



find_informative_expression_regions <- function(input_profile,
                                                i,
                                                #tf = NULL, 
                                                #target = NULL,
                                                info = NULL,
                                                expr_data = NULL,
                                                model = NULL,
                                                layer = NULL
                                                ){
  tf <- info$TF[i]
  target <- info$Target[i]
  
 #input_profile <- array(input_profile, dim = c(1,2, 232, 1) )
  
  
  act_levels <- get_class_activation_values(input_profile, 
                                            model, 
                                            layer = "separable_conv2d_2")
  norm_acts <- normalise_activations(act_levels)
  peaks <- find_biggest_peak_indices(norm_acts) ## in this extract the smooth curve for debugging

  if( is.na(peaks$start) & is.na(peaks$end) ){ #found no peak
    return(list(informative_region = NA, peaks = peaks)
    )
  }

  exprs <- get_expression_subprofile(expr_data, 
                                     tf, 
                                     target, 
                                     peaks$start, peaks$end,
                                     sub_only = TRUE)
  return(list( 
               informative_region = exprs,
               peaks = peaks
              )
         )
}