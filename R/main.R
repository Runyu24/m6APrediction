# Functions reused from Practical 5. Add documentation in your new package.
#' Encode DNA Sequences
#'
#' Converts DNA sequences into a data frame of categorical nucleotide
#' features, with one column for each nucleotide position.
#'
#' @param dna_strings A character vector of DNA sequences of equal length.
#'   Sequences should contain A, T, C, or G.
#' @return A data frame containing factor columns named nt_pos1, nt_pos2,
#'   and so on.
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Sites for Multiple Sequences
#'
#' Uses a fitted random forest model to predict m6A probabilities
#' and classification statuses for multiple input sequences.
#'
#' @param ml_fit A fitted random forest model.
#' @param feature_df A data frame containing gc_content, RNA_type,
#'   RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation and DNA_5mer.
#' @param positive_threshold A numeric classification threshold.
#'   Defaults to 0.5. Scores strictly greater than this value are
#'   classified as Positive.
#' @return A data frame containing the original input features and
#'   two additional columns: predicted_m6A_prob and
#'   predicted_m6A_status.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' model_path <- system.file("extdata", "rf_fit.rds",
#'                           package = "m6APrediction")
#' model <- readRDS(model_path)
#' data_path <- system.file("extdata", "m6A_input_example.csv",
#'                          package = "m6APrediction")
#' example_data <- read.csv(data_path)
#' prediction_multiple(model, example_data)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict an Individual m6A Site
#'
#' Uses a fitted random forest model to predict the m6A probability
#' and classification status for a single input sequence.
#'
#' @param ml_fit A fitted random forest model.
#' @param gc_content The GC content of the input sequence.
#' @param RNA_type The RNA type, such as mRNA or lincRNA.
#' @param RNA_region The RNA region, such as CDS or intron.
#' @param exon_length The log2-transformed exon length.
#' @param distance_to_junction The log2-transformed distance to the
#'   nearest junction.
#' @param evolutionary_conservation The evolutionary conservation score.
#' @param DNA_5mer A five-nucleotide DNA sequence.
#' @param positive_threshold A numeric classification threshold.
#'   Defaults to 0.5. Scores strictly greater than this value are
#'   classified as Positive.
#' @return A named character vector containing predicted_m6A_prob
#'   and predicted_m6A_status.
#' @export
#' @examples
#' model_path <- system.file("extdata", "rf_fit.rds",
#'                           package = "m6APrediction")
#' model <- readRDS(model_path)
#' data_path <- system.file("extdata", "m6A_input_example.csv",
#'                          package = "m6APrediction")
#' example_data <- read.csv(data_path)
#' prediction_single(
#'   ml_fit = model,
#'   gc_content = example_data$gc_content[1],
#'   RNA_type = as.character(example_data$RNA_type[1]),
#'   RNA_region = as.character(example_data$RNA_region[1]),
#'   exon_length = example_data$exon_length[1],
#'   distance_to_junction = example_data$distance_to_junction[1],
#'   evolutionary_conservation = example_data$evolutionary_conservation[1],
#'   DNA_5mer = as.character(example_data$DNA_5mer[1])
#' )
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
