#' Process Raw SNP Data
#'
#' This function checks for format consistency,
#' handles double-letter IUPAC notation, and computes basic summary statistics.
#'
#' @param snp_matrix Matrix with SNP data in specified format
#'
#' @return An S3 object containing:
#' \describe{
#'  \item{snp_data}{A data frame with SNP data converted to single-letter IUPAC notation}
#'  \item{snp_summary}{A named integer vector of SNP call counts}
#'  \item{marker_status}{A data frame with the fraction of polymorphic, monomorphic, and failed markers}
#'  \item{allele_df}{A data frame with allele information for each marker}
#'  \item{marker_names}{A character vector containing the names of all markers}
#' }
#' @export

snphub <- function(snp_matrix) {
  # TODO: implement check function for input snp matrix
  unique_vals <- unique(snp_matrix)

  # Check whether "failed" is used for "--/-" and stop immediately
  if (any(tolower(unique_vals) == "failed")) {
    message <- paste(
      "Error: 'failed' is used in the data instead of the allowed '--' or '-'.\n",
      "This issue may be caused by missing or improperly coded genotype information.\n",
      "Please check your input data file and ensure that '--' or '-' is used for missing data.\n",
      "Process cancelled due to invalid data format."
    )
    stop(message)
  }

  # Count single and double letter notations
  single_letter <- sum(nchar(unique_vals) == 1)
  double_letter <- sum(nchar(unique_vals) == 2)

  # Stop if double and single-letter IUPAC notation is detected
  if (double_letter > 0 & single_letter > 0) {
    message <- paste(
      "Error: Mixture of single- and double-letter IUPAC notation detected!\n",
      "Please check your input data file and ensure that SNP calls are either in single- or double-letter notation.\n",
      "Process cancelled due to invalid data format."
    )
    stop(message)
  }

  # If double-letter notation is detected, check characters and then convert it to single-letter notation
  else if (double_letter > 0 & single_letter == 0) {
    check_unrecognized(unique_vals, allowed_vals_double)

    # Apply the IUPAC map to convert double-letter codes
    processed_snp_data <- remap_snps(snp_matrix, iupac_map)

    print("Double-letter IUPAC notation detected and converted.")
  } else if (double_letter == 0 & single_letter > 0) {
    # If only single-letter IUPAC notation is detected, return the raw data
    print("Single-letter IUPAC notation detected.")
    processed_snp_data <- snp_matrix
  }

  # Final check for unrecognized characters
  # TODO: checking matrix twice should be avoided
  check_unrecognized(unique(processed_snp_data), allowed_vals_single)

  # Calculate SNP call statistics
  snp_summary <- count_snp_calls(processed_snp_data)

  # Calculate SNP allele statistics and get allele data frame
  allele_results <- calc_allelic_stats(processed_snp_data)
  marker_stats <- allele_results$marker_allele_stats
  allele_df <- allele_results$allele_df

  # Create the S3 object using the create_snphub_object function
  snp_object <- create_snphub_object(
    snp_data = processed_snp_data,
    snp_summary = snp_summary,
    marker_stats = marker_stats,
    allele_df = allele_df,
    marker_names = rownames(snp_matrix),
  )

  return(snp_object)
}
