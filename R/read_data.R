#' Import and Process Raw SNP Data
#'
#' This function imports raw SNP data from a CSV file, checks for format consistency,
#' handles double-letter IUPAC notation, and computes basic summary statistics.
#'
#' @param file_path Path to the CSV file containing SNP data.
#'
#' @return An S3 object containing:
#' \describe{
#'  \item{snp_data}{A data frame with SNP data converted to single-letter IUPAC notation}
#'  \item{snp_summary}{A named integer vector of SNP call counts}
#'  \item{marker_status}{A data frame with the fraction of polymorphic, monomorphic, and failed markers}
#'  \item{allele_df}{A data frame with allele information for each marker}
#'  \item{marker_names}{A character vector containing the names of all markers}
#'  \item{output_dir}{The output directory to which data and figures are exported}
#' }
#' @export

read_raw_snp_data <- function(file_path, plot = FALSE) {
  # Function to import, check and convert raw SNP data into correct format
  raw_snp_df <- as.matrix(read.csv(file_path, sep = ",", quote = "", row.names = 1))

  # Extract the unique values
  unique_vals <- unique(raw_snp_df)

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
    processed_snp_data <- remap_snps(raw_snp_df, iupac_map)

    print("Double-letter IUPAC notation detected and converted.")
  } else if (double_letter == 0 & single_letter > 0) {
    # If only single-letter IUPAC notation is detected, return the raw data
    print("Single-letter IUPAC notation detected.")
    processed_snp_data <- raw_snp_df
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

  if (plot) {
    # Generate SNP distribution plot
    barplot_snp_file <- file.path(dirname(file_path), sub("\\.csv$", "_SNP_Distribution.png", basename(file_path)))
    plot_snp_distribution(snp_summary, output = barplot_snp_file, title = "SNP Call Distribution")
    # Generate Nei_H vs PIC boxplots
    boxplot_snp_file <- file.path(dirname(file_path), sub("\\.csv$", "_Nei_PIC_Boxplot.png", basename(file_path)))
    plot_nei_pic_boxplot(allele_df, output = boxplot_snp_file, title = "Nei's vs PIC")
    # Generate marker status barplots
    barplot_marker_file <- file.path(dirname(file_path), sub("\\.csv$", "_Marker_Status_Boxplot.png", basename(file_path)))
    plot_marker_stats(marker_stats, output = barplot_marker_file, title = "Marker Statistics")
    cat("Plots saved to:", dirname(file_path), "\n")
  }

  # Create the S3 object using the create_snphub_object function
  snp_object <- create_snphub_object(
    snp_data = processed_snp_data,
    snp_summary = snp_summary,
    marker_stats = marker_stats,
    allele_df = allele_df,
    marker_names = rownames(raw_snp_df),
    output_dir = dirname(file_path)
  )

  print(snphub_summary(snp_object))

  return(snp_object)
}
