#' Import and Process Raw SNP Data 
#'
#' This function imports raw SNP data from a CSV file, checks for format consistency,
#' handles double-letter IUPAC notation, and computes basic summary statistics.
#'
#' @param file_path Path to the CSV file containing SNP data. 
#'
#' @return A list containing:
#' \describe{
#'	\item{data}{A data frame with processed SNP data using single-letter IUPAC notation}
#'	\item{summary}{A list of SNP call counts}
#'}
#' @export

# Load helper functions and variables
#source("helpers.R")

read_raw_snp_data <- function(file_path, plot = TRUE) {
	
	# Function to import, check and convert raw SNP data into correct format
	
	raw_snp_df <- read.csv(file_path, sep = ",", quote = "")

	# Extract the unique values and exclude first row with marker names
	unique_vals <- unique(as.vector(as.matrix(raw_snp_df[,-1])))

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
				
		# Create a copy of the raw data for processing
		processed_snp_data <- raw_snp_df
		
		# Apply the IUPAC map to convert double-letter codes
		dat <- names(iupac_map)
		processed_snp_data[,-1] <- apply(processed_snp_data[,-1], c(1,2), function(x) {
			ifelse(x %in% dat, iupac_map[[x]], x)
		})
		
		print("Double-letter IUPAC notation detected and converted.")
			
	} else if (double_letter == 0 & single_letter > 0) {
		# If only single-letter IUPAC notation is detected, return the raw data
		print("Single-letter IUPAC notation detected.")
		processed_snp_data <- raw_snp_df
	}
	
	# Extract the unique values once again
	unique_vals <- unique(as.vector(as.matrix(processed_snp_data[,-1])))
	
	# Final check for unrecognized characters
	check_unrecognized(unique_vals, allowed_vals_single)
	
	# Compute basic statistics
	num_markers <- nrow(processed_snp_data) # Number of SNPs (rows)
	num_genotypes <- ncol(processed_snp_data) - 1 # Number of genotypes (excluding marker column)
	
	# Calculate SNP call statistics
	snp_summary <- count_snp_calls(as.matrix(processed_snp_data[,-1]))
	
	# Print summary information
	cat("\n--- SNP Data Summary ---\n")
	cat("Number of markers (SNPs):", num_markers, "\n")
	cat("Number of genotypes (samples):", num_genotypes, "\n")
	cat("\n--- SNP Call Counts ---\n")
	cat("A:", snp_summary[["A"]], "T:", snp_summary[["T"]], "C:", snp_summary[["C"]], "G:", snp_summary[["G"]], "\n")
	cat("R:", snp_summary[["R"]], "Y:", snp_summary[["Y"]], "S:", snp_summary[["S"]], "W:", snp_summary[["W"]], "K:", snp_summary[["K"]], "M:", snp_summary[["M"]], "\n")
	cat("Heterozygous SNPs:", snp_summary[["Hets"]], "\n")
	cat("Missing (Failed) Calls:", snp_summary[["-"]], "\n")
	cat("Total SNP Calls:", snp_summary[["Total"]], "\n")
	cat("------------------------\n")
	
	if (plot) {
    # Generate plot file path in the same directory
    plot_file <- file.path(dirname(file_path), sub("\\.csv$", "_snp_distribution.png", basename(file_path)))
    plot_snp_distribution(snp_summary, output = plot_file, title = "SNP Call Distribution")
    cat("Plot saved to:", plot_file, "\n")
	}
	
	return(list(data = processed_snp_data, summary = snp_summary))

}
