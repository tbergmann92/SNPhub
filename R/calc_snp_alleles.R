#' Analyse Allelic Status of SNPs 
#'
#' This function checks whether a SNP call is bi-allelic (polymorphic), monomorphic, or completely failed.
#'
#' @param processed_snp_data SNP data frame that has been converted to single-letter IUPAC by read_raw_snp_data(). 
#'
#' @return A list containing:
#' \describe{
#'	\item{}{A data frame the fraction of polymorphic, monomorphic, and failed SNP calls}
#'}
#' @export

# Load helper functions and variables
#source("helpers.R")

calc_snp_alleles <- function(processed_snp_data) {
	
	allele_a <- c()
	allele_b <- c()
	pol_counter <- 0
	mon_counter <- 0
	fail_counter <- 0
	
	# List to remove Hets and failed calls
	to_remove <- c("R", "Y", "S", "W", "K", "M", "H", "-")
	
	# Apply function across rows
	apply(processed_snp_data, 1, function(snps) {
		snps <- as.character(snps[-1]) # Skip first column (marker names)
		#print(str(snps))
		
		# Get unique alleles in order of appearance (order of appearance is important!)
		unique_snps <- unique(snps)
		#print(unique_snps)
		
		# Discard Hets and failed calls (logic is about to only count the "main" SNP calls
		main_snps <- unique_snps[!unique_snps %in% to_remove]
		#print(main_snps)
			
		if (length(main_snps) == 2) {
			allele_a <<- c(allele_a, main_snps[1])
			allele_b <<- c(allele_b, main_snps[2])
			pol_counter <<- pol_counter + 1
		} else if (length(main_snps) == 1) {
			allele_a <<- c(allele_a, main_snps[1])
			allele_b <<- c(allele_b, NA)
			mon_counter <<- mon_counter + 1
		} else {
			allele_a <<- c(allele_a, NA)
			allele_b <<- c(allele_b, NA)
			fail_counter <<- fail_counter + 1
		}	
	})
	
	total_markers <- nrow(processed_snp_data[,-1])
	
	# Move to helpers.R
	calc_percent <- function(count, total) round((count/total) * 100, 2)
	
	snp_call_stats <- data.frame(
		Type = c("Polymorphic", "Monomorphic", "Failed"),
		Fraction = c(
			calc_percent(pol_counter, total_markers),
			calc_percent(mon_counter, total_markers),
			calc_percent(fail_counter, total_markers)
			)
		)
		
	return(snp_call_stats)
	
	#return(list(
    #status = snp_call_stats
    #allele_a = allele_a,
    #allele_b = allele_b
	#))
	
}

