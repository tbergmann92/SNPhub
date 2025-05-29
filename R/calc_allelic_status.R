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

calc_allelic_status <- function(processed_snp_data) {

	allele_a <- c()
	allele_b <- c()
	allele_ab <- c()
	allele_a_count <- c()
	allele_b_count <- c()
	allele_ab_count <- c()
	failed_count <- c()
	allele_a_freq <- c()
	allele_b_freq <- c()
	allele_ab_freq <- c()
	failed_freq <- c()
	pol_counter <- 0
	mon_counter <- 0
	fail_counter <- 0
	pic_values <- c()
	
	# List to remove Hets and failed calls
	to_remove <- chars("RYSWKMH-")
	iupac_het_codes <- chars("RYSWKMH")

	# Apply function across rows
	apply(processed_snp_data, 1, function(snps) {
		snps <- as.character(snps) 

		# Get unique alleles in order of appearance (order of appearance is important!)
		unique_snps <- unique(snps)

		# Discard Hets and failed calls (logic is about to only count the "main" SNP calls)
		main_snps <- unique_snps[!unique_snps %in% to_remove]

		# Identify the main alleles for each marker
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
		
		# Identify the heterozygous call if there is one
		het_call <- unique_snps[unique_snps %in% iupac_het_codes]
		allele_ab <<- c(allele_ab, if (length(het_call) > 0) paste(het_call) else NA)

		# Determine allele counts and frequencies
		allele_a_count <<- c(allele_a_count, sum(snps == main_snps[1], na.rm = TRUE))
		allele_b_count <<- c(allele_b_count, sum(snps == main_snps[2], na.rm = TRUE))
		allele_ab_count <<- c(allele_ab_count, sum(snps %in% iupac_het_codes, na.rm = TRUE))
		failed_count <<- c(failed_count, sum(snps == "-", na.rm = TRUE))
		allele_a_freq <<- c(allele_a_freq, round((sum(snps == main_snps[1], na.rm = TRUE) / length(snps)), digits=3))
		allele_b_freq <<- c(allele_b_freq, round((sum(snps == main_snps[2], na.rm = TRUE) / length(snps)), digits=3))
		allele_ab_freq <<- c(allele_ab_freq, round((sum(snps %in% iupac_het_codes, na.rm = TRUE) / length(snps)), digits=3))
		failed_freq <<- c(failed_freq, round((sum(snps == "-", na.rm = TRUE) / length(snps)), digits=3))

	})

	# Create a data frame with marker names and alleles
	allele_df <- data.frame(
		Marker = rownames(processed_snp_data),
		Allele_A = allele_a,
		Allele_B = allele_b,
		Allele_AB = allele_ab,
		Allele_A_Count = allele_a_count,
		Allele_B_Count = allele_b_count,
		Allele_AB_Count = allele_ab_count,
		Failed_Count = failed_count,
		Allele_A_Freq = allele_a_freq,
		Allele_B_Freq = allele_b_freq,
		Allele_AB_Freq = allele_ab_freq,
		Failed_Freq = failed_freq
		)

	# Calculate SNP allele stats
	total_markers <- nrow(processed_snp_data[,-1])
	marker_allele_stats <- data.frame(
		Type = c("Polymorphic", "Monomorphic", "Failed"),
		Fraction = c(
			calc_percent(pol_counter, total_markers),
			calc_percent(mon_counter, total_markers),
			calc_percent(fail_counter, total_markers)
			)
		)

	# Add whether a marker is polymorphic (Pol), monomorphic (Mon), or failed (Fail)
	allele_df$Status <- with(allele_df, ifelse(
		is.na(Allele_A) & is.na(Allele_B), "Fail",
		ifelse(!is.na(Allele_A) & is.na(Allele_B), "Mon", "Pol")
	))

	# 	Calculate heterozygosity (He) and PIC value for all marker
	allele_df <- calc_polymorphism(allele_df)
	
	
	return(list(marker_allele_stats = marker_allele_stats, 
				allele_df = allele_df
	))
	
}
