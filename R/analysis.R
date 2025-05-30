# Analysis Functions for SNPhub

#'Calculate percentage
#'
#' @param count Numerical value
#' @param total Maximum value
#'
#' @export

calc_percent <- function(count, total) round((count/total) * 100, 2)

#' Count SNP Calls
#'
#' Takes a matrix of single-letter SNP calls and counts the frequency of each type.
#'
#' @param matrix_data A character matrix with SNP genotypes.
#' @return A named list with counts for each SNP type, heterozygous calls, and missing data.
#' @export

count_snp_calls <- function(matrix_data) {
	SNP <- "ATCGRYSWKM-"
	HETS <- "RYSWM"
	snp_counts <- sapply(chars(SNP), USE.NAMES = TRUE, function(x) sum(matrix_data == x, na.rm = TRUE))
	snp_counts[["Hets"]] <- sum(snp_counts[chars(HETS)])
	snp_counts[["Total"]] <- sum(snp_counts)
	return(snp_counts)
}

#' Analyse Allelic Status of SNPs 
#'
#' This function checks whether a SNP call is bi-allelic (polymorphic), monomorphic, or completely failed.
#' It summarizes allele counts for each SNP and calculates the frequencies.
#' It invokes the calc_polymorphism() function and generates the allele_df frame.
#' The allele_df data frame contains all essential information about each marker.
#'
#' @param processed_snp_data SNP data frame that has been converted to single-letter IUPAC by read_raw_snp_data(). 
#'
#' @return A list containing:
#' \describe{
#'	\item{}{A data frame the fraction of polymorphic, monomorphic, and failed SNP calls}
#'	\item{}{The allele_df data frame that contains all essential marker metrics}
#'}
#' @export

calc_allelic_stats <- function(processed_snp_data) {

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

#' Determining the Informativeness of a Marker 
#'
#' Calculates Nei's genetic diversity (expected heterozygosity, He) and the
#' Polymorphism Information Content (PIC) for each marker
#'
#' @param allele_df Data frame with allele count columns and frequencies for each marker 
#'
#' @return A modified version of the input data frame with three new columns:
#' \itemize{
#'	\item \code{Nei_H}:Nei's genetic diveristy (expected heterozgosity)
#'  \item \code{PIC}: Polymorphism Information Content
#'  \item \code{MAF}: Minor Allele Frequency (only for polymorphic markers)
#' }
#'
#' @export

calc_polymorphism <- function(df) {
	
	# Summarise the total calls
	total_calls <- df$Allele_A_Count + df$Allele_B_Count + df$Allele_AB_Count
	
	# Compute the corrected (total) allele frequencies using AB as half to A and to B
	p = (df$Allele_A_Count + 0.5 * df$Allele_AB_Count) / total_calls
	q = (df$Allele_B_Count + 0.5 * df$Allele_AB_Count) / total_calls
	
	# Nei's H: expected heterozygosity
	H <- round((1 - (p^2 + q^2)), digits=3)
	
	# PIC for biallelic markers (Botstein, 1980)
	PIC <- round(1 - (p^2 + q^2) - 2 * (p^2) * (q^2), digits = 3)
	
	# MAF: minor allele frequency
	MAF <- round(pmin(p, q), 3)
	
	# Assign and mask failed sites
	df$Nei_H <- H
	df$Nei_H[df$Failed_Freq == 1 | is.nan(H)] <- NA
	
	df$PIC <- PIC
	df$PIC[df$Failed_Freq == 1 | is.nan(PIC)] <- NA
	
	df$MAF <- ifelse(df$Status == "Pol", round(MAF, 3), NA)
	
	return(df)
}