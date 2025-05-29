# helper functions for SNPhub

#'Check for unrecognized SNP values
#'
#' @param unique_vals Vector of unique SNP values from the matrix
#' @param allowed_vals Vector of allowed characters
#'
#' @export
check_unrecognized <- function(unique_vals, allowed_vals) {
	unrecognized_chars <- unique_vals[!unique_vals %in% allowed_vals]
	if (length(unrecognized_chars) > 0) {
		print("Error: Unknown characters identified in the matrix!")
		print(paste("Identified characters:", paste(unrecognized_chars, collapse = ", ")))
		stop("Allowed characters: ", paste(allowed_vals, collapse = ", "))
	}
}

#'Calculate percentage
#'
#' @param count Numerical value
#' @param total Maximum value
#'
#' @export
calc_percent <- function(count, total) round((count/total) * 100, 2)

#'Plot SNP Call Distribution
#'
#' @param snp_counts A named list with SNP call counts
#' @param output Path to the output PNG file
#' @param title Title of the barplot
#'
#' @export
plot_snp_distribution <- function(snp_counts, output, title) {
	# Remove "Total" and optionally "Hets" from plotting
	plot_counts <- snp_counts[!names(snp_counts) %in% c("Total")]
	
	# Convert named list to a dataframe
    snp_df <- data.frame(
        Type = names(plot_counts),
        Fraction = unlist(plot_counts) / sum(unlist(plot_counts)) * 100  # Convert to percentage
    )

    # Set up plot colors and margins
    bar_colors <- "#6799a8"  # Bar fill color
    text_color <- "#120f0f"  # Text color

    # Save the plot as a PNG file
    png(output, width = 2500, height = 2500, res = 300)

    # Create the barplot
    bp <- barplot(
        snp_df$Fraction,
        names.arg = snp_df$Type,
        col = bar_colors,
        border = "black",
        ylim = c(0, 100),
        ylab = "Fraction [%]",
        main = title,
        las = 1 # Rotate x-axis labels for better readability
    )

    # Add text labels on top of bars
    text(
        x = bp, y = snp_df$Fraction + 2, 
        labels = round(snp_df$Fraction, 1), 
        col = text_color, font = 2, cex = 0.9
    )

    # Add grid lines for better readability
    grid(nx = NA, ny = NULL, col = "#9c9a9a", lty = "dotted")

    dev.off()  # Close the PNG device
}

#'Plot Polymorphism Information Content
#'
#' @param df The allele_df data frame from calc_allelic_status
#' @param output Path to the output PNG file
#' @param title Title of the boxplot
#'
#' @export
plot_nei_pic_boxplot <- function(df, output, title) {
  
  # Remove NA values
  df_clean <- df[!is.na(df$Nei_H) & !is.na(df$PIC),]
 
  # Prepare data
  nei_h <- df_clean$Nei_H
  pic   <- df_clean$PIC
  box_data <- list("Nei_H" = nei_h, "PIC" = pic)
  
  # Set up PNG output
  png(output, width = 1500, height = 2500, res = 300)
  
  # Set plot margins
  par(mar = c(6, 6, 4, 2))

  # Create boxplot
  boxplot(box_data,
		  names = c("Nei's H", "PIC"),
          main = title,
          ylab = "Value",
          col = "#6799a8",
          border = "black",
          cex.main = 1.5,
          cex.lab = 1.4,
          cex.axis = 1.2,
          las = 1)
		  
  # Overlay points (jittered horizontally)
  set.seed(1)
  points(jitter(rep(1, length(nei_h)), amount = 0.1), nei_h,
	pch = 16, col = "#120f0f50", cex = 0.8)
  points(jitter(rep(2, length(pic)), amount = 0.1), pic,
	pch = 16, col = "#120f0f50", cex = 0.8)
  
  # Add horizontal grid lines
  grid(nx = NA, ny = NULL, col = "#9c9a9a", lty = "dotted")

  dev.off()
}

#'Plot Marker Status
#'
#' @param df The data frame containing the fraction of polymorphic, monomorphic, and failed marker
#' @param output Path to the output PNG file
#' @param title Title of the barplot
#'
#' @export
plot_marker_status <- function(df, output, title) {

    # Set up plot colors and margins
    bar_colors <- "#6799a8"  # Bar fill color
    text_color <- "#120f0f"  # Text color

    # Save the plot as a PNG file
    png(output, width = 1500, height = 2500, res = 300)

    # Create the barplot
    bp <- barplot(
        df$Fraction,
        names.arg = df$Type,
        col = bar_colors,
        border = "black",
        ylim = c(0, 100),
        ylab = "Fraction [%]",
        main = title,
        las = 1 # Rotate x-axis labels for better readability
    )

    # Add text labels on top of bars
    text(
        x = bp, y = df$Fraction + 2, 
        labels = round(df$Fraction, 1), 
        col = text_color, font = 2, cex = 0.9
    )

    # Add grid lines for better readability
    grid(nx = NA, ny = NULL, col = "#9c9a9a", lty = "dotted")

    dev.off()  # Close the PNG device
}

#' Convert string to individual chars
#' @param str A string to be converted into chars
#' @return A vector of chars
#' @export
#' @examples
#' chars("ATGC")
#' # Returns: c("A", "T", "G", "C")
chars <- function(str) strsplit(str, "")[[1]]

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

# allowed_vals_single: Single letter IUPAC values 
allowed_vals_single <- c("-", "A", "C", "G", "T", "K", "M", "R", "S", "W", "Y")

# allowed_vals_double: Double letter IUPAC values
allowed_vals_double <- c("AA", "CC", "GG", "TT", 
						 "AC", "CA",
						 "AG", "GA",
						 "AT", "TA",
						 "CG", "GC",
						 "CT", "TC",
						 "GT", "TG",
						 "--")

# iupac_map: Mapping for double-letter IUPAC notation to single-letter notation 
iupac_map <- c(
	"AA" = "A", "CC" = "C", "GG" = "G", "TT" = "T",
	"AC" = "M", "CA" = "M",
	"AG" = "R", "GA" = "R",
	"AT" = "W", "TA" = "W",
	"CG" = "S", "GC" = "S",
	"CT" = "Y", "TC" = "Y",
	"GT" = "K", "TG" = "K",
	"--" = "-"
)

#' Determining the Informativeness of a Marker 
#'
#' This function calculates the heterozygosity value (He) of a marker after Nei (1978)
#' and the Polymorphism Information Content (PIC) of a marker after Botstein (1980). 
#'
#' @param allele_df Data frame that contains allele counts and frequencies for each marker. 
#'
#' @return Adds two additional columns to the allele data frame:
#' \describe{
#'	\item{}{The column contains the He and PIC value for each marker in the allele data frame}
#'}
#' @export

# Load helper functions and variables
#source("helpers.R")

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