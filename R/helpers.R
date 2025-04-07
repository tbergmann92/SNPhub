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
        las = 2 # Rotate x-axis labels for better readability
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

#' Count SNP Calls
#'
#' Takes a matrix of single-letter SNP calls and counts the frequency of each type.
#'
#' @param matrix_data A character matrix with SNP genotypes.
#' @return A named list with counts for each SNP type, heterozygous calls, and missing data.
#' @export
count_snp_calls <- function(matrix_data) {
	snp_counts <- list(
		A = sum(matrix_data == "A", na.rm = TRUE),
		T = sum(matrix_data == "T", na.rm = TRUE),
		C = sum(matrix_data == "C", na.rm = TRUE),
		G = sum(matrix_data == "G", na.rm = TRUE),
		R = sum(matrix_data == "R", na.rm = TRUE),
		Y = sum(matrix_data == "Y", na.rm = TRUE),
		S = sum(matrix_data == "S", na.rm = TRUE),
		W = sum(matrix_data == "W", na.rm = TRUE),
		K = sum(matrix_data == "K", na.rm = TRUE),
		M = sum(matrix_data == "M", na.rm = TRUE),
		Failed = sum(matrix_data == "-", na.rm = TRUE)
	)
		snp_counts$Hets <- sum(snp_counts$R, snp_counts$Y, snp_counts$S, 
							   snp_counts$W, snp_counts$K, snp_counts$M)
		snp_counts$Total <- sum(unlist(snp_counts))
		
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