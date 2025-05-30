# Plot Functions for SNPhub

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
plot_marker_stats <- function(df, output, title) {

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
