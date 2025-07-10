#' Create a SNPhub S3 Object
#'
#' @param raw_data Raw or processed SNP data frame.
#' @param summary SNP call summary statistics.
#' @param alleles Allele info (optional).
#' @param output_dir Output directory (optional).
#'
#' @return An S3 object of class "snp_data"
#' @export

create_snphub_object <- function(snp_data, snp_summary, marker_stats = NULL, output_dir = NULL,
                                 marker_names = NULL, allele_df = NULL) {
  obj <- list(
    snp_data = snp_data,
    snp_summary = snp_summary,
    marker_stats = marker_stats,
    allele_df = allele_df,
    marker_names = marker_names,
    output_dir = output_dir,
    filtered_snp_data = NULL
  )
  class(obj) <- "snp_data"
  return(obj)
}

#' Print Method for SNP Stats
#'
#' @param x The S3 objectcreated by read_raw_snp_data()
#'
#' @export

snphub_summary <- function(x, ...) {
  cat("\n--- Data Input ---\n")
  cat("Number of markers:", nrow(x$snp_data), "\n")
  cat("Number of genotypes (samples):", ncol(x$snp_data) - 1, "\n")
  cat("\n--- SNP Calls ---\n")
  cat("A:", x$snp_summary[["A"]], "T:", x$snp_summary[["T"]], "C:", x$snp_summary[["C"]], "G:", x$snp_summary[["G"]], "\n")
  cat("R:", x$snp_summary[["R"]], "Y:", x$snp_summary[["Y"]], "S:", x$snp_summary[["S"]], "W:", x$snp_summary[["W"]], "K:", x$snp_summary[["K"]], "M:", x$snp_summary[["M"]], "\n")
  cat("Heterozygous SNP Calls:", x$snp_summary[["Hets"]], "\n")
  cat("Missing (Failed) Calls:", x$snp_summary[["-"]], "\n")
  cat("Total SNP Calls:", x$snp_summary[["Total"]], "\n")
  cat("Median Heterozygosity (He):", median(x$allele_df$Nei_H, na.rm = TRUE), "\n")
  cat("Median Polymorphism Information Content (PIC):", median(x$allele_df$PIC, na.rm = TRUE), "\n")
  cat("\n--- After filtering ---\n")
  cat("Number of markers:", nrow(x$filtered_snp_data), "\n")
  cat("Number of genotypes (samples):", ncol(x$filtered_snp_data) - 1, "\n")
  cat("\n-----------------------\n")
  if (!is.null(x$output_dir)) {
    cat("Output saved to:", x$output_dir, "\n")
  }
}
