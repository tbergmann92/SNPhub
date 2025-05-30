# Helpers for SNPhub

#' Check for unrecognized SNP values
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

#' Convert string to individual characters
#'
#' @param str A string to be converted into chars
#' @return A vector of characters
#'
#' @export
#' @examples
#' chars("ATGC")
#' # Returns: c("A", "T", "G", "C")

chars <- function(str) strsplit(ifelse(is.na(str), "", str), "")[[1]]