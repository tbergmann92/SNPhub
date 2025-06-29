# Helpers for SNPhub

#' Check for unrecognized SNP values
#'
#' @param unique_vals Vector of unique SNP values from the matrix
#' @param allowed_vals Vector of allowed characters

check_unrecognized <- function(unique_vals, allowed_vals) {
  unrecognized_chars <- unique_vals[!unique_vals %in% allowed_vals]
  if (length(unrecognized_chars) > 0) {
    print("Error: Unknown characters identified in the matrix!")
    print(paste("Identified characters:", paste(unrecognized_chars, collapse = ", ")))
    stop("Allowed characters: ", paste(allowed_vals, collapse = ", "))
  }
}

#' Calculate percentage
#'
#' @param count Numerical value
#' @param total Maximum value

calc_percent <- function(count, total) round((count / total) * 100, 2)

#' Convert string to individual characters
#'
#' @param str A string to be converted into chars
#' @return A vector of characters
#'
#' @examples
#' chars("ATGC")
#' # Returns: c("A", "T", "G", "C")
chars <- function(str) strsplit(ifelse(is.na(str), "", str), "")[[1]]

#' Function to assert a set of conditions
#' @param statement vector of comparison to proof
#' @param err_message string of error message when not all comparison are TRUE
#' @return character string used as error message or nothing if all TRUE
#' @keywords internal
assert <- function(statement, err_message = NULL) {
  if (!all(statement)) stop(err_message)
}

#' Function to replace SNPs values of an matrix by a named vector
#' @param mat original matrix
#' @param remap_list named vector of SNPs to replace the original matrix values with
#' @keywords internal
#' @examples
#' remap_snps(matrix(c("AA", "AG", "--", "GG"), ncol = 2), c("AA" = "A"))
#' # Returns: matrix(c("A", "AG", "--", "GG"), ncol = 2)
remap_snps <- function(mat, remap_list) {
  # First check if matrix values are in remap_list, if so, mark only those for replacement - replace
  in_remap_list <- mat %in% names(remap_list)
  mat[in_remap_list] <- as.character(remap_list[as.character(mat[in_remap_list])])
  return(mat)
}
