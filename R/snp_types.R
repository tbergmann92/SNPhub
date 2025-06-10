#' Allowed single-letter IUPAC codes
#'
#' A character vector containing all valid single-letter IUPAC nucleotide codes
#' @format A character vector with 11 elements
#' @examples
#' allowed_vals_single 
#' 
#' "A" = Adenine, "C" = Cytosine, "G" = Guanine, "T" = Thymine
#' "K" = G/T, "M" = A/C, "R" = A/G, "S" = G/C, "W" = A/T, "Y" = C/T
#' "-" = missing data
#'

allowed_vals_single <- c("-", "A", "C", "G", "T", "K", "M", "R", "S", "W", "Y")

#' Allowed double-letter IUPAC
#'
#' A character vector containing the valid double-letter IUPAC
#'
#' @format A character vector with 17 elements
#' @examples
#' allowed_vals_double 
#' 
#' "AA" = Adenine, "CC" = Cytosine, "GG" = Guanine, "TT" = Thymine
#' "K" = GT|TG, "M" = AC|CA, "R" = AG|GA, "S" = GC|CG, "W" = AT|TA, "Y" = CT|TC
#' "--" = missing data
#'

allowed_vals_double <- c("AA", "CC",
                         "GG", "TT",
                         "AC", "CA",
                         "AG", "GA",
                         "AT", "TA",
                         "CG", "GC",
                         "CT", "TC",
                         "GT", "TG",
                         "--")

#' IUPAC conversion table (double-letter to single-letter)
#'
#' An named character vector that maps valid double-letter IUPAC genotypes to
#' corresponding single-letter codes
#'
#' Used during genotype standardization in 'read_raw_snp_data()'
#'
#' @format Named character vector of length 17
#' @examples
#' iupac_map[["AG"]] 
#' # Returns "R"
#'
 
iupac_map <- c(
    "AA" = "A", "CC" = "C",
    "GG" = "G", "TT" = "T",
    "AC" = "M", "CA" = "M",
    "AG" = "R", "GA" = "R",
    "AT" = "W", "TA" = "W",
    "CG" = "S", "GC" = "S",
    "CT" = "Y", "TC" = "Y",
    "GT" = "K", "TG" = "K",
    "--" = "-"
)