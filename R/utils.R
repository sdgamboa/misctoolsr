utils::globalVariables(c("expression1", "expression2"))
#' Example SummarizedExperment
#'
#' \code{\link{exampleSE}} creates an example of SummarizedExperiment (toy data).
#'
#' @return
#' A SummarizedExperiment.
#'
#' @importFrom stats rnorm
#' @importFrom magrittr set_rownames
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr mutate
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#'
#' library(misctoolsr)
#' se <- exampleSE()
#' se
#'
exampleSE <- function() {
    dimensions <- c(892, 76)
    m <- matrix(stats::rnorm(dimensions[1] * dimensions[2]),
                nrow = dimensions[1] , ncol = dimensions[2],
                dimnames = list(paste0("row", 1:dimensions[1]), paste0("col", 1:dimensions[2])))
    col_data <- data.frame(
        samples = colnames(m),
        condition = c(rep("control", dimensions[2] / 2), c(rep("treatment", dimensions[2] / 2)))
    ) %>%
        magrittr::set_rownames(colnames(m)) %>%
        S4Vectors::DataFrame()

    set.seed(1234)
    unchanged <- sample(rownames(m), 805)
    changed <- rownames(m)[!rownames(m) %in% unchanged]
    up <- sample(changed, 30)
    down <- changed[!changed %in% up]
    # rm(.Random.seed, envir=globalenv())

    row_data <- data.frame(
        taxa = rownames(m)
    ) %>%
        dplyr::mutate(
            expression1 = dplyr::case_when(taxa %in% unchanged ~ "not_differential", taxa %in% changed ~ "differential"),
            expression2 = dplyr::case_when(taxa %in% unchanged ~ "not_differential", taxa %in% up ~ "up", taxa %in% down ~ "down"),
            differential = ifelse(expression1 == "differential", TRUE, FALSE),
            up = ifelse(expression2 == "up", TRUE, FALSE),
            down = ifelse(expression2 == "down", TRUE, FALSE)
        ) %>%
        magrittr::set_rownames(rownames(m)) %>%
        S4Vectors::DataFrame()

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(assay1 = m),
        colData = col_data,
        rowData = row_data
    )
    return(se)
}

#' Example Signatures
#'
#' \code{\link{exampleSignatures}} creates an example of signatures (toy data).
#'
#' @return
#' A named list of signatures.
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @examples
#' library(misctoolsr)
#' sigs <- exampleSignatures()
#'
exampleSignatures <- function() {

    dimensions <- c(892, 76)
    m <- matrix(stats::rnorm(dimensions[1] * dimensions[2]),
                nrow = dimensions[1] , ncol = dimensions[2],
                dimnames = list(paste0("row", 1:dimensions[1]), paste0("col", 1:dimensions[2])))

    set.seed(1234)
    unchanged <- sample(rownames(m), 805)
    changed <- rownames(m)[!rownames(m) %in% unchanged]
    # set.seed(2930)
    up <- sample(changed, 30)
    down <- changed[!changed %in% up]
    # rm(.Random.seed, envir=globalenv())

    list_of_sigs <- list(
        sig1 = c(sample(up, 13, replace = FALSE),
                 sample(c(unchanged), 186, replace = FALSE),
                 sample(paste0("row", 893:1150), sample(1:100, 1), replace = FALSE)),
        sig2 = c(sample(up, 10, replace = FALSE),
                 sample(down, 38, replace = FALSE),
                 sample(unchanged, 277, replace = FALSE),
                 sample(paste0("row", 893:1150), sample(1:100, 1), replace = FALSE)),
        sig3 = c(sample(up, 7, replace = FALSE),
                 sample(down, 7, replace = FALSE),
                 sample(unchanged, 277, replace = FALSE),
                 sample(paste0("row", 893:1150), sample(1:100, 1), replace = FALSE))

        )

    return(list_of_sigs)

}
