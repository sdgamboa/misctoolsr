utils::globalVariables(c("p_value", "odds_ratio", "upper_ci", "lower_ci", "n_sig_genes"))
#' Create contingency Table
#'
#' \code{contingencyTable} creates a list of 2x2 contingency tables for Fisher's exact test.
#'
#' @param x
#' A SummarizedExperiment with rowData or a data frame.
#' @param cols
#' A character vector containing the names of the columns (must be logical) with the target (TRUE) and reference sets (FALSE).
#' @param signatures
#' A named list with signature sets.
#' @param features
#' A character vector of length 1. In case that `x` is a SummarizedExperiment, this argument will be
#' used to create the features column. In case that `x` is a data frame, this argument must indicate
#' the name of the column containing the features names (taxa, genes, etc.).
#'
#' @return
#' A list of 2x2 contingenc tables (xtabs class).
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom tibble as_tibble
#' @importFrom magrittr set_names
#' @importFrom dplyr bind_cols
#' @importFrom stats xtabs
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @family enrichment functions
#' @seealso
#' \code{\link{contingencyTable}};
#' \code{\link{enrichmentTestFisher}};
#' \code{\link{ora}}
#'
#' @examples
#' library(misctoolsr)
#' se <- exampleSE()
#' sigs <- exampleSignatures()
#' contabs <- contingencyTable(x = se, cols = c("up", "down"), signatures = sigs)
#' contabs
#'
#' df <- as.data.frame(SummarizedExperiment::rowData(se))
#' contabs2 <- contingencyTable(x = df, features = "taxa", cols = c("up", "down"), signatures = sigs)
#'
#'
contingencyTable <- function(x, cols, signatures, features = "rownames") {

    if (any(class(x) %in% c("SummarizedExperiment", "TreeSummarizedExperiment"))) {
        df <- SummarizedExperiment::rowData(x) %>%
            tibble::as_tibble(rownames = features)
    } else if (any(class(x) %in% c("data.frame"))) {
        df <- x
    }

    # coerce columns into factors
    for (i in seq_along(cols)) {
        df[cols[i]] <- factor(df[[cols[i]]], levels = c(TRUE, FALSE), labels = c(cols[i], paste0("not_", cols[i])))
    }

    # create a new dataset with columns (factor) for each signature set
    counter <- 1
    signatures_tbl <- lapply(signatures, function(y) {
        fct_col <- factor(df[[features]] %in% y, levels = c(TRUE, FALSE), labels = c( names(signatures)[counter], paste0("not_", names(signatures)[counter])))
        counter <<- counter + 1
        fct_col
        }) %>%
        magrittr::set_names(paste0("contabXYZ.", names(signatures))) %>%
        tibble::as_tibble()

    # combine datasets
    df <- dplyr::bind_cols(df, signatures_tbl)

    # create contingency tables
    vctr_len <- length(cols) * length((signatures))
    output <- vector("list", vctr_len)
    counter <- 1

    for (i in seq_along(cols)) {
        for (j in seq_along(names(signatures))){

            form <- paste0("~ ", cols[i], " + ", names(signatures_tbl)[j])
            contigency_table <- stats::xtabs(form, data = df)
            output[[counter]] <- contigency_table
            names(output)[counter] <- paste0(cols[i], "_x_" , names(signatures)[j])
            counter <- counter + 1
        }
    }

    return(output)
}

#' Enrichment with Fisher's Exact Test (Hypergeometric)
#'
#' \code{enrichmentTestFisher} performs a one-sided overrepresentation enrichment analysis
#' using Fisher's exact test (hypergeometric). The function also calculates the odds.ratio
#' using the \code{epitools::oddsratio.wald} function.
#'
#' @param contabs
#' A iist of contingency tables (xtabs objects); output from \code{\link{contingencyTable}}.
#'
#' @return
#' A data frame; a report containing the p-value and odds ratio for each set x signature combination
#' in the input list  of contingency tables
#'
#' @importFrom stats fisher.test
#' @importFrom epitools oddsratio.wald
#' @importFrom tibble as_tibble
#' @importFrom magrittr set_names
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#'
#' @export
#' @family enrichment functions
#' @seealso
#' \code{\link{contingencyTable}};
#' \code{\link{enrichmentTestFisher}};
#' \code{\link{ora}}
#'
#' @examples
#' library(misctoolsr)
#' se <- exampleSE()
#' sigs <- exampleSignatures()
#' contabs <- contingencyTable(x = se, cols = c("up", "down"), signatures = sigs)
#' contabs
#'
#' enrichment_result <- enrichmentTestFisher(contabs)
#' enrichment_result
#'
enrichmentTestFisher <- function(contabs) {

    output <-  lapply(contabs, function(y) {
        p_value <- stats::fisher.test(y, alternative = "g")$p.value
        odds_ratio <- suppressWarnings(epitools::oddsratio.wald(y + 0.5)$measure[2,1])
        upper_ci <- exp(log(odds_ratio) + 1.96 * sqrt(sum(1 / (y + 1) )))
        lower_ci <- exp(log(odds_ratio) - 1.96 * sqrt(sum(1 / (y + 1) )))
        n_sig_genes <- y[1,1]
        c(p_value, odds_ratio, upper_ci, lower_ci, n_sig_genes)
    }) %>%
        as.data.frame() %>%
        t() %>%
        tibble::as_tibble(rownames = "rownames") %>%
        magrittr::set_colnames(c("set_x_signature",  "p_value", "odds_ratio", "upper_ci", "lower_ci", "n_sig_genes")) %>%
        dplyr::arrange(p_value, odds_ratio)

    return(output)

}

#' Overrepresentation Analysis
#'
#' \code{ora} performs an overrepresentation enrichment analysis (Fisher's test; hypergeometric)
#' given a data set and a list of signatures. This function is a wrapper around
#' \code{\link{contingencyTable}} and \code{\link{enrichmentTestFisher}}.
#'
#' @param x
#' A SummarizedExperiment with rowData or a data frame.
#' @param cols
#' A character vector containing the names of the columns (must be logical) with the target (TRUE) and reference sets (FALSE).
#' @param signatures
#' A named list with signature sets.
#' @param features
#' A character vector of length 1. In case that `x` is a SummarizedExperiment, this argument will be
#' used to create the features column. In case that `x` is a data frame, this argument must indicate
#' the name of the column containing the features names (taxa, genes, etc.).
#' @param alpha
#' Minimum p-value to filter results. Defult is 0.05.
#'
#' @return
#' A data frame; a report containing the p-value and odds ratio for each set x signature combination
#' in the input list  of contingency tables
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @family enrichment functions
#' @seealso
#' \code{\link{contingencyTable}};
#' \code{\link{enrichmentTestFisher}};
#' \code{\link{ora}}
#'
#' @examples
#' library(misctoolsr)
#' se <- exampleSE()
#' sigs <- exampleSignatures()
#' ora_results <- ora(x = se, cols = c("up", "down"), signatures = sigs)
ora <- function(x, cols, signatures, features = "rownames", alpha = 0.05) {
    contab <- contingencyTable(x = x, features = features, cols = cols, signatures = signatures)
    output <- enrichmentTestFisher(contab) %>%
        dplyr::filter(p_value <= alpha)
    return(output)
}
