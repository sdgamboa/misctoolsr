contigencyTable <- function(x, features, cols, signatures) {

    if (!any(class(x) %in% c("SummarizedExperiment", "TreeSummarizedExperiment"))) {
        stop("Input data must be of class (Tree)SummarizedExperiment", call. = FALSE)
    }

    row_data <- SummarizedExperiment::rowData(x) %>%
        tibble::as_tibble(rownames = features)

    # coerce columns into factors
    for (i in seq_along(cols)) {
        row_data[cols[i]] <- factor(row_data[[cols[i]]], levels = c(TRUE, FALSE), labels = c(cols[i], paste0("not_", cols[i])))
    }

    # create a new dataset with columns (factor) for each signature set
    counter <- 1
    signatures_tbl <- lapply(signatures, function(y) {
        fct_col <- factor(row_data[[features]] %in% y, levels = c(TRUE, FALSE), labels = c( names(signatures)[counter], paste0("not_", names(signatures)[counter])))
        counter <<- counter + 1
        fct_col
        }) %>%
        magrittr::set_names(paste0("contabXYZ.", names(signatures))) %>%
        dplyr::as_tibble()

    # combine datasets
    row_data <- dplyr::bind_cols(row_data, signatures_tbl)

    # create contingency tables
    vctr_len <- length(cols) * length((signatures))
    output <- vector("list", vctr_len)
    counter <- 1

    for (i in seq_along(cols)) {
        for (j in seq_along(names(signatures))){

            form <- paste0("~ ", cols[i], " + ", names(signatures_tbl)[j])
            contigency_table <- stats::xtabs(form, data = row_data)
            output[[counter]] <- contigency_table
            names(output)[counter] <- paste0(cols[i], "_x_" , names(signatures)[j])
            counter <- counter + 1
        }
    }

    return(output)
}

enrichmentTest <- function(x) {

    output <-  lapply(x, function(y) {
        p_value <- fisher.test(y, alternative = "g")$p.value
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

ora <- function(x, features, cols, signatures, alpha = 0.05) {
    contab <- contigencyTable(x = x, features = features, cols = cols, signatures = signatures)
    output <- enrichmentTest(contab) %>%
        dplyr::filter(p_value <= alpha)
    return(output)
}
