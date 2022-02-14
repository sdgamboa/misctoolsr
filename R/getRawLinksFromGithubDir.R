#' Get raw links from GitHub directory
#'
#' /code{getRawLinksFromGithubDir} get the list of links leading to raw files
#' contained in a GitHub repository (recursively).
#'
#' @param repo Character string; user and repo names.
#' @param dir Character string; directory name.
#' @param ext Character string; the extension of the files. Default: ".*"
#' (every file in the directory).
#'
#' @return A character vector with the raw links of all the files.
#'
#' @export
#'
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#'
#' @examples
#' my_repo <- "bovee/fattyacids"
#' my_dir <- "reference_data"
#' my_ext <- "tsv"
#' raw_links <- getRawLinksFromGithubDir(
#'     repo = my_repo, dir = my_dir, ext = my_ext
#' )
getRawLinksFromGithubDir <- function(repo, dir, ext = ".*") {
    api_url <- 
        paste0("https://api.github.com/repos/", repo, "/git/trees/master")
    raw_url <- 
        paste0("https://raw.githubusercontent.com/", repo, "/master/")
    response <- httr::GET(api_url, query = list(recursive = 1))
    data <- jsonlite::fromJSON(rawToChar(response$content), flatten = TRUE)
    regex <- paste0(dir, "/.*", ext)
    paste0(raw_url, grep(regex, data$tree$path, value = TRUE))
}
