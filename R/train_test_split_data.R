#' Split data into training and test sets
#'
#' @param dall Data table with item response data
#' @param test_individuals_p Proportion of individuals to include in test set
#' @param test_items_p Proportion of items to include in test set
#' @import data.table
#' @return A list with training and test data in from of two data.tables
train_test_split_data <- function(dall, test_individuals_p = 0.5, test_items_p = 0.5) {
    stopifnot(c("pid", "item_label", "oid") %in% colnames(dall))
    is_test_individual <- is_test_item <- oid <- NULL
    require(data.table)

    dtest <- data.table(pid = unique(dall$pid))
    dtest[, is_test_individual := rbinom(.N, 1, test_individuals_p)]
    dtest <- merge(dtest, unique(subset(dall, select = c("pid", "item_label"))), by = "pid")
    dtest[, is_test_item := rbinom(.N, 1, test_items_p)]
    dtest <- subset(dtest, is_test_individual == 1L & is_test_item == 1L, select = c("pid", "item_label"))
    dtest <- merge(dtest, dall, by = c("pid", "item_label"))

    setkey(dtest, oid)
    dtrain <- dall[!dtest, on = "oid"]

    return(list(train = dtrain, test = dtest))
}
