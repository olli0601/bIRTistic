#' Read and preprocess Colombia study data
#'
#' This function reads the Colombia baseline and endline data, cleans variable names,
#' processes metadata, transforms outcome labels, and prepares the data for analysis.
#' It returns participant-level data in long format along with item metadata.
#'
#' @param file_data Character string. Path to the CSV data file containing Colombia
#'   study data with baseline and endline measurements.
#'
#' @return A list with three elements:
#'   \item{dp}{data.table in long format with participant outcomes by item and timepoint}
#'   \item{dit}{data.table with item metadata including item types and labels}
#'   \item{dmeta}{data.table with participant-level covariates and metadata}
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#' file_data <- "path/to/Colombia_data_baseline_endline_itemised.csv"
#' tmp <- read_data_colombia(file_data)
#' dp <- copy(tmp$dp)
#' dit <- copy(tmp$dit)
#' dmeta <- copy(tmp$dmeta)
#' }
#'
#' @export
read_data_colombia <- function(file_data) {
    require(data.table)
    # Suppress R CMD check notes about data.table's non-standard evaluation
    submission_date <- time_label <- pid <- n <- pid_new <- f_label <- fid <-
        time <- y <- y_label <- item_label <- item_type <- item_high_label <-
        group_label <- item_label_short <- group_label_long <- endpoint_measure <- NULL

    dp <- as.data.table(read.csv(file_data))

    setnames(
        dp,
        c("SID", "Timepoint", "staff_name"),
        c("pid", "time_label", "f_label")
    )

    # Separate out meta data/ covariates
    col_meta <- c(
        "age", "sex", "household_adults", "household_children", "child_impairment",
        "education", "moved", "maritalstat", "income", "income.adj",
        "income.per.person", "outfhelp", "ngp", "services", "services_FOOD",
        "services_HOUSING_SUBS", "services_CHILDCARE", "services_COUNSELING",
        "stressmeals", "ruv", "time_since_death",
        "Months.since.caregiver.death", "Months.since.caregiver.death.v2",
        "Months.since.caregiver.death.v3"
    )
    dmeta <- subset(dp, select = c("pid", "time_label", col_meta))

    # Keep core outcome data
    set(dp, NULL, col_meta, NULL)

    # Clean up outcome labels
    setnames(
        dp,
        c(
            "CAREGIVER_MENTAL_HEALTH", "nervous", "hopeless", "restless",
            "sad", "effort", "worthless"
        ),
        c(
            "CG-MH_agg", "CG-MH_nervous", "CG-MH_hopeless", "CG-MH_restless",
            "CG-MH_sad", "CG-MH_effort", "CG-MH_worthless"
        )
    )
    setnames(
        dp,
        c("PHYSICAL_EMOTIONAL_VIOLENCE", "physic_punish", "scream"),
        c("CG-VIO_agg", "CG-VIO_ph-punish", "CG-VIO_scream")
    )
    setnames(
        dp,
        c("POSITIVE_PARENTING", "praise", "play"),
        c("CG-POS_agg", "CG-POS_praise", "CG-POS_play")
    )
    setnames(
        dp,
        c("CHILD_MONITORING", "safe_time", "child_safe"),
        c(
            "CG-MONITOR-CHI_agg", "CG-MONITOR-CHI_safe-time",
            "CG-MONITOR-CHI_child-safe"
        )
    )
    setnames(
        dp,
        c("PARENTAL_INVOLVEMENT", "help_learn", "child_problems"),
        c("CG-INVOLVE_agg", "CG-INVOLVE_help-learn", "CG-INVOLVE_child-problems")
    )
    setnames(
        dp,
        c("CHILD_BEHAVIOURAL_ISSUES", "angry", "unhappy", "no_interest"),
        c(
            "CHI-BEHAVIOUR_agg", "CHI-BEHAVIOUR_angry", "CHI-BEHAVIOUR_unhappy",
            "CHI-BEHAVIOUR_no-interest"
        )
    )
    setnames(
        dp,
        c("DEPRESSION", "SELFCARE", "RESILIENCE", "NONVIOLENT_DISCIPLINE"),
        c("CG-DEPRESSION", "CG-SELFCARE", "CG-RESILIENCE", "CG-NONVIOLENT-DISCIPLINE")
    )

    # Clean up date
    set(
        dp, NULL, "submission_date",
        dp[, as.Date(submission_date, format = "%m/%d/%y")]
    )

    # Set time id
    dp[, time := as.integer(time_label == "Endline")]

    # Remove participant with double endline - remove last record
    dp <- subset(
        dp,
        !(pid %in% c("otmar20231963") & submission_date == "2024-12-12")
    )

    # Remove participants with only baseline records
    # as we are interested in pre-post comparison
    tmp <- dp[, list(n = length(submission_date)), by = "pid"]
    tmp <- subset(tmp, n == 2, select = pid)
    dp <- merge(tmp, dp, by = "pid")

    # Set participant id
    tmp <- data.table(pid = dp[, sort(unique(pid))])
    tmp[, pid_new := seq_len(nrow(tmp))]
    dp <- merge(dp, tmp, by = "pid")
    setnames(dp, c("pid", "pid_new"), c("pid_label", "pid"))

    # Set facilitator id
    tmp <- data.table(f_label = dp[, sort(unique(f_label))])
    tmp[, fid := seq_len(nrow(tmp))]
    dp <- merge(dp, tmp, by = "f_label")

    # Bring table into long format
    dp <- data.table::melt(
        dp,
        id.vars = c(
            "time", "time_label", "pid", "pid_label", "fid",
            "f_label", "submission_date", "d_year"
        ),
        variable.name = "item_label",
        value.name = "y"
    )

    # Remove NA's
    dp <- subset(dp, !is.na(y))

    # Define character values for y
    dp[, y_label := NA_character_]
    tmp <- dp[, which(grepl("^CG-MH_.*", item_label) & !grepl("agg", item_label))]
    set(
        dp,
        tmp,
        "y_label",
        c(
            "a - none of the time", "b - a little of the time",
            "c - some of the time", "d - most of the time",
            "e - all of the time"
        )[dp[tmp, y] + 1L]
    )
    tmp <- dp[, which(is.na(y_label) & !grepl("agg", item_label))]
    set(dp, tmp, "y_label", dp[tmp, paste0(as.character(y), " of 7 days")])

    # Show value mins and maxs
    dp[, list(ymin = min(y), ymax = max(y)), by = "item_label"]

    dit <- unique(subset(dp, select = "item_label"))
    dit[, item_type := ifelse(
        grepl("CG-MH", item_label), "categorical", "out-of-7"
    )]
    dit[, item_high_label := ifelse(
        grepl("CG-MH|CG-DEPRESSION|CG-VIO|CHI-BEHAVIOUR", item_label),
        "lower_is_better",
        "higher_is_better"
    )]
    dit[, group_label := factor(gsub("([^_]+)_([^_]+)", "\\1", item_label))]
    dit[, item_label_short := gsub("([^_]+)_([^_]+)", "\\2", item_label)]
    set(dit, dit[, which(grepl("^CG", item_label_short))], "item_label_short", "")
    dit[, group_label_long := group_label]
    tmp <- dit[, levels(group_label_long)]
    tmp <- sapply(tmp, function(x) {
        switch(x,
            "CG-MH" = "Caregiver mental health",
            "CG-VIO" = "Caregiver exercising physical or emotional violence",
            "CG-MONITOR-CHI_agg" = "Child monitoring",
            "CG-INVOLVE" = "Caregiver involvement",
            "CHI-BEHAVIOUR" = "Child behavioural issues",
            "CG-DEPRESSION" = "Caregiver depression",
            "CG-SELFCARE" = "Caregiver self-care",
            "CG-RESILIENCE" = "Caregiver resilience",
            "CG-POS" = "Caregiver positive parenting",
            "CG-MONITOR-CHI" = "Caregiver monitoring child",
            "CG-NONVIOLENT-DISCIPLINE" = "Caregiver exercising nonviolent discipline",
            x
        )
    })
    levels(dit$group_label_long) <- tmp
    dit[, endpoint_measure := sapply(item_type, function(x) {
        switch(x,
            "categorical" = "events occurring most or all of the time",
            "out-of-7" = "mean days in week",
            x
        )
    })]

    list(dp = dp, dit = dit, dmeta = dmeta)
}
