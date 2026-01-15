#' Read and preprocess Ukraine study data
#'
#' This function reads the Ukraine baseline and endline data, cleans variable names,
#' processes metadata, transforms outcome labels, and prepares the data for analysis.
#' It returns participant-level data in long format along with item metadata.
#' Note: Mental health outcomes use different scales in Ukraine vs Colombia and won't map directly.
#'
#' @param file_data_Ukraine Character string. Path to the CSV data file containing Ukraine
#'   study data with baseline and endline measurements in wide format.
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
#' file_data_Ukraine <- "path/to/Ukraine_Hope_Groups_Baseline_Endline_Wide_Aug6.csv"
#' tmp <- read_data_ukraine(file_data_Ukraine)
#' dp <- copy(tmp$dp)
#' dit <- copy(tmp$dit)
#' dmeta <- copy(tmp$dmeta)
#' }
#'
#' @export
read_data_ukraine <- function(file_data_Ukraine) {
    require(data.table)
    # Suppress R CMD check notes about data.table's non-standard evaluation
    submission_date <- time_label <- pid <- n <- pid_new <- f_label <- fid <-
        time <- y <- y_label <- item_label <- item_type <- item_high_label <-
        group_label <- item_label_short <- group_label_long <- endpoint_measure <- NULL

    # Read in data from Ukraine
    dp <- as.data.table(read.csv(file_data_Ukraine))
    
    # Split wide format data into baseline (.x suffix) and endline (.y suffix)
    tmp <- subset(dp, select = c('UniqueID', colnames(dp)[grepl('.x$', colnames(dp))]))
    tmp2 <- subset(dp, select = c('UniqueID', colnames(dp)[grepl('.y$', colnames(dp))]))
    
    # Remove suffixes from column names
    setnames(tmp, old = colnames(tmp)[grepl('.x$', colnames(tmp))],
             new = gsub('.x$', '', colnames(tmp)[grepl('.x$', colnames(tmp))]))
    setnames(tmp2, old = colnames(tmp2)[grepl('.y$', colnames(tmp2))],
             new = gsub('.y$', '', colnames(tmp2)[grepl('.y$', colnames(tmp2))]))
    
    # Remove scale columns
    tmp <- subset(tmp, select = !grepl('^Scale_', colnames(tmp)))
    
    # Verify column consistency between baseline and endline
    stopifnot(setdiff(colnames(tmp), colnames(tmp2)) == character(0))
    stopifnot(setdiff(colnames(tmp2), colnames(tmp)) == character(0))
    
    # Combine baseline and endline data
    dp <- rbind(tmp, tmp2)
    
    # Standardize column names
    setnames(
        dp,
        c("UniqueID", "Timepoint", "f_name", "SubmissionDate"),
        c("pid", "time_label", "f_label", "submission_date")
    )

    # Separate out meta data/ covariates
    col_meta <- c("marital_status", "living_partner", "served_partnered",
        "demo_sex_labelled", "age_range_labelled", "edu_level_grped.labelled", "income_labelled", "country",
        "displacement_status", "shelter_now", "children", "health_disability2", "under_12months",
        "btwn_1_3_yrs", "btwn_4_7_yrs", "btwn_8_12yrs", "assistance.mhpss",
        "partner_sharing", "partner_conflict", "facilitator_relationship",
        "spillover_frequency", "spillover_book", "past_programs", "spiritual_strength", "resources_afterHG", 
        "life_worse_me", "life_worse_family", "shelter_now_clean", "training", "training_type")
    dmeta <- subset(dp, select = c("pid", "time_label", col_meta))

    # Remove non primary outcome data
    set(dp, NULL, col_meta, NULL)
    set(dp, NULL, c("Physical_Emotional_Violence7",
        "Positive_Parenting7", "Parental_Involvement7", "Parental_Monitoring7",
        "Resilience7", "Child_Wellbeing7", "IPV_Prevention7", "Key"), NULL)
    set(dp, NULL, c("CESD_hopeful", "Format_Final", 
        "IPV_Prevention", "overall_session_completion", "partner_conflict_num", 
        "PHQ4_add_ins1", "PHQ4_down_numeric", 
        "PHQ4_down_numeric_weight", "PHQ4_interest_numeric", 
        "PHQ4_interest_numeric_weight", "PHQ4_nervous_numeric", 
        "PHQ4_nervous_numeric_weight", "PHQ4_total", "PHQ4_worry_numeric", 
        "PHQ4_worry_numeric_weight", "report_attend_all_sessions", "sexual_viol_prevention", 
        "sexual_viol_prevention_num"), NULL)

    # Rename outcome labels - Violence
    setnames(
        dp,
        c("Physical_Emotional_Violence", "ICAST_PA_object", "ICAST_EA_scream"),
        c("CG-VIO_agg", "CG-VIO_ph-punish", "CG-VIO_scream")
    )

    # Rename up outcome labels - Mental Health
    # TODO: only the following last 4 items map to categorical 0 - 3
    # TODO: I don t know how they map to the Colombbia labels
    setnames(
         dp,
         c(
             "PHQ4_total_weight", "PHQ4_nervous", "PHQ4_interest", "PHQ4_worry", "PHQ4_down"
         ),
         c(
             "CG-MH_agg", "CG-MH_nervous", "CG-MH_effort", "CG-MH_hopeless", "CG-MH_sad"
         )
    )
    set(dp, NULL, c("PHQ4_anxious", "PHQ4_depress"), NULL)

    # Clean up outcome labels - Positive Parenting
    # TODO please check, I replaced here "PSSS_learn" with "APQ_PP_compliment", is this as it should be?
    setnames(
        dp,
        c("Positive_Parenting", "APQ_PP_compliment", "APQ_I_play"),
        c("CG-POS_agg", "CG-POS_praise", "CG-POS_play")
    )

    # Rename outcome labels - Parental Monitoring
    setnames(
        dp,
        c("Parental_Monitoring", "PPPS_accompained", "risk_rider"),
        c("CG-MONITOR-CHI_agg", "CG-MONITOR-CHI_safe-time", "CG-MONITOR-CHI_child-safe")
    )

    # Rename outcome labels - Parental Involvement
    setnames(
         dp,
         c("Parental_Involvement", "PSSS_learn", "share_problems"),
         c("CG-INVOLVE_agg", "CG-INVOLVE_help-learn", "CG-INVOLVE_child-problems")
    )

    # Clean up outcome labels - Child Behaviour    
    setnames(
         dp,
         c("Child_Wellbeing", "CABI_E_angry", "unhappy_internal", "CABI_I_interest"),
         c("CHI-BEHAVIOUR_agg", "CHI-BEHAVIOUR_angry", "CHI-BEHAVIOUR_unhappy", "CHI-BEHAVIOUR_no-interest")
    )

    # Clean up outcome labels - Self-care and Discipline
    setnames(
        dp,
        c("CESD_depressed", "selfcare", "grieve", "PARYC_SL_calmly"),
        c("CG-DEPRESSION", "CG-SELFCARE", "CG-RESILIENCE", "CG-NONVIOLENT-DISCIPLINE")
    )

    # TODO: "Resilience" does not map to 0-7 days as in Colombia
    set(dp, NULL, "Resilience", NULL)

    # Clean up date - handles ISO 8601 format automatically
    set(
        dp, NULL, "submission_date",
        dp[, as.Date(submission_date)]
    )

    # Set time id
    dp[, time := as.integer(time_label == "Endline")]

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

    # Recode mental health outcomes to 0-3 scale 
    for(x in c("CG-MH_nervous", "CG-MH_effort", "CG-MH_hopeless", "CG-MH_sad"))
    {
        set(dp, NULL, x, as.integer(gsub('([a-z]+)_([a-z]+)_([0-9])','\\3',(dp[[x]]))))
    }
    
    set(dp, NULL, 'CG-MH_agg', dp[, as.integer(`CG-MH_agg`)])

    # Bring table into long format
    dp <- data.table::melt(
        dp,
        id.vars = c(
            "time", "time_label", "pid", "pid_label", "fid",
            "f_label", "submission_date", "treat"
        ),
        variable.name = "item_label",
        value.name = "y"
    )

    # Remove NA's
    dp <- subset(dp, !is.na(y))

    # Define character values for y
    dp[, y_label := NA_character_]
    tmp <- dp[, which(grepl("^CG-MH_.*", item_label) & !grepl("agg", item_label))]
    stopifnot(length(tmp) > 0L)
    set(
        dp,
        tmp,
        "y_label",
        c(
            "a - not at all", "b - several days",
            "c - more than half of the time", "d - nearly every day"
        )[dp[tmp, y] + 1L]
    )
    tmp <- dp[, which(is.na(y_label) & !grepl("agg", item_label))]
    set(dp, tmp, "y_label", dp[tmp, paste0(as.character(y), " of 7 days")])

    # Show value mins and maxs
    dp[, list(ymin = min(y), ymax = max(y)), by = "item_label"]

    # Create item metadata table
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
