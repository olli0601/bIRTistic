"""
Data loading functions for bIRTistic.

This module provides functions to load and preprocess IRT survey data,
ported from the original R implementations.
"""

from typing import Dict
import re
import pandas as pd
import numpy as np
from pathlib import Path


def read_data_colombia(file_data: str) -> Dict[str, pd.DataFrame]:
    """
    Read and preprocess Colombia study data.
    
    This function reads the Colombia baseline and endline data, cleans variable names,
    processes metadata, transforms outcome labels, and prepares the data for analysis.
    It returns participant-level data in long format along with item metadata.
    
    Parameters
    ----------
    file_data : str
        Path to the CSV data file containing Colombia study data
        with baseline and endline measurements.
    
    Returns
    -------
    dict
        Dictionary with three DataFrames:
        - 'dp': Participant data in long format (outcomes by item and timepoint)
        - 'dit': Item metadata (item types, labels, categories)
        - 'dmeta': Participant metadata (covariates, demographics)
    
    Examples
    --------
    >>> data = read_data_colombia("path/to/Colombia_data.csv")
    >>> dp = data['dp']
    >>> dit = data['dit']
    >>> dmeta = data['dmeta']
    """
    # Validate file exists
    if not Path(file_data).exists():
        raise FileNotFoundError(f"Data file not found: {file_data}")
    
    # Read CSV
    dp = pd.read_csv(file_data)
    
    # convert spaces to dots in column names
    dp.columns = [col.replace(' ', '.') for col in dp.columns]
    
    # Rename columns
    dp = dp.rename(columns={
        'SID': 'pid',
        'Timepoint': 'time_label',
        'staff_name': 'f_label'
    })
    
    # Separate out metadata/covariates
    col_meta = [
        'age', 'sex', 'household_adults', 'household_children', 'child_impairment',
        'education', 'moved', 'maritalstat', 'income', 'income.adj',
        'income.per.person', 'outfhelp', 'ngp', 'services', 'services_FOOD',
        'services_HOUSING_SUBS', 'services_CHILDCARE', 'services_COUNSELING',
        'stressmeals', 'ruv', 'time_since_death',
        'Months.since.caregiver.death', 'Months.since.caregiver.death.v2',
        'Months.since.caregiver.death.v3'
    ]
    dmeta = dp[['pid', 'time_label'] + col_meta].copy()
    dmeta = dmeta.rename(columns={'pid': 'pid_label'})
    
    # Keep core outcome data (drop metadata columns)
    dp = dp.drop(columns=col_meta)
    
    # Clean up outcome labels
    dp = dp.rename(columns={
        'CAREGIVER_MENTAL_HEALTH': 'CG-MH_agg',
        'nervous': 'CG-MH_nervous',
        'hopeless': 'CG-MH_hopeless',
        'restless': 'CG-MH_restless',
        'sad': 'CG-MH_sad',
        'effort': 'CG-MH_effort',
        'worthless': 'CG-MH_worthless',
        'PHYSICAL_EMOTIONAL_VIOLENCE': 'CG-VIO_agg',
        'physic_punish': 'CG-VIO_ph-punish',
        'scream': 'CG-VIO_scream',
        'POSITIVE_PARENTING': 'CG-POS_agg',
        'praise': 'CG-POS_praise',
        'play': 'CG-POS_play',
        'CHILD_MONITORING': 'CG-MONITOR-CHI_agg',
        'safe_time': 'CG-MONITOR-CHI_safe-time',
        'child_safe': 'CG-MONITOR-CHI_child-safe',
        'PARENTAL_INVOLVEMENT': 'CG-INVOLVE_agg',
        'help_learn': 'CG-INVOLVE_help-learn',
        'child_problems': 'CG-INVOLVE_child-problems',
        'CHILD_BEHAVIOURAL_ISSUES': 'CHI-BEHAVIOUR_agg',
        'angry': 'CHI-BEHAVIOUR_angry',
        'unhappy': 'CHI-BEHAVIOUR_unhappy',
        'no_interest': 'CHI-BEHAVIOUR_no-interest',
        'DEPRESSION': 'CG-DEPRESSION',
        'SELFCARE': 'CG-SELFCARE',
        'RESILIENCE': 'CG-RESILIENCE',
        'NONVIOLENT_DISCIPLINE': 'CG-NONVIOLENT-DISCIPLINE'
    })
    
    # Clean up date
    dp['submission_date'] = pd.to_datetime(dp['submission_date'], format='%m/%d/%y')
    
    # Set time id
    dp['time'] = (dp['time_label'] == 'Endline').astype(int)
    
    # Remove participant with double endline - remove last record
    dp = dp[~((dp['pid'] == 'otmar20231963') & (dp['submission_date'] == '2024-12-12'))]
    
    # Select participants who have both baseline and endline records
    participant_counts = dp.groupby('pid')['submission_date'].count().reset_index(name='n')
    participants_with_both = participant_counts[participant_counts['n'] == 2]['pid']
    dp = dp[dp['pid'].isin(participants_with_both)]
    
    # Set participant id (sequential)
    pid_mapping = pd.DataFrame({
        'pid': sorted(dp['pid'].unique())
    })
    pid_mapping['pid_new'] = range(1, len(pid_mapping) + 1)
    dp = dp.merge(pid_mapping, on='pid')
    dp = dp.rename(columns={'pid': 'pid_label', 'pid_new': 'pid'})
    
    # Set facilitator id (sequential)
    f_mapping = pd.DataFrame({
        'f_label': sorted(dp['f_label'].unique())
    })
    f_mapping['fid'] = range(1, len(f_mapping) + 1)
    dp = dp.merge(f_mapping, on='f_label')
    
    # Convert aggregated columns to integer (truncate like R's as.integer())
    agg_cols = ['CG-INVOLVE_agg', 'CHI-BEHAVIOUR_agg', 'CG-MONITOR-CHI_agg',
                'CG-MH_agg', 'CG-VIO_agg', 'CG-POS_agg']
    for col in agg_cols:
        if col in dp.columns:
            dp[col] = np.floor(dp[col]).astype('Int64')
    
    # Bring table into long format
    id_vars = ['time', 'time_label', 'pid', 'pid_label', 'fid', 
               'f_label', 'submission_date', 'd_year']
    dp = pd.melt(
        dp,
        id_vars=id_vars,
        var_name='item_label',
        value_name='y'
    )
    
    # Remove NA's
    dp = dp.dropna(subset=['y'])
    
    # Define character values for y
    dp['y_label'] = pd.Series(dtype='object')  # Initialize as object dtype
    
    # For CG-MH items (not aggregates)
    cg_mh_mask = (
        dp['item_label'].str.contains('^CG-MH_', regex=True) & 
        ~dp['item_label'].str.contains('agg')
    )
    y_labels = [
        'a - none of the time',
        'b - a little of the time',
        'c - some of the time',
        'd - most of the time',
        'e - all of the time'
    ]
    dp.loc[cg_mh_mask, 'y_label'] = dp.loc[cg_mh_mask, 'y'].apply(
        lambda y: y_labels[int(y)] if pd.notna(y) and 0 <= int(y) < len(y_labels) else np.nan
    )
    
    # For other items (not aggregates)
    other_mask = dp['y_label'].isna() & ~dp['item_label'].str.contains('agg')
    dp.loc[other_mask, 'y_label'] = dp.loc[other_mask, 'y'].apply(
        lambda y: f"{int(y)} of 7 days" if pd.notna(y) else np.nan
    )
    
    # Create item metadata table
    dit = dp[['item_label']].drop_duplicates().sort_values('item_label').reset_index(drop=True)
    
    # Set item_type
    dit['item_type'] = np.where(
        dit['item_label'].str.contains('CG-MH'),
        'categorical',
        'out-of-7'
    )
    
    # Set item_high_label
    dit['item_high_label'] = np.where(
        dit['item_label'].str.contains('CG-MH|CG-DEPRESSION|CG-VIO|CHI-BEHAVIOUR', regex=True),
        'lower_is_better',
        'higher_is_better'
    )
    
    # Extract group_label from item_label
    dit['group_label'] = dit['item_label'].str.replace(r'([^_]+)_([^_]+)', r'\1', regex=True)
    
    # Extract item_label_short
    dit['item_label_short'] = dit['item_label'].str.replace(r'([^_]+)_([^_]+)', r'\2', regex=True)
    # Set to NaN if starts with 'CG'
    dit.loc[dit['item_label_short'].str.startswith('CG'), 'item_label_short'] = np.nan
    
    # Create group_label_long with full names
    group_mapping = {
        'CG-MH': 'Caregiver mental health',
        'CG-VIO': 'Caregiver exercising physical or emotional violence',
        'CG-MONITOR-CHI_agg': 'Child monitoring',
        'CG-INVOLVE': 'Caregiver involvement',
        'CHI-BEHAVIOUR': 'Child behavioural issues',
        'CG-DEPRESSION': 'Caregiver depression',
        'CG-SELFCARE': 'Caregiver self-care',
        'CG-RESILIENCE': 'Caregiver resilience',
        'CG-POS': 'Caregiver positive parenting',
        'CG-MONITOR-CHI': 'Caregiver monitoring child',
        'CG-NONVIOLENT-DISCIPLINE': 'Caregiver exercising nonviolent discipline'
    }
    dit['group_label_long'] = dit['group_label'].map(group_mapping).fillna(dit['group_label'])
    
    # Set endpoint_measure
    dit['endpoint_measure'] = dit['item_type'].map({
        'categorical': 'events occurring most or all of the time',
        'out-of-7': 'mean days in week'
    })
    
    # Set cat_length
    dit['cat_length'] = dit['item_type'].map({
        'categorical': 5,
        'out-of-7': 8
    }).astype('Int64')
    
    # Set item_type_id
    type_mapping = pd.DataFrame({
        'item_type': sorted(dit['item_type'].unique())
    })
    type_mapping['item_type_id'] = range(1, len(type_mapping) + 1)
    dit = dit.merge(type_mapping, on='item_type')
    
    # Reset indices
    dp = dp.reset_index(drop=True)
    dit = dit.reset_index(drop=True)
    dmeta = dmeta.reset_index(drop=True)
    
    return {'dp': dp, 'dit': dit, 'dmeta': dmeta}


def read_data_ukraine(file_data: str) -> Dict[str, pd.DataFrame]:
    """
    Read and preprocess Ukraine study data.
    
    This function reads the Ukraine baseline and endline data from wide format,
    cleans variable names, processes metadata, transforms outcome labels, and
    prepares the data for analysis. It returns participant-level data in long
    format along with item metadata.
    
    Note: Mental health outcomes use different scales in Ukraine vs Colombia
    and won't map directly. Ukraine uses PHQ-4 (0-3 scale), Colombia uses
    different instruments (0-4 scale).
    
    Parameters
    ----------
    file_data : str
        Path to the CSV data file containing Ukraine study data
        with baseline and endline measurements in wide format.
        Baseline columns have '.x' suffix, endline have '.y' suffix.
    
    Returns
    -------
    dict
        Dictionary with three DataFrames:
        - 'dp': Participant data in long format (outcomes by item and timepoint)
        - 'dit': Item metadata (item types, labels, categories)
        - 'dmeta': Participant metadata (covariates, demographics)
    
    Examples
    --------
    >>> data = read_data_ukraine("path/to/Ukraine_data.csv")
    >>> dp = data['dp']
    >>> dit = data['dit']
    >>> dmeta = data['dmeta']
    """
    # Validate file exists
    if not Path(file_data).exists():
        raise FileNotFoundError(f"Data file not found: {file_data}")
    
    # Read CSV
    dp = pd.read_csv(file_data)
    
    # R's read.csv() converts spaces and special characters to dots in column names
    # Replicate R's make.names() behavior: spaces, /, and other special chars -> dots
    dp.columns = [col.replace(' ', '.').replace('/', '.') for col in dp.columns]
    
    # Split wide format data into baseline (.x suffix) and endline (.y suffix)
    baseline_cols = [col for col in dp.columns if col.endswith('.x')] + ['UniqueID']
    endline_cols = [col for col in dp.columns if col.endswith('.y')] + ['UniqueID']
    
    baseline = dp[baseline_cols].copy()
    endline = dp[endline_cols].copy()
    
    # Remove suffixes from column names
    baseline.columns = [col.replace('.x', '') for col in baseline.columns]
    endline.columns = [col.replace('.y', '') for col in endline.columns]
    
    # Remove Scale_ columns from baseline
    baseline = baseline[[col for col in baseline.columns if not col.startswith('Scale_')]]
    
    # Verify column consistency between baseline and endline
    assert set(baseline.columns) == set(endline.columns), \
        "Column mismatch between baseline and endline data"
    
    # Combine baseline and endline
    dp = pd.concat([baseline, endline], ignore_index=True)
    
    # Standardize column names
    dp = dp.rename(columns={
        'UniqueID': 'pid',
        'Timepoint': 'time_label',
        'f_name': 'f_label',
        'SubmissionDate': 'submission_date'
    })
    
    # Separate out metadata/covariates
    col_meta = [
        "marital_status", "living_partner", "served_partnered",
        "demo_sex_labelled", "age_range_labelled", "edu_level_grped.labelled", 
        "income_labelled", "country", "displacement_status", "shelter_now", 
        "children", "health_disability2", "under_12months", "btwn_1_3_yrs", 
        "btwn_4_7_yrs", "btwn_8_12yrs", "assistance.mhpss", "partner_sharing", 
        "partner_conflict", "facilitator_relationship", "spillover_frequency", 
        "spillover_book", "past_programs", "spiritual_strength", 
        "resources_afterHG", "life_worse_me", "life_worse_family", 
        "shelter_now_clean", "training", "training_type"
    ]
    dmeta = dp[['pid', 'time_label'] + col_meta].copy()
    
    # Remove non-primary outcome data
    dp = dp.drop(columns=col_meta)
    
    # Remove additional columns
    cols_to_drop = [
        "Physical_Emotional_Violence7", "Positive_Parenting7", "Parental_Involvement7",
        "Parental_Monitoring7", "Resilience7", "Child_Wellbeing7", "IPV_Prevention7", "Key",
        "Format_Final", "Resilience", "IPV_Prevention", "overall_session_completion",
        "partner_conflict_num", "PHQ4_add_ins1", "PHQ4_down_numeric",
        "PHQ4_down_numeric_weight", "PHQ4_interest_numeric", "PHQ4_interest_numeric_weight",
        "PHQ4_nervous_numeric", "PHQ4_nervous_numeric_weight", "PHQ4_total",
        "PHQ4_worry_numeric", "PHQ4_worry_numeric_weight", "report_attend_all_sessions",
        "sexual_viol_prevention", "sexual_viol_prevention_num"
    ]
    dp = dp.drop(columns=[col for col in cols_to_drop if col in dp.columns])
    
    # Rename outcome labels - Violence
    dp = dp.rename(columns={
        "Physical_Emotional_Violence": "CG-VIO_agg",
        "ICAST_PA_object": "CG-VIO_ph-punish",
        "ICAST_EA_scream": "CG-VIO_scream"
    })
    
    # Rename outcome labels - Mental Health
    # Note: MH variables do not map to Colombia labels - different survey scales
    dp = dp.rename(columns={
        "PHQ4_total_weight": "CG-MH_agg",
        "PHQ4_nervous": "CG-MH_nervous",
        "PHQ4_interest": "CG-MH_effort",
        "PHQ4_worry": "CG-MH_hopeless",
        "PHQ4_down": "CG-MH_sad"
    })
    
    # Remove unused MH columns
    dp = dp.drop(columns=[col for col in ["PHQ4_anxious", "PHQ4_depress"] if col in dp.columns])
    
    # Rename outcome labels - Positive Parenting
    dp = dp.rename(columns={
        "Positive_Parenting": "CG-POS_agg",
        "APQ_PP_compliment": "CG-POS_praise",
        "APQ_I_play": "CG-POS_play"
    })
    
    # Rename outcome labels - Parental Monitoring
    dp = dp.rename(columns={
        "Parental_Monitoring": "CG-MONITOR-CHI_agg",
        "PPPS_accompained": "CG-MONITOR-CHI_safe-time",
        "risk_rider": "CG-MONITOR-CHI_child-safe"
    })
    
    # Rename outcome labels - Parental Involvement
    dp = dp.rename(columns={
        "Parental_Involvement": "CG-INVOLVE_agg",
        "PSSS_learn": "CG-INVOLVE_help-learn",
        "share_problems": "CG-INVOLVE_child-problems"
    })
    
    # Rename outcome labels - Child Behaviour
    dp = dp.rename(columns={
        "Child_Wellbeing": "CHI-BEHAVIOUR_agg",
        "CABI_E_angry": "CHI-BEHAVIOUR_angry",
        "unhappy_internal": "CHI-BEHAVIOUR_unhappy",
        "CABI_I_interest": "CHI-BEHAVIOUR_no-interest"
    })
    
    # Rename outcome labels - Self-care and Discipline
    dp = dp.rename(columns={
        "CESD_depressed": "CG-DEPRESSION",
        "selfcare": "CG-SELFCARE",
        "CESD_hopeful": "CG-RESILIENCE",
        "PARYC_SL_calmly": "CG-NONVIOLENT-DISCIPLINE"
    })
    
    # Clean up date - pd.to_datetime handles ISO 8601 format automatically
    dp['submission_date'] = pd.to_datetime(dp['submission_date'])
    
    # Set time id
    dp['time'] = (dp['time_label'] == 'Endline').astype(int)
    
    # Remove participants with only baseline records (we need pre-post comparison)
    participant_counts = dp.groupby('pid')['submission_date'].count().reset_index(name='n')
    participants_with_both = participant_counts[participant_counts['n'] == 2]['pid']
    dp = dp[dp['pid'].isin(participants_with_both)].copy()
    
    # Set participant id (sequential)
    pid_mapping = pd.DataFrame({'pid': sorted(dp['pid'].unique())})
    pid_mapping['pid_new'] = range(1, len(pid_mapping) + 1)
    dp = dp.merge(pid_mapping, on='pid')
    dp = dp.rename(columns={'pid': 'pid_label', 'pid_new': 'pid'})
    
    # Set facilitator id (sequential)
    f_mapping = pd.DataFrame({'f_label': sorted(dp['f_label'].unique())})
    f_mapping['fid'] = range(1, len(f_mapping) + 1)
    dp = dp.merge(f_mapping, on='f_label')
    
    # Recode mental health outcomes to 0-3 scale
    # These come as strings like "phq4_nervous_0", "phq4_nervous_1", etc.
    mh_cols = ["CG-MH_nervous", "CG-MH_effort", "CG-MH_hopeless", "CG-MH_sad"]
    for col in mh_cols:
        if col in dp.columns:
            # Extract trailing digit from strings like "phq4_nervous_2" -> "2"
            dp[col] = dp[col].astype(str).str.extract(r'([0-9])$')[0].astype('Int64')
    
    # Convert aggregate to integer (truncate, don't round)
    if 'CG-MH_agg' in dp.columns:
        dp['CG-MH_agg'] = np.floor(dp['CG-MH_agg']).astype('Int64')
    
    # Bring table into long format
    id_vars = ['time', 'time_label', 'pid', 'pid_label', 'fid', 
               'f_label', 'submission_date', 'treat']
    dp = pd.melt(
        dp,
        id_vars=id_vars,
        var_name='item_label',
        value_name='y'
    )
    
    # Remove NA's
    dp = dp.dropna(subset=['y'])
    
    # Define character values for y
    dp['y_label'] = pd.Series(dtype='object')
    
    # For CG-MH items (not aggregates) - Ukraine uses 0-3 scale
    cg_mh_mask = (
        dp['item_label'].str.contains('^CG-MH_', regex=True) & 
        ~dp['item_label'].str.contains('agg')
    )
    y_labels_ukraine = [
        'a - not at all',
        'b - several days',
        'c - more than half of the time',
        'd - nearly every day'
    ]
    dp.loc[cg_mh_mask, 'y_label'] = dp.loc[cg_mh_mask, 'y'].apply(
        lambda y: y_labels_ukraine[int(y)] if pd.notna(y) and 0 <= int(y) < 4 else np.nan
    )
    
    # For other items (not aggregates)
    other_mask = dp['y_label'].isna() & ~dp['item_label'].str.contains('agg')
    dp.loc[other_mask, 'y_label'] = dp.loc[other_mask, 'y'].apply(
        lambda y: f"{int(y)} of 7 days" if pd.notna(y) else np.nan
    )
    
    # Create item metadata table
    dit = dp[['item_label']].drop_duplicates().sort_values('item_label').reset_index(drop=True)
    
    # Set item_type
    dit['item_type'] = np.where(
        dit['item_label'].str.contains('CG-MH'),
        'categorical',
        'out-of-7'
    )
    
    # Set item_high_label
    dit['item_high_label'] = np.where(
        dit['item_label'].str.contains('CG-MH|CG-DEPRESSION|CG-VIO|CHI-BEHAVIOUR', regex=True),
        'lower_is_better',
        'higher_is_better'
    )
    
    # Extract group_label and item_label_short
    dit['group_label'] = dit['item_label'].str.replace(r'([^_]+)_([^_]+)', r'\1', regex=True)
    dit['item_label_short'] = dit['item_label'].str.replace(r'([^_]+)_([^_]+)', r'\2', regex=True)
    
    # Clean up item_label_short for CG items
    dit.loc[dit['item_label_short'].str.startswith('CG'), 'item_label_short'] = ''
    dit['item_label_short'] = dit['item_label_short'].replace('', np.nan)
    
    # Create group_label_long
    group_mapping = {
        'CG-MH': 'Caregiver mental health',
        'CG-VIO': 'Caregiver exercising physical or emotional violence',
        'CG-MONITOR-CHI_agg': 'Child monitoring',
        'CG-INVOLVE': 'Caregiver involvement',
        'CHI-BEHAVIOUR': 'Child behavioural issues',
        'CG-DEPRESSION': 'Caregiver depression',
        'CG-SELFCARE': 'Caregiver self-care',
        'CG-RESILIENCE': 'Caregiver resilience',
        'CG-POS': 'Caregiver positive parenting',
        'CG-MONITOR-CHI': 'Caregiver monitoring child',
        'CG-NONVIOLENT-DISCIPLINE': 'Caregiver exercising nonviolent discipline'
    }
    dit['group_label_long'] = dit['group_label'].map(group_mapping).fillna(dit['group_label'])
    
    # Set endpoint_measure
    dit['endpoint_measure'] = dit['item_type'].map({
        'categorical': 'events occurring most or all of the time',
        'out-of-7': 'mean days in week'
    })
    
    # Set cat_length - Ukraine: 4 for categorical (vs Colombia: 5)
    dit['cat_length'] = dit['item_type'].map({
        'categorical': 4,
        'out-of-7': 8
    }).astype('Int64')

    # Set item_type_id
    type_mapping = pd.DataFrame({
        'item_type': sorted(dit['item_type'].unique())
    })
    type_mapping['item_type_id'] = range(1, len(type_mapping) + 1)
    dit = dit.merge(type_mapping, on='item_type')

    # Reset indices
    dp = dp.reset_index(drop=True)
    dit = dit.reset_index(drop=True)
    dmeta = dmeta.reset_index(drop=True)

    return {'dp': dp, 'dit': dit, 'dmeta': dmeta}


def _common_dit(items, item_type, high_low, group_of, group_long, endpoint, cat_length):
    """Assemble a `dit` item-metadata table in the common bIRTistic format."""
    dit = pd.DataFrame({'item_label': list(items)})
    dit['item_type'] = item_type
    dit['item_high_label'] = [high_low(i) for i in dit['item_label']]
    dit['group_label'] = [group_of(i) for i in dit['item_label']]
    dit['item_label_short'] = dit['item_label']
    dit['group_label_long'] = [group_long.get(g, g) for g in dit['group_label']]
    dit['endpoint_measure'] = endpoint
    dit['cat_length'] = pd.array([cat_length] * len(dit), dtype='Int64')
    tmap = pd.DataFrame({'item_type': sorted(dit['item_type'].unique())})
    tmap['item_type_id'] = range(1, len(tmap) + 1)
    return dit.merge(tmap, on='item_type').reset_index(drop=True)


# Mycelium acceptance items (1-7 Likert): positive constructs vs reverse-direction
# risk/disgust. Recoded duplicates (hearR, envrR), composites (Env, EnvD), the Q17
# battery and design/demographic columns are excluded from the item set.
_MYCELIUM_HIGH = ['INT1', 'INT2', 'ATT1', 'ATT2', 'SOC1', 'SOC2', 'SOC3',
                  'HeaB', 'EnvB', 'NATURAL', 'Familiarity']
_MYCELIUM_LOW = ['DISG1', 'DISG2', 'DISG3', 'DISG4', 'HeaR', 'EnvR']
_MYCELIUM_GROUP_LONG = {'INT': 'Intention to eat', 'ATT': 'Attitude', 'SOC': 'Social norms',
                        'DISG': 'Disgust', 'HeaB': 'Health benefit', 'EnvB': 'Environmental benefit',
                        'HeaR': 'Health risk', 'EnvR': 'Environmental risk',
                        'NATURAL': 'Naturalness', 'Familiarity': 'Familiarity'}


def read_data_mycelium(file_data: str) -> Dict[str, pd.DataFrame]:
    """Read the mycelium novel-food acceptance survey (doc §3.14) into common format.

    Cross-sectional 3x3 (processing x substrate) survey, UK Prolific N=449, item-level
    1-7 Likert acceptance constructs (Zenodo 10628634).

    Parameters
    ----------
    file_data : str
        Path to ``Mycelium.csv``.

    Returns
    -------
    dict
        'dp'  : long format (time, time_label, pid, pid_label, fid, f_label,
                submission_date, treat, item_label, y, y_label);
        'dit' : item metadata; 'dmeta': age, gender, education, condition, processing, substrate.
    """
    import re
    if not Path(file_data).exists():
        raise FileNotFoundError(f"Data file not found: {file_data}")
    raw = pd.read_csv(file_data, encoding='utf-8-sig', low_memory=False)
    raw.columns = [str(c).replace('﻿', '').replace('ï»¿', '').strip() for c in raw.columns]
    raw = raw.reset_index(drop=True)
    raw['pid'] = np.arange(1, len(raw) + 1)

    items = [c for c in _MYCELIUM_HIGH + _MYCELIUM_LOW if c in raw.columns]
    cond = pd.to_numeric(raw['Condition'], errors='coerce') if 'Condition' in raw.columns else pd.Series(np.nan, index=raw.index)
    pid2cond = dict(zip(raw['pid'], cond))

    # dmeta: design + demographics
    dmeta = raw[['pid']].copy()
    for src, out in [('Age', 'age'), ('Gender', 'gender'), ('Education', 'education'),
                     ('Condition', 'treat'), ('Process', 'processing'), ('Source', 'substrate')]:
        if src in raw.columns:
            dmeta[out] = pd.to_numeric(raw[src], errors='coerce')

    # dp long
    long = raw[['pid'] + items].melt(id_vars='pid', var_name='item_label', value_name='y')
    long['y'] = pd.to_numeric(long['y'], errors='coerce')
    long = long.dropna(subset=['y'])
    dp = pd.DataFrame({
        'time': 0, 'time_label': 'survey',
        'pid': long['pid'].to_numpy(), 'pid_label': long['pid'].astype(str).to_numpy(),
        'fid': np.nan, 'f_label': np.nan, 'submission_date': pd.NaT,
        'treat': long['pid'].map(pid2cond).to_numpy(),
        'item_label': long['item_label'].to_numpy(), 'y': long['y'].astype(float).to_numpy(),
        'y_label': long['y'].astype(int).astype(str).to_numpy(),
    })

    dit = _common_dit(
        items, 'likert-7',
        high_low=lambda i: 'lower_is_better' if i in _MYCELIUM_LOW else 'higher_is_better',
        group_of=lambda i: re.match(r'([A-Za-z]+)', i).group(1),
        group_long=_MYCELIUM_GROUP_LONG,
        endpoint='mean 7-point acceptance rating', cat_length=7)
    return {'dp': dp.reset_index(drop=True), 'dit': dit, 'dmeta': dmeta.reset_index(drop=True)}


# REFUGE-ED MSPSS: drop the four subscale-mean / total columns, keep the 12 items.
_REFUGE_MSPSS_DROP = {'MSPSS_Mean', 'MSPSS_SO', 'MSPSS_Fam', 'MSPSS_Fri'}


def read_data_refuge_ed(file_data: str) -> Dict[str, pd.DataFrame]:
    """Read REFUGE-ED Youth Baseline & Endline (doc §3.9) into common format.

    Refugee/migrant youth, baseline+endline, item-level 1-7 MSPSS perceived-social-support
    items (Zenodo 10908209). The four MSPSS subscale means/totals are excluded.

    Parameters
    ----------
    file_data : str
        Path to ``Youth Baseline & Endline .xlsx`` (requires ``openpyxl``).

    Returns
    -------
    dict
        'dp' (long, time 0=baseline / 1=endline), 'dit', 'dmeta' (country, site, gender, age).
    """
    if not Path(file_data).exists():
        raise FileNotFoundError(f"Data file not found: {file_data}")
    xl = pd.ExcelFile(file_data)
    bl = xl.parse('Baseline Data').dropna(how='all'); bl['time'] = 0
    el = xl.parse('Endline Data').dropna(how='all'); el['time'] = 1
    raw = pd.concat([bl, el], ignore_index=True)
    pc = 'Participant Code'

    def clean_item(c):
        s = pd.to_numeric(raw[c], errors='coerce').dropna()
        return len(s) > 0 and s.min() >= 1 and s.max() <= 7
    items = [c for c in raw.columns
             if str(c).startswith('MSPSS_') and c not in _REFUGE_MSPSS_DROP and clean_item(c)]
    raw['pid'] = raw[pc].astype('category').cat.codes + 1

    meta_cols = {'Country Code': 'country', 'Pilot Site': 'site', 'Gender': 'gender', 'Y_Age': 'age'}
    dmeta = raw[['pid', pc, 'time'] + [c for c in meta_cols if c in raw.columns]].rename(
        columns={pc: 'pid_label', **meta_cols})

    long = raw[['pid', pc, 'time'] + items].melt(
        id_vars=['pid', pc, 'time'], var_name='item_label', value_name='y')
    long['y'] = pd.to_numeric(long['y'], errors='coerce')
    long = long.dropna(subset=['y'])
    dp = pd.DataFrame({
        'time': long['time'].to_numpy(),
        'time_label': long['time'].map({0: 'Baseline', 1: 'Endline'}).to_numpy(),
        'pid': long['pid'].to_numpy(), 'pid_label': long[pc].astype(str).to_numpy(),
        'fid': np.nan, 'f_label': np.nan, 'submission_date': pd.NaT, 'treat': 0.0,
        'item_label': long['item_label'].to_numpy(), 'y': long['y'].astype(float).to_numpy(),
        'y_label': long['y'].astype(int).astype(str).to_numpy(),
    })

    # MSPSS has three subscales (Family/Friends/Significant Other); group by subscale
    def _subscale(i):                                          # 'MSPSS_Fam5' -> 'Fam'
        return re.sub(r'\d+$', '', i.replace('MSPSS_', ''))
    dit = _common_dit(
        items, 'likert-7',
        high_low=lambda i: 'higher_is_better',                 # more perceived support = better
        group_of=lambda i: f'Perceived social support: {_subscale(i)}',
        group_long={},                                         # identity: group_label_long == group_label
        endpoint='mean 7-point perceived-support rating', cat_length=7)
    dit['item_label_short'] = [re.sub(r'^MSPSS_(Fam|Fri|SO)', '', c) for c in dit['item_label']]  # subscale in group; keep item number
    return {'dp': dp.reset_index(drop=True), 'dit': dit, 'dmeta': dmeta.reset_index(drop=True)}


# Temporal-dynamics symptom scales (doc §3.5). All score higher = worse (lower_is_better).
_TEMPORAL_SCALES = {'phq9': 'Patient Health Questionnaire (depression)',
                    'gad7': 'Generalised Anxiety Disorder scale',
                    'isi': 'Insomnia Severity Index',
                    'pss': 'Perceived Stress Scale'}


def read_data_temporal_dynamics(file_data: str) -> Dict[str, pd.DataFrame]:
    """Read the temporal-dynamics psychological item bank (doc §3.5) into common format.

    Cross-sectional; four self-report symptom scales at the item level (Zenodo 10423537).

    Parameters
    ----------
    file_data : str
        Directory holding phq9.csv, gad7.csv, isi.csv, pss.csv
        (each: export_id, score, question1.., time1..).

    Returns
    -------
    dict with 'dp', 'dit' (per-scale cat_length inferred from the data), 'dmeta'.
    """
    d = Path(file_data)
    if not d.exists():
        raise FileNotFoundError(f"Directory not found: {file_data}")
    dps, dits = [], []
    for scale, longname in _TEMPORAL_SCALES.items():
        fp = d / f"{scale}.csv"
        if not fp.exists():
            continue
        raw = pd.read_csv(fp)
        qcols = [c for c in raw.columns if str(c).startswith('question')]
        sub = raw[['export_id'] + qcols].copy()
        for c in qcols:
            sub[c] = pd.to_numeric(sub[c], errors='coerce')
        cat_len = int(np.nanmax(sub[qcols].to_numpy())) + 1        # infer K (values 0..K-1)
        ren = {c: f"{scale}_q{int(str(c).replace('question', ''))}" for c in qcols}
        long = sub.rename(columns=ren).melt(
            id_vars='export_id', var_name='item_label', value_name='y').dropna(subset=['y'])
        dps.append(long)
        di = pd.DataFrame({'item_label': list(ren.values())})
        di['item_type'] = scale
        di['item_high_label'] = 'lower_is_better'
        di['group_label'] = scale.upper()
        di['item_label_short'] = [l.split('_')[-1] for l in di['item_label']]
        di['group_label_long'] = longname
        di['endpoint_measure'] = 'mean symptom item score'
        di['cat_length'] = pd.array([cat_len] * len(di), dtype='Int64')
        dits.append(di)
    long = pd.concat(dps, ignore_index=True)
    pmap = {e: i + 1 for i, e in enumerate(sorted(long['export_id'].unique()))}
    dp = pd.DataFrame({
        'time': 0, 'time_label': 'survey',
        'pid': long['export_id'].map(pmap).to_numpy(),
        'pid_label': long['export_id'].astype(str).to_numpy(),
        'fid': np.nan, 'f_label': np.nan, 'submission_date': pd.NaT, 'treat': np.nan,
        'item_label': long['item_label'].to_numpy(), 'y': long['y'].astype(float).to_numpy(),
        'y_label': long['y'].astype(int).astype(str).to_numpy(),
    })
    dit = pd.concat(dits, ignore_index=True)
    tmap = pd.DataFrame({'item_type': sorted(dit['item_type'].unique())})
    tmap['item_type_id'] = range(1, len(tmap) + 1)
    dit = dit.merge(tmap, on='item_type').reset_index(drop=True)
    dmeta = pd.DataFrame({'pid': list(pmap.values()), 'pid_label': list(pmap.keys())})
    return {'dp': dp.reset_index(drop=True), 'dit': dit, 'dmeta': dmeta}


_CHATGPT_SURVEY = {'Pre-Intervention Survey on Critical Approach to AI.xlsx': 0,
                   'Post-Intervention Survey on Critical Approach to AI.xlsx': 1}


def read_data_chatgpt_rct(file_data: str) -> Dict[str, pd.DataFrame]:
    """Read the ChatGPT-vs-expert-feedback RCT survey (doc §3.7) into common format.

    Two-arm RCT; the six 'critical approach to AI' 1-7 Likert items, pre (time 0) and
    post (time 1). The partial-credit key-feature test items are not loaded here.

    Parameters
    ----------
    file_data : str
        Directory holding the Pre-/Post-Intervention Survey xlsx files (Zenodo 13769970).
        Requires openpyxl.

    Returns
    -------
    dict with 'dp' (time 0=pre / 1=post, treat = arm), 'dit', 'dmeta'.
    """
    d = Path(file_data)
    if not d.exists():
        raise FileNotFoundError(f"Directory not found: {file_data}")
    frames = []
    for fname, t in _CHATGPT_SURVEY.items():
        fp = d / fname
        if not fp.exists():
            continue
        raw = pd.ExcelFile(fp).parse('Survey')
        cols = list(raw.columns)
        items = cols[4:]                                    # 6 statement items (after ID/Group/Gender/Repeat)
        sub = raw[[cols[0], cols[1]] + items].copy()
        sub.columns = ['ID', 'treat'] + [f"AI_q{i + 1}" for i in range(len(items))]
        sub['time'] = t
        frames.append(sub)
    alld = pd.concat(frames, ignore_index=True)
    qs = [c for c in alld.columns if c.startswith('AI_q')]
    long = alld.melt(id_vars=['ID', 'treat', 'time'], value_vars=qs,
                     var_name='item_label', value_name='y')
    long['y'] = pd.to_numeric(long['y'], errors='coerce')
    long = long.dropna(subset=['y'])
    dp = pd.DataFrame({
        'time': long['time'].to_numpy(),
        'time_label': long['time'].map({0: 'Pre', 1: 'Post'}).to_numpy(),
        'pid': pd.to_numeric(long['ID'], errors='coerce').astype('Int64').to_numpy(),
        'pid_label': long['ID'].astype(str).to_numpy(),
        'fid': np.nan, 'f_label': np.nan, 'submission_date': pd.NaT,
        'treat': pd.to_numeric(long['treat'], errors='coerce').to_numpy(),
        'item_label': long['item_label'].to_numpy(), 'y': long['y'].astype(float).to_numpy(),
        'y_label': long['y'].astype(int).astype(str).to_numpy(),
    })
    dit = _common_dit(
        qs, 'likert-7', high_low=lambda i: 'higher_is_better', group_of=lambda i: 'AI-attitude',
        group_long={'AI-attitude': 'Critical approach to AI'},
        endpoint='mean 7-point attitude rating', cat_length=7)
    dmeta = alld[['ID', 'treat']].drop_duplicates().rename(columns={'ID': 'pid'}).reset_index(drop=True)
    return {'dp': dp.reset_index(drop=True), 'dit': dit, 'dmeta': dmeta}


def read_data_pisa(pisa_dir: str, country: str,
                   cycles=(2012, 2015, 2018, 2022)) -> Dict[str, pd.DataFrame]:
    """Read the PISA Math extract (doc §3.11) for one country into common format.

    Cross-sectional international assessment; content-id-matched Math items across
    cycles (see ``data_web_extracting.build_pisa_math_extract``). Baseline-anchored:
    ``cycles[0]`` (2012) is the baseline, each later cycle a subsequent interim.

    Parameters
    ----------
    pisa_dir : str
        Directory holding ``pisa_math_{cycle}.parquet`` (+ ``common_math.json``).
    country : str
        PISA CNT code (e.g. ``'USA'``).

    Returns
    -------
    dict with 'dp' (long: cycle, pid, item_label, y) and 'dit' (item metadata,
    per-item ``cat_length`` = max category + 1 across the country's cycles).
    """
    frames = []
    for yr in cycles:
        fp = Path(pisa_dir) / f"pisa_math_{yr}.parquet"
        if not fp.exists():
            raise FileNotFoundError(f"missing PISA extract: {fp}")
        d = pd.read_parquet(fp)
        d = d[d['CNT'] == country]
        frames.append(d)
    long = pd.concat(frames, ignore_index=True)
    idc = next(c for c in long.columns if 'STU' in c.upper())
    dp = pd.DataFrame({
        'cycle': long['cycle'].to_numpy(),
        'pid_label': (long['cycle'].astype(str) + '_' + long[idc].astype(str)).to_numpy(),
        'item_label': long['item'].to_numpy(),
        'y': long['y'].astype(int).to_numpy(),
    })
    dp['pid'] = pd.factorize(dp['pid_label'])[0] + 1
    kmax = dp.groupby('item_label')['y'].max()
    items = sorted(dp['item_label'].unique())
    dit = _common_dit(
        items, 'out-of-7',
        high_low=lambda i: 'higher_is_better',            # higher score = better performance
        group_of=lambda i: 'PISA Math',
        group_long={'PISA Math': 'PISA Math'},
        endpoint='mean item score (baseline vs cycle)', cat_length=1)
    dit['cat_length'] = dit['item_label'].map(lambda i: int(kmax[i]) + 1).astype('Int64')
    return {'dp': dp.reset_index(drop=True), 'dit': dit}
