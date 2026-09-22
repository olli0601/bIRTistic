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


def _dilution_ladder_labels(K, below_detection, continuous=False):
    """Scientifically accurate, reciprocal-dilution category labels for a serial 2-fold antibody
    titre ladder (HAI, neutralisation ID50, ...) shared across the immunogenicity applications.

    Ordered category ``k`` maps to titre ``below_detection * 2**k`` (the coding used by the
    loaders, ``k = round(log2(titre / below_detection))``). For a genuine discrete serial-dilution
    assay (HAI), ``below_detection`` is the imputed value for a sample that fails the LOWEST real
    dilution (half of it), so category 0 is the below-detection floor:
        k=0 -> ``"<1:{2*below_detection}"`` (e.g. ``<1:10``),  k>=1 -> ``"1:{below_detection*2**k}"``
    giving ``['<1:10','1:10','1:20','1:40','1:80', ...]`` (``1:40`` == k3 == the CHMP
    seroprotection threshold). For a CONTINUOUS readout binned into 2-fold classes (neutralisation
    ID50), pass ``continuous=True``: there is no below-detection imputation, so category 0 is the
    assay's lowest reciprocal dilution itself (``"1:{below_detection}"``) and every class is
    ``"1:{below_detection*2**k}"``."""
    if continuous:
        return [f"1:{int(round(below_detection * (2 ** k)))}" for k in range(K)]
    lab = [f"<1:{int(round(2 * below_detection))}"]
    lab += [f"1:{int(round(below_detection * (2 ** k)))}" for k in range(1, K)]
    return lab


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
        'Timepoint': 'group_label',
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
    dmeta = dp[['pid', 'group_label'] + col_meta].copy()
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
    dp['group'] = (dp['group_label'] == 'Endline').astype(int)
    
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
    id_vars = ['group', 'group_label', 'pid', 'pid_label', 'fid', 
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
    
    # Extract construct from item_label
    dit['construct'] = dit['item_label'].str.replace(r'([^_]+)_([^_]+)', r'\1', regex=True)
    
    # Extract item_label_short
    dit['item_label_short'] = dit['item_label'].str.replace(r'([^_]+)_([^_]+)', r'\2', regex=True)
    # Set to NaN if starts with 'CG'
    dit.loc[dit['item_label_short'].str.startswith('CG'), 'item_label_short'] = np.nan
    
    # Create construct_long with full names
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
    dit['construct_long'] = dit['construct'].map(group_mapping).fillna(dit['construct'])
    
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
        'Timepoint': 'group_label',
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
    dmeta = dp[['pid', 'group_label'] + col_meta].copy()
    
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
    dp['group'] = (dp['group_label'] == 'Endline').astype(int)
    
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
    id_vars = ['group', 'group_label', 'pid', 'pid_label', 'fid', 
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
    
    # Extract construct and item_label_short
    dit['construct'] = dit['item_label'].str.replace(r'([^_]+)_([^_]+)', r'\1', regex=True)
    dit['item_label_short'] = dit['item_label'].str.replace(r'([^_]+)_([^_]+)', r'\2', regex=True)
    
    # Clean up item_label_short for CG items
    dit.loc[dit['item_label_short'].str.startswith('CG'), 'item_label_short'] = ''
    dit['item_label_short'] = dit['item_label_short'].replace('', np.nan)
    
    # Create construct_long
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
    dit['construct_long'] = dit['construct'].map(group_mapping).fillna(dit['construct'])
    
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
    dit['construct'] = [group_of(i) for i in dit['item_label']]
    dit['item_label_short'] = dit['item_label']
    dit['construct_long'] = [group_long.get(g, g) for g in dit['construct']]
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
        'dp'  : long format (time, group_label, pid, pid_label, fid, f_label,
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
        'group': 0, 'group_label': 'survey',
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
    bl = xl.parse('Baseline Data').dropna(how='all'); bl['group'] = 0
    el = xl.parse('Endline Data').dropna(how='all'); el['group'] = 1
    raw = pd.concat([bl, el], ignore_index=True)
    pc = 'Participant Code'

    def clean_item(c):
        s = pd.to_numeric(raw[c], errors='coerce').dropna()
        return len(s) > 0 and s.min() >= 1 and s.max() <= 7
    items = [c for c in raw.columns
             if str(c).startswith('MSPSS_') and c not in _REFUGE_MSPSS_DROP and clean_item(c)]
    raw['pid'] = raw[pc].astype('category').cat.codes + 1

    meta_cols = {'Country Code': 'country', 'Pilot Site': 'site', 'Gender': 'gender', 'Y_Age': 'age'}
    dmeta = raw[['pid', pc, 'group'] + [c for c in meta_cols if c in raw.columns]].rename(
        columns={pc: 'pid_label', **meta_cols})

    long = raw[['pid', pc, 'group'] + items].melt(
        id_vars=['pid', pc, 'group'], var_name='item_label', value_name='y')
    long['y'] = pd.to_numeric(long['y'], errors='coerce')
    long = long.dropna(subset=['y'])
    dp = pd.DataFrame({
        'group': long['group'].to_numpy(),
        'group_label': long['group'].map({0: 'Baseline', 1: 'Endline'}).to_numpy(),
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
        group_long={},                                         # identity: construct_long == construct
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
        di['construct'] = scale.upper()
        di['item_label_short'] = [l.split('_')[-1] for l in di['item_label']]
        di['construct_long'] = longname
        di['endpoint_measure'] = 'mean symptom item score'
        di['cat_length'] = pd.array([cat_len] * len(di), dtype='Int64')
        dits.append(di)
    long = pd.concat(dps, ignore_index=True)
    pmap = {e: i + 1 for i, e in enumerate(sorted(long['export_id'].unique()))}
    dp = pd.DataFrame({
        'group': 0, 'group_label': 'survey',
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
        sub['group'] = t
        frames.append(sub)
    alld = pd.concat(frames, ignore_index=True)
    qs = [c for c in alld.columns if c.startswith('AI_q')]
    long = alld.melt(id_vars=['ID', 'treat', 'group'], value_vars=qs,
                     var_name='item_label', value_name='y')
    long['y'] = pd.to_numeric(long['y'], errors='coerce')
    long = long.dropna(subset=['y'])
    dp = pd.DataFrame({
        'group': long['group'].to_numpy(),
        'group_label': long['group'].map({0: 'Pre', 1: 'Post'}).to_numpy(),
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


def read_data_immport_flu(xlsx_path, study, endline_day=None, min_start_dilution=5.0, arm=None):
    """ImmuneSpace/ImmPort influenza HAI -> paired PCM frame for one study (SDY).

    Strains are the items; day 0 (Baseline) vs the post-vaccination visit (Endline,
    default = the largest of {28,27,24,21} present) are the paired time axis; the HAI
    titre is mapped to an ORDERED category k = round(log2(titre / min_start_dilution))
    (2-fold serial dilution ladder), with a single, study-wide K across strains (they
    share one assay -> one item_type_id, no K-mixing, cf. PISA). Keeps only
    participants with both time-points (paired). Endpoint via get_endpoints as the
    'out-of-7' expected-category ratio (higher titre is better).

    `arm` (optional, case-insensitive substring of the Arm `Name`) restricts to one
    vaccination arm for a head-to-head study (e.g. SDY269 LAIV vs TIV): the per-arm
    paired seroresponse (GMFR/SPR) is then compared across arms. Note the arms may be
    assayed on DIFFERENT strain panels (only shared strains are directly comparable).

    Returns dict(dp=dp1, dit=dit) in the same shape the SVI producers consume.
    """
    import re
    asy = pd.read_excel(xlsx_path, sheet_name='Assays')
    ev = pd.read_excel(xlsx_path, sheet_name='Events')
    hai = asy[(asy['Study ID'] == study)
              & asy['Assay Subtype'].astype(str).str.contains('hemagglut', case=False)].copy()
    if hai.empty:
        raise ValueError(f'{study}: no HAI rows in {xlsx_path}')
    if arm is not None:                                  # restrict to one vaccination arm (head-to-head)
        arms = pd.read_excel(xlsx_path, sheet_name='Arms')
        part = pd.read_excel(xlsx_path, sheet_name='Participants')
        anames = arms[arms['Study ID'] == study][['Arm ID', 'Name']]
        pa = part[part['Study ID'] == study][['Participant ID', 'Arm ID']].merge(anames, on='Arm ID', how='left')
        keep = set(pa.loc[pa['Name'].astype(str).str.contains(arm, case=False, na=False), 'Participant ID'])
        if not keep:
            raise ValueError(f"{study}: no participants in arm matching {arm!r}; "
                             f"arm names present: {sorted(anames['Name'].astype(str).unique())}")
        hai = hai[hai['Participant ID'].isin(keep)].copy()
    evd = ev[['Event ID', 'Start', 'Participant ID']].rename(columns={'Start': 'day'})
    h = hai.merge(evd, on=['Event ID', 'Participant ID'], how='left')
    h['titre'] = pd.to_numeric(h['Value'], errors='coerce')
    h = h.dropna(subset=['titre', 'day'])
    if endline_day is None:
        present = set(int(d) for d in h['day'].unique())
        endline_day = next((d for d in (28, 27, 24, 21) if d in present), max(present - {0}))
    h = h[h['day'].isin([0, int(endline_day)])].copy()
    # titre -> ordered category on the 2-fold ladder; single study-wide K
    k = np.clip(np.rint(np.log2(h['titre'] / float(min_start_dilution))), 0, None).astype(int)
    h['k'] = k
    K = int(h['k'].max()) + 1
    # short strain label: the designation inside 'Influenza ... virus (...)'
    def _strain(s):
        m = re.search(r'virus \((.+)\)\s*$', str(s))
        return (m.group(1) if m else str(s)).strip()
    h['item_label'] = h['Target Entity Subtype'].map(_strain)
    # keep only paired participants (both day 0 and endline), 1 row per (pid,strain,day)
    h = h.sort_values(['Participant ID', 'item_label', 'day']).drop_duplicates(
        ['Participant ID', 'item_label', 'day'], keep='last')
    tp = h.groupby('Participant ID')['day'].nunique()
    paired = tp[tp >= 2].index
    h = h[h['Participant ID'].isin(paired)].copy()
    # build dp1
    h['group'] = (h['day'] == int(endline_day)).astype(int)
    h['pid_label'] = h['Participant ID'].astype(str)
    h['pid'] = pd.factorize(h['pid_label'])[0] + 1
    dp1 = pd.DataFrame({
        'pid': h['pid'].to_numpy(), 'pid_label': h['pid_label'].to_numpy(),
        'group': h['group'].to_numpy(),
        'group_label': h['group'].map({0: 'Baseline', 1: 'Endline'}).to_numpy(),
        'fid': np.nan, 'f_label': np.nan, 'submission_date': pd.NaT, 'treat': 0.0,
        'item_label': h['item_label'].to_numpy(), 'y': h['k'].to_numpy(),
    })
    # category label = the reciprocal HAI titre on the 2-fold dilution ladder; k0 (value =
    # min_start_dilution = below the lowest 1:10 dilution) reads "<1:10", then 1:10, 1:20, 1:40
    # (== k3 == seroprotection), ... The plot fixes the discrete x-axis order to 0..K-1 via limits.
    import json, re as _re
    titre_labels = _dilution_ladder_labels(K, min_start_dilution)   # ['<1:10','1:10','1:20',...]
    def _short_strain(name):                            # "A/Puerto Rico/8/1934(H1N1)" -> "A/Puerto Rico/34 (H1N1)"
        m = _re.match(r'([AB])/(.+?)/(?:\d+/)?(\d{2,4})\s*(\(.+\))?$', str(name))
        if not m:
            return str(name)
        typ, loc, yr, sub = m.groups()
        return f"{typ}/{loc}/{yr[-2:]}{(' ' + sub) if sub else ''}"
    dp1['y_stan'] = dp1['y'] + 1                        # 1..K model input
    dp1['y_label'] = dp1['y'].map(lambda k: titre_labels[int(k)])
    # 'categorical' (ordered HAI dilution levels) -> endpoints are computed via rho_specs
    # (seroprotection rate + GMT fold-rise), not the 'out-of-7' expected-score branch
    dp1['item_type'] = 'categorical'; dp1['item_type_id'] = 1
    cat_labels = json.dumps(titre_labels)               # complete 0..K-1 map -> no NaN
    strains = sorted(dp1['item_label'].unique())
    dit = pd.DataFrame({'item_label': strains})
    dit['item_type'] = 'categorical'; dit['item_type_id'] = 1; dit['cat_length'] = K
    dit['item_label_short'] = np.nan                    # item_label_long == construct_long (short strain)
    dit['construct'] = 'HAI titre'
    dit['construct_long'] = dit['item_label'].map(_short_strain)  # short strain per facet row
    dit['item_high_label'] = 'higher_is_better'
    dit['endpoint_measure'] = f'{study} HAI: titre (log2 dilution)'
    dit['cat_labels'] = cat_labels
    return {'dp': dp1, 'dit': dit, 'K': K, 'endline_day': int(endline_day)}


def read_data_immport_flu_headhead(xlsx_path, study, arms=('LAIV', 'TIV'),
                                   min_start_dilution=5.0, endline_day=None):
    """ImmPort influenza HAI HEAD-TO-HEAD -> ONE JOINT paired PCM frame over both vaccine arms.

    Rather than fitting each arm separately, this pools both arms into a single PCM so the arms
    are estimated in one go and their endpoints are correlated within a posterior draw (the shared
    strain's arm-difference then needs no random draw-pairing). The **item** is the (strain, arm)
    pair, distinguished in `item_label` as ``"<strain> [<arm>]"`` (with `strain`/`arm` kept as
    their own columns); the paired day-0/endline axis is `group` (0/1) as in
    :func:`read_data_immport_flu`. Every item shares one `item_type_id` and a single, study-wide
    `K` = max over arms (so the incremental thresholds live on one ladder -> no K-mixing), and
    the participant abilities $\\theta_i$ are shared across arms (participants are disjoint people,
    re-indexed to a common `pid` space). Each participant only responds to its own arm's strains,
    so the joint design is an incomplete block — standard for the PCM (only observed responses
    contribute). The (strain, arm, time) difficulty each get their own `item_group_id` downstream
    (built by the producer from `item_label` x `group`), which is exactly the non-reduced structure
    `get_endpoints_per_draw` consumes to emit SPR/GMFR per arm-strain and the cross-arm difference.

    **Schema convention.** `item_label` stays the EXACT response item (the strain), so a strain
    assayed in both arms (A/Uruguay) is ONE item, not two. The vaccine arm folds into the flexible
    condition axis together with the paired time-point: `group_label` in
    {LAIV_baseline, LAIV_endline, TIV_baseline, TIV_endline} (`group` its 0-based id), and the full
    non-reduced (strain x arm x time) structure is the CROSS `item_group_id` = item_label x group
    (built by the producer). Two helper columns decompose the condition for the endpoint machinery:
    `arm` (the stratum) and `phase` (0/1 baseline/endline — the within-arm paired contrast the
    endpoints pivot on). The fit is driven by `item_group_id` (the same partition either way) and
    `x_formula="~ phase - 1"`, so it is invariant to this relabelling.

    Returns dict(dp, dit, K, arms, shared) with `shared` = strains assayed in BOTH arms (the ones
    with a direct strain-for-strain arm contrast; for SDY269 = A/Uruguay/716/2007(H3N2))."""
    import json
    per = []
    for arm in arms:
        d = read_data_immport_flu(xlsx_path, study, endline_day=endline_day,
                                  min_start_dilution=min_start_dilution, arm=arm)
        per.append((arm, d))
    K = max(d['K'] for _, d in per)
    cat_labels = json.dumps(_dilution_ladder_labels(K, min_start_dilution))  # ['<1:10','1:10',...]
    dps, dits, pid_off = [], [], 0
    for arm, d in per:
        dp = d['dp'].copy(); dit = d['dit'].copy()
        # item_label stays the pure strain; the arm becomes part of the flexible condition axis
        dp['arm'] = arm
        dp['phase'] = dp['group']; dp['phase_label'] = dp['group_label']    # 0/1 baseline/endline
        dp['group_label'] = arm + '_' + dp['phase_label'].str.lower()       # e.g. TIV_baseline
        dp['pid'] = dp['pid'] + pid_off; dp['pid_label'] = arm + ':' + dp['pid_label'].astype(str)
        pid_off = int(dp['pid'].max())
        dit['cat_length'] = K; dit['cat_labels'] = cat_labels               # global K across all items
        dps.append(dp); dits.append(dit)
    dp1 = pd.concat(dps, ignore_index=True)
    # numeric condition id (0..3) from the arm x phase label; item stays 5 distinct strains
    order = [f'{a}_{p}' for a in arms for p in ('baseline', 'endline')]
    dp1['group'] = pd.Categorical(dp1['group_label'], categories=order, ordered=True).codes
    dp1['item_id'] = pd.factorize(dp1['item_label'])[0] + 1
    dit = pd.concat(dits, ignore_index=True).drop_duplicates('item_label').reset_index(drop=True)
    shared = sorted(set.intersection(*[set(d['dp']['item_label']) for _, d in per]))
    ref, foc = arms[0], arms[1]                                  # e.g. LAIV (reference), TIV
    # rho definitions declared UPFRONT with the application: each rho carries its own short
    # `rho_label` and pretty `rho_label_long` (used verbatim in plots), the get_endpoints params,
    # its H1 threshold `h1`, whether it is a per-arm `level` (restricted to `arm`) or a cross-arm
    # `diff`, and `is_level` for the composite "at least one met" rule. Built one-by-one downstream
    # but all indexed on the same draw, then joined -> supports joint decisions.
    # endpoint-major column order: both arms' SPR, then both arms' GMFR, then the two diffs
    rho_specs = []
    _lvl = [('spr', 'seroprotection rate  P(titre ≥ 1:40) at endline', 0.70,
             {'reduction': 'threshold', 'threshold': 3, 'compare': 'endline'}),
            ('gmfr', 'GMT fold-rise (endline / baseline)', 2.5,
             {'reduction': 'mean', 'compare': 'fold_log2'})]
    for tag, long, h1, params in _lvl:
        for a in arms:
            rho_specs.append({'rho_id': len(rho_specs) + 1, 'rho_label': f'{a}_{tag}',
                              'rho_label_long': f'{a} — {long}', 'kind': 'level', 'arm': a,
                              'is_level': True, 'h1': h1, **params})
    rho_specs += [
        {'rho_id': len(rho_specs) + 1, 'rho_label': f'{foc}-{ref}_spr_diff',
         'rho_label_long': f'{foc} − {ref} — seroprotection-rate difference (shared strain)',
         'kind': 'diff', 'is_level': False, 'h1': 0.0,
         'reduction': 'threshold', 'threshold': 3, 'compare': 'endline',
         'across': 'arm', 'across_within': 'item_label', 'across_values': [ref, foc],
         'across_compare': 'diff'},
        {'rho_id': len(rho_specs) + 2, 'rho_label': f'{foc}-{ref}_gmfr_diff',
         'rho_label_long': f'{foc} − {ref} — GMT fold-rise difference (shared strain)',
         'kind': 'diff', 'is_level': False, 'h1': 0.0,
         'reduction': 'mean', 'compare': 'fold_log2',
         'across': 'arm', 'across_within': 'item_label', 'across_values': [ref, foc],
         'across_compare': 'diff'},
    ]
    return {'dp': dp1, 'dit': dit, 'K': K, 'arms': list(arms), 'shared': shared,
            'rho_specs': rho_specs, 'endline_day': int(per[0][1]['endline_day'])}


def read_data_cavd_nab(nab_parquet, demo_parquet=None, min_start_dilution=5.0,
                       baseline_day=None, endline_day=None):
    """CAVD DataSpace neutralising-antibody (NAb) -> paired PCM frame for one HVTN study.

    The direct HIV analogue of read_data_immport_flu: the neutralised virus/isolate is the
    item, the pre/post visits are the paired time axis, and the ID50 neutralisation titre is
    a serial-dilution readout mapped to an ordered category k = round(log2(ID50 / start)),
    with a single study-wide K across isolates (one item_type_id -> no K-mixing). Keeps
    paired participants only. `demo_parquet` (optional) attaches arm/treat (study_group:
    1=vaccine, 2=placebo). Titres at the below-detection code (== start) map to k=0.

    Both parquet files are produced by data_web_extracting.download_cavd_dataspace, e.g.
    cavd_vtn505_NAb.parquet / cavd_vtn505_Demographics.parquet. Returns dict(dp, dit, K).
    """
    n = pd.read_parquet(nab_parquet)
    n = n.dropna(subset=['titer_ID50', 'visit_day', 'virus', 'SubjectId']).copy()
    days = sorted(n['visit_day'].unique())
    b = baseline_day if baseline_day is not None else days[0]
    e = endline_day if endline_day is not None else days[-1]
    n = n[n['visit_day'].isin([b, e])].copy()
    # ID50 -> ordered category on the 2-fold ladder; single study-wide K
    k = np.clip(np.rint(np.log2(n['titer_ID50'] / float(min_start_dilution))), 0, None).astype(int)
    n['k'] = k
    K = int(n['k'].max()) + 1
    n['item_label'] = n['virus'].astype(str)
    # 1 row per (subject, virus, visit); keep paired subjects (both visits)
    n = n.sort_values(['SubjectId', 'item_label', 'visit_day']).drop_duplicates(
        ['SubjectId', 'item_label', 'visit_day'], keep='last')
    tp = n.groupby('SubjectId')['visit_day'].nunique()
    n = n[n['SubjectId'].isin(tp[tp >= 2].index)].copy()
    # optional arm/treat
    treat = pd.Series(0.0, index=n.index); arm = pd.Series(np.nan, index=n.index)
    if demo_parquet is not None:
        dm = pd.read_parquet(demo_parquet)[['SubjectId', 'study_group', 'study_arm_summary']].drop_duplicates('SubjectId')
        dm = dm.set_index('SubjectId')
        treat = n['SubjectId'].map((dm['study_group'].astype(str) == '1').astype(float)).fillna(0.0)
        arm = n['SubjectId'].map(dm['study_arm_summary'])
    n['pid_label'] = n['SubjectId'].astype(str)
    n['pid'] = pd.factorize(n['pid_label'])[0] + 1
    dp1 = pd.DataFrame({
        'pid': n['pid'].to_numpy(), 'pid_label': n['pid_label'].to_numpy(),
        'group': (n['visit_day'] == e).astype(int).to_numpy(),
        'group_label': np.where(n['visit_day'] == e, 'Endline', 'Baseline'),
        'fid': n['SubjectId'].map(lambda s: s).to_numpy(), 'f_label': arm.to_numpy(),
        'submission_date': pd.NaT, 'treat': treat.to_numpy(),
        'item_label': n['item_label'].to_numpy(), 'y': n['k'].to_numpy(),
    })
    import json
    titre_labels = _dilution_ladder_labels(K, min_start_dilution)   # reciprocal ID50: <1:10, 1:10, ...
    dp1['fid'] = np.nan                                  # within-participant paired: no facilitator grouping
    dp1['y_stan'] = dp1['y'] + 1                         # 1..K model input
    dp1['y_label'] = dp1['y'].map(lambda k: titre_labels[int(k)])
    dp1['item_type'] = 'out-of-7'; dp1['item_type_id'] = 1
    isolates = sorted(dp1['item_label'].unique())
    dit = pd.DataFrame({'item_label': isolates})
    dit['item_type'] = 'out-of-7'; dit['item_type_id'] = 1; dit['cat_length'] = K
    dit['item_label_short'] = dit['item_label']
    dit['construct'] = 'HIV NAb ID50'; dit['construct_long'] = 'HIV-1 neutralisation ID50 (reciprocal titre)'
    dit['item_high_label'] = 'higher_is_better'
    dit['endpoint_measure'] = 'mean log2 ID50 neutralisation titre (endline vs baseline)'
    dit['cat_labels'] = json.dumps(titre_labels)
    return {'dp': dp1, 'dit': dit, 'K': K, 'baseline_day': int(b), 'endline_day': int(e)}


def read_data_cavd_bama_endline(bama_parquet, demo_parquet, endline_day=None, n_cat=3,
                                min_mfi=1.0):
    """CAVD DataSpace BAMA (binding-antibody) -> BETWEEN-ARM endline PCM frame for one HVTN
    study. HIV vaccine trials have no informative paired baseline (HIV-naive floor), so the
    estimand is the single-timepoint vaccine-vs-placebo contrast, encoded MYCELIUM-style:
    placebo -> time 0 ('Baseline'), vaccine -> time 1 ('Endline'); each subject in one arm.
    The antigen is the item; the background-subtracted magnitude mfi_delta is binned into
    n_cat ordered categories PER ANTIGEN (self-normalised qcut, so each antigen spans all
    n_cat levels). Antigens where either arm fails to span all n_cat categories are dropped
    (the codebase requires every item_time to observe every category). rho_j = relative
    vaccine-vs-placebo shift, higher_is_better. Returns dict(dp, dit, K, kept, dropped).

    Arm from Demographics.study_group (1=vaccine, 2=placebo). bama/demo parquet from
    data_web_extracting.download_cavd_dataspace (cavd_vtn505_BAMA/_Demographics.parquet).
    """
    b = pd.read_parquet(bama_parquet)
    dm = pd.read_parquet(demo_parquet)[['SubjectId', 'study_group']].drop_duplicates('SubjectId')
    b = b.merge(dm, on='SubjectId', how='left')
    b = b.dropna(subset=['mfi_delta', 'antigen', 'SubjectId', 'study_group']).copy()
    if 'visit_day' in b.columns and b['visit_day'].notna().any():
        e = endline_day if endline_day is not None else int(b['visit_day'].dropna().max())
        b = b[b['visit_day'] == e]
    else:
        e = endline_day
    b['arm'] = np.where(b['study_group'].astype(str) == '1', 'vaccine', 'placebo')
    b['mag'] = np.log10(np.clip(pd.to_numeric(b['mfi_delta'], errors='coerce'), min_mfi, None))
    b = b.dropna(subset=['mag'])
    b = b.drop_duplicates(['SubjectId', 'antigen'], keep='last')
    K = int(n_cat)
    b['y'] = b.groupby('antigen')['mag'].transform(
        lambda s: pd.qcut(s.rank(method='first'), K, labels=False)).astype(int)
    # keep antigens where BOTH arms span all K categories (item_time coverage requirement)
    def _cov(g):
        return all(g.loc[g.arm == a, 'y'].nunique() == K for a in ('vaccine', 'placebo'))
    kept = sorted(ag for ag, g in b.groupby('antigen') if _cov(g))
    dropped = sorted(set(b['antigen'].unique()) - set(kept))
    b = b[b['antigen'].isin(kept)].copy()
    b['pid_label'] = b['SubjectId'].astype(str)
    b['pid'] = pd.factorize(b['pid_label'])[0] + 1
    dp1 = pd.DataFrame({
        'group': (b['arm'] == 'vaccine').astype(int).to_numpy(),
        'group_label': np.where(b['arm'] == 'vaccine', 'vaccine', 'placebo'),   # display = arm; contrast keyed on numeric time
        'pid': b['pid'].to_numpy(), 'pid_label': b['pid_label'].to_numpy(),
        'fid': np.nan, 'f_label': b['arm'].to_numpy(), 'submission_date': pd.NaT,
        'treat': (b['arm'] == 'vaccine').astype(float).to_numpy(),
        'item_label': b['antigen'].astype(str).to_numpy(), 'y': b['y'].to_numpy(),
    })
    dp1['y_stan'] = dp1['y'] + 1                        # 1..K model input
    dp1['y_label'] = dp1['y'].astype(str)
    dp1['item_type'] = 'out-of-7'; dp1['item_type_id'] = 1
    antigens = sorted(dp1['item_label'].unique())
    dit = pd.DataFrame({'item_label': antigens})
    dit['item_type'] = 'out-of-7'; dit['item_type_id'] = 1; dit['cat_length'] = K
    dit['item_label_short'] = dit['item_label'].str.slice(0, 18)
    dit['construct'] = 'HIV BAMA IgG'; dit['construct_long'] = 'HIV-1 binding IgG (log10 MFI, ordinal)'
    dit['item_high_label'] = 'higher_is_better'
    dit['endpoint_measure'] = 'mean ordinal binding-IgG level (vaccine vs placebo, endline)'
    return {'dp': dp1, 'dit': dit, 'K': K, 'endline_day': e, 'kept': kept, 'dropped': dropped}


def read_data_immport_covid_neut(xlsx_path, study='SDY1764', group='age',
                                 min_start_dilution=4.0, timepoint='last'):
    """ImmuneSpace/ImmPort COVID-19 serum-neutralisation -> BETWEEN-GROUP PCM frame (§3.20).
    SARS-CoV-2 neutralising ID50 titre is the ordinal readout; the SUBPOPULATION is the group
    axis (MYCELIUM/CAVD-BAMA style: group A -> time 0 'Baseline', group B -> time 1 'Endline';
    each participant in one group). `group`: 'age' (pediatric <18 vs adult) or 'severity'
    (severe {ARDS, MIS-C} vs mild {non-MIS-C, convalescent}). Neut Value is log10(ID50), so
    titre=10**Value and k=round(log2(titre/min_start)) on the 2-fold ladder (single study-wide
    K). One measurement per participant (`timepoint`='last' takes the latest day). Endpoint
    rho = relative group-B-vs-A shift in mean log2 titre, higher_is_better. dict(dp, dit, K).
    """
    asy = pd.read_excel(xlsx_path, sheet_name='Assays')
    ev = pd.read_excel(xlsx_path, sheet_name='Events')
    par = pd.read_excel(xlsx_path, sheet_name='Participants')
    arm = pd.read_excel(xlsx_path, sheet_name='Arms')
    n = asy[(asy['Study ID'] == study)
            & asy['Assay Subtype'].astype(str).str.contains('neutral', case=False)].copy()
    if n.empty:
        raise ValueError(f'{study}: no neutralisation assay rows')
    evd = ev[['Event ID', 'Start', 'Participant ID']].rename(columns={'Start': 'day'})
    n = n.merge(evd, on=['Event ID', 'Participant ID'], how='left')
    n['val'] = pd.to_numeric(n['Value'], errors='coerce')
    n = n.dropna(subset=['val'])
    # one measurement per participant (latest / earliest day)
    n = n.sort_values(['Participant ID', 'day'])
    n = n.groupby('Participant ID', as_index=False).last() if timepoint == 'last'         else n.groupby('Participant ID', as_index=False).first()
    n['titre'] = 10.0 ** n['val']
    n['k'] = np.clip(np.rint(np.log2(n['titre'] / float(min_start_dilution))), 0, None).astype(int)
    K = int(n['k'].max()) + 1
    # subpopulation group from age / arm-severity
    p = par[par['Study ID'] == study][['Participant ID', 'Arm ID', 'Age']].merge(
        arm[['Arm ID', 'Description']], on='Arm ID', how='left')
    if group == 'age':
        p['grp'] = np.where(pd.to_numeric(p['Age'], errors='coerce') < 18, 'pediatric', 'adult')
        gA, gB = 'adult', 'pediatric'
    elif group == 'severity':
        dl = p['Description'].astype(str).str.lower()
        sev = dl.str.contains('ards') | dl.str.contains('mis-c') & ~dl.str.contains('without')
        p['grp'] = np.where(sev, 'severe', 'mild'); gA, gB = 'mild', 'severe'
    else:
        raise ValueError("group must be 'age' or 'severity'")
    n = n.merge(p[['Participant ID', 'grp']], on='Participant ID', how='left').dropna(subset=['grp'])
    n = n[n['grp'].isin([gA, gB])]
    n['pid_label'] = n['Participant ID'].astype(str)
    n['pid'] = pd.factorize(n['pid_label'])[0] + 1
    item = 'SARS-CoV-2 serum-neutralisation ID50 (reciprocal dilution)'
    dp1 = pd.DataFrame({
        'group': (n['grp'] == gB).astype(int).to_numpy(),
        'group_label': np.where(n['grp'] == gB, gB, gA),   # display = group name; contrast keyed on numeric time
        'pid': n['pid'].to_numpy(), 'pid_label': n['pid_label'].to_numpy(),
        'fid': np.nan, 'f_label': n['grp'].to_numpy(), 'submission_date': pd.NaT,
        'treat': (n['grp'] == gB).astype(float).to_numpy(),
        'item_label': item, 'y': n['k'].to_numpy(),
    })
    # ID50 is a CONTINUOUS neutralisation readout binned into 2-fold classes -> reciprocal serum
    # dilution labels 1:4, 1:8, ... (no below-detection floor; k0 is the assay's lowest ID50)
    titre_labels = _dilution_ladder_labels(K, min_start_dilution, continuous=True)
    dp1['y_stan'] = dp1['y'] + 1
    dp1['y_label'] = dp1['y'].map(lambda k: titre_labels[int(k)])
    dp1['item_type'] = 'out-of-7'; dp1['item_type_id'] = 1
    import json
    cat_labels = json.dumps(titre_labels)
    dit = pd.DataFrame({'item_label': [item]})
    dit['item_type'] = 'out-of-7'; dit['item_type_id'] = 1; dit['cat_length'] = K
    dit['item_label_short'] = np.nan                    # legend == facet strip == the item name
    dit['construct'] = item
    dit['construct_long'] = item                        # facet strip label == legend
    dit['item_high_label'] = 'higher_is_better'
    dit['endpoint_measure'] = 'titres among participants'
    dit['cat_labels'] = cat_labels
    return {'dp': dp1, 'dit': dit, 'K': K, 'group': group, 'groups': (gA, gB)}
