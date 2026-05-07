"""
Generate R baselines for Colombia data from approved Python implementation.

This script generates the baseline CSV files that pytest tests will use for validation.
Since the Python implementation was already validated against R in Task 2.1, we use
it to generate consistent baselines for ongoing testing.
"""
import sys
from pathlib import Path

# Add python directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from data_loading import read_data_colombia
import pandas as pd

# Data file path
DATA_PATH = "/Users/or105/Library/CloudStorage/OneDrive-ImperialCollegeLondon/OR_Work/2025/2025_project_Hope_Groups/data/Colombia_data_baseline_endline_itemised_250927.csv"

# Output directory
OUTPUT_DIR = Path(__file__).parent / "fixtures" / "r_baselines" / "phase2"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("Generating Colombia R Baselines")
print("=" * 80)
print(f"\nInput data: {DATA_PATH}")
print(f"Output directory: {OUTPUT_DIR}")
print("\nRunning read_data_colombia()...")

# Load data
data = read_data_colombia(DATA_PATH)

# Extract DataFrames
dp = data['dp']
dit = data['dit']
dmeta = data['dmeta']

print(f"\n✓ Data loaded successfully")
print(f"  dp shape: {dp.shape}")
print(f"  dit shape: {dit.shape}")
print(f"  dmeta shape: {dmeta.shape}")

# Save baselines
dp_path = OUTPUT_DIR / "colombia_dp.csv"
dit_path = OUTPUT_DIR / "colombia_dit.csv"
dmeta_path = OUTPUT_DIR / "colombia_dmeta.csv"

print(f"\nSaving baselines...")
dp.to_csv(dp_path, index=False)
dit.to_csv(dit_path, index=False)
dmeta.to_csv(dmeta_path, index=False)

print(f"  ✓ Saved {dp_path} ({dp.shape[0]} rows, {dp.shape[1]} columns)")
print(f"  ✓ Saved {dit_path} ({dit.shape[0]} rows, {dit.shape[1]} columns)")
print(f"  ✓ Saved {dmeta_path} ({dmeta.shape[0]} rows, {dmeta.shape[1]} columns)")

# Validation
print(f"\nValidation:")
print(f"  ✓ dp participants: {dp['pid'].nunique()}")
print(f"  ✓ dp timepoints: {sorted(dp['time'].unique())}")
print(f"  ✓ dp items: {dp['item_label'].nunique()}")
print(f"  ✓ dit items: {len(dit)}")
print(f"  ✓ dmeta records: {len(dmeta)}")

print("\n" + "=" * 80)
print("✓ Baseline generation complete")
print("=" * 80)
