# Migration to zarr Format

**Date:** 2026-05-07

## Summary

Updated bIRTistic Python code to use **zarr** instead of NetCDF for MCMC draws storage.

## Benchmark Results

Format comparison for xarray-compatible storage:

| Format   | Full Load  | Compute    | Overall Speedup |
| -------- | ---------- | ---------- | --------------- |
| **zarr** | **0.008s** | **0.020s** | **4.2x faster** |
| NetCDF   | 0.073s     | 0.048s     | baseline        |
| HDF5     | 0.080s     | 0.047s     | 0.9x            |

**✅ zarr is 8.9x faster for loading, 2.4x faster for computations**

## Files Updated

### Core Library
- ✅ `python/model_fitting.py`
  - Updated all 3 model fitting functions
  - Changed `.nc` → `.zarr` 
  - Changed `az.from_netcdf()` → `az.from_zarr()`
  - Changed `idata.to_netcdf()` → `idata.to_zarr()`

- ✅ `python/analysis.py`
  - Added zarr support alongside NetCDF
  - Supports both `.zarr` and `.nc` (backward compatible)

### Scripts
- ✅ `scripts-py/Colombia_analysis_ordered_logit_ADVI_vs_HMC.py`
- ✅ `scripts-py/Colombia_analysis_ordered_logit_ADVI_vs_HMC_FAST.py`

### Benchmark Scripts
- ✅ `scripts-py/benchmark_storage_formats.py` - Comprehensive comparison
- `scripts-py/benchmark_zarr_simple.py` - Original test (kept for reference)

## Migration Path

### For New Runs
- ✅ Automatic - new runs will create `.zarr` files

### For Existing Data
**Option 1:** Keep existing NetCDF files (backward compatible)
- `analysis.py` supports both `.zarr` and `.nc`
- Old analyses will continue to work

**Option 2:** Convert to zarr (recommended for performance)
```python
import arviz as az
import shutil
from pathlib import Path

# Convert a single file
nc_file = "path/to/file_draws.nc"
zarr_file = nc_file.replace('.nc', '.zarr')

idata = az.from_netcdf(nc_file)
if Path(zarr_file).exists():
    shutil.rmtree(zarr_file)
idata.to_zarr(zarr_file)
print(f"✓ Converted {nc_file} → {zarr_file}")
```

## Benefits

1. **9x faster loading** - Critical for iterative analysis
2. **2.4x faster compute** - Better for xarray operations
3. **Same file size** - No storage penalty
4. **Chunked storage** - Better for parallel operations
5. **xarray native** - Maintains full compatibility
6. **Data integrity** - Verified identical results

## Dependencies

- ✅ `zarr>=2.5,<3` - Added to pixi.toml
- ✅ `h5netcdf>=1.8.1,<2` - Added for HDF5 comparison

## Notes

- zarr stores as a directory (not a single file)
- Git should ignore `*.zarr/` directories
- zarr supports cloud storage (S3, GCS) for future scaling
