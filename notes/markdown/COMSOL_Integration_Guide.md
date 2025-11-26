# COMSOL Integration Guide

## Overview

This guide explains how to run high-fidelity COMSOL simulations to validate the preliminary optimization results.

## Files Created

1. **`src/comsol/COMSOLInterface.m`** - COMSOL LiveLink wrapper
   - Connects to COMSOL server
   - Loads and modifies models
   - Runs simulations and extracts results

2. **`src/comsol/HighFidelityRunner.m`** - Automated batch runner
   - Loads candidates from preliminary optimization
   - Runs COMSOL simulations on all candidates
   - Compares MATLAB predictions vs COMSOL results
   - Generates validation report

3. **`test_comsol_connection.m`** - Connection test script
4. **`run_overnight_comsol.m`** - Overnight batch run script

## Prerequisites

1. **COMSOL 6.3** with LiveLink for MATLAB license
2. **COMSOL Server** running (or auto-start enabled)
3. **Template COMSOL model** linked to SolidWorks

## Quick Start

### Step 1: Test COMSOL Connection

```matlab
run('test_comsol_connection.m')
```

This will:
- Check if LiveLink is installed
- Try to connect to COMSOL server
- Report any issues

### Step 2: Add LiveLink to MATLAB Path

If connection fails, add COMSOL LiveLink to your MATLAB path:

```matlab
addpath('C:\Program Files\COMSOL\COMSOL63\Multiphysics\mli')
```

Or run MATLAB through COMSOL:
```
comsol mphserver matlab
```

### Step 3: Run Overnight Simulations

1. Edit `run_overnight_comsol.m`:
   - Set `OPTIMIZATION_DIR` to your latest optimization
   - Set `COMSOL_MODEL` to your template model path
   - Adjust `MAX_CANDIDATES` based on simulation time

2. Run the script:
```matlab
run('run_overnight_comsol.m')
```

## Parameter Mapping

The system automatically maps MATLAB parameters to COMSOL LiveLink parameters:

| MATLAB Parameter | COMSOL Parameter | Units |
|-----------------|------------------|-------|
| `theta_deg` | `LL_theta` | deg |
| `t_TEC_um` | `LL_t_TEC` | µm |
| `k_r` | `LL_k_r` | - |
| `I_current_A` | `I0` | A |
| `q_flux_W_m2` | `q` | W/m² |
| `r_chip_mm` | `LL_r_chip` | mm |
| ... | ... | ... |

## Output Files

After running, you'll find in `output/comsol_validation/<timestamp>/`:

- `result_candidate_XX.mat` - Individual simulation results
- `validation_report.json` - Overall comparison report
- `temperature_profiles/*.png` - Temperature field plots
- `comparison_plots/*.png` - MATLAB vs COMSOL comparison

## Estimated Runtime

- **Per simulation**: ~5-15 minutes (depends on mesh complexity)
- **50 candidates**: ~4-12 hours
- **81 feasible designs (full run)**: ~7-20 hours

## Troubleshooting

### COMSOL Server Not Running
```
comsolmphserver -port 2036
```

### LiveLink Not Found
Ensure `mli` folder is in MATLAB path:
```matlab
which mphstart
% Should return path to mphstart.m
```

### Model Loading Fails
- Check model path is correct
- Ensure model is compatible with COMSOL 6.3
- Verify SolidWorks LiveLink is working

### Memory Issues
- Reduce `MAX_CANDIDATES`
- Increase system RAM
- Use mesh refinement carefully

## Best Results Summary

From preliminary optimization (2025-11-25_21-50-48):

| Rank | N | θ (deg) | t_TEC (µm) | k_r | I_opt (mA) | T_max (°C) |
|------|---|---------|------------|-----|------------|------------|
| 1 | 3 | 30 | 200 | 1.2 | 96.9 | 81.6 |
| 2 | 3 | 25 | 200 | 1.2 | 90.4 | 84.7 |
| 3 | 3 | 30 | 250 | 1.2 | 111.9 | 85.6 |

**Key finding**: N=3 stages is optimal, achieving T_max = 81.6°C (below 100°C target)
