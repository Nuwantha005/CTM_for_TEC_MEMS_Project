# Radial TEC Preliminary Optimization (MATLAB + COMSOL)

This repository contains a compact thermal-network model, optimization workflows, and COMSOL integration scripts for a **radial multi-stage thermoelectric cooler (TEC)** used in MEMS/chip cooling studies.

The project is organized around:
- A MATLAB solver for the radial TEC resistor network.
- Multiple optimization strategies (local, global, multi-objective, physics-constrained).
- Parameter sweeps and diagnostic studies.
- Optional high-fidelity validation with COMSOL LiveLink.

## What This Project Solves

Given chip geometry, heat flux, coolant temperature, material properties, and TEC geometry parameters, the code estimates:
- Node temperatures across silicon and TEC stages.
- Heat flow in/out (`Q_in`, `Q_out`).
- Feasible design regions and optimized parameter sets.

A summary of historical optimization outcomes is documented in:
- `OPTIMIZATION_RESULTS_SUMMARY.md`

## Repository Structure

- `main.m`: top-level entry script that adds paths and lists available workflows.
- `src/`: core implementation.
  - `src/core/`: solver, geometry, material model, and optimizer classes.
  - `src/config/`: JSON defaults + centralized optimization variable configuration.
  - `src/comsol/`: COMSOL LiveLink wrapper and high-fidelity batch runner.
  - `src/utils/`: result logging/plotting helpers and config path utility.
  - `src/tests/`: MATLAB verification script.
- `scripts/`: runnable workflows.
  - `scripts/sweeps/`: parameter/stage/wedge sweep runners.
  - `scripts/optimization/`: optimization pipelines and post-analysis.
  - `scripts/comsol/`: COMSOL connection, validation, and verification scripts.
  - `scripts/analysis/`, `scripts/diagnosis/`: model analysis/debug studies.
  - `scripts/solidworks/`: SolidWorks parameter/design table assets.
- `data/`: material-property files and COMSOL/SolidWorks parameter artifacts.
- `notes/`: research notes, markdown/latex papers, and references.
- `output/`: generated run outputs (solver logs, plots, optimization results, COMSOL data).

## Core Model Components

Main classes in `src/core/`:
- `RadialTECSolver.m`: iterative solve loop, convergence tracking, result export.
- `ThermalNetwork.m`: system assembly/solve for the compact thermal network.
- `TECGeometry.m`: stage geometry and resistance/conductance helper calculations.
- `MaterialProperties.m`: temperature-dependent property interpolation from `data/material_props/*`.
- `TECOptimizer.m`: local optimization (`fmincon`/fallback) around solver objective.
- `DesignOptimizer.m`: staged design sweep + current optimization + ranking/plotting.
- `PhysicsConstrainedOptimizer.m`: optimization with feasibility/physics constraints.
- `StageSweeper.m`, `ParametricSweeper.m`: scripted sweep frameworks.

## Requirements

### MATLAB

Tested design expects MATLAB with:
- Base MATLAB (classdef, plotting, json functions)
- Optimization Toolbox (`fmincon`)
- Global Optimization Toolbox (`ga`, `particleswarm`, `gamultiobj`) for global/multi-objective scripts
- Parallel Computing Toolbox for `run_thickness_temp_optimization_parallel.m`

If some toolboxes are unavailable, many scripts still run in reduced mode (for example local fallbacks), but global/multi-objective workflows require their toolboxes.

### COMSOL (optional)

For COMSOL workflows:
- COMSOL Multiphysics with LiveLink for MATLAB
- COMSOL server reachable (commonly port `2036`)
- `.mph` model file matching the expected parameter names (see `src/comsol/COMSOLInterface.m` mapping)

## Quick Start

From MATLAB in the project root:

```matlab
run('main.m')
```

Then run one of the common scripts:

```matlab
run('scripts/sweeps/run_solver.m')
run('scripts/sweeps/run_stage_sweep.m')
run('scripts/optimization/run_optimization.m')
run('scripts/optimization/run_global_optimization.m')
```

A light validation run:

```matlab
run('scripts/sweeps/run_quick_test.m')
```

## Configuration

Primary config files:
- `src/config/default_params.json`: default solver geometry/material/boundary values.
- `src/config/test_params.json`: lighter test configuration.
- `src/config/optimization_variables.m`: central control for optimization variables, bounds, fixed integers (`N`, `M`), and boundary conditions.

Helper:
- `src/utils/get_config_path.m` resolves config paths robustly from any working directory.

## Typical Workflows

### 1) Run the thermal solver

```matlab
run('scripts/sweeps/run_solver.m')
```

Outputs are saved under timestamped folders in `output/solverRuns/`.

### 2) Sweep stage count / parameters

```matlab
run('scripts/sweeps/run_stage_sweep.m')
run('scripts/sweeps/run_sweep.m')
run('scripts/sweeps/run_wedge_sweep.m')
```

### 3) Run design/optimization studies

```matlab
run('scripts/optimization/run_design_optimization.m')
run('scripts/optimization/run_global_optimization.m')
run('scripts/optimization/run_multiobjective_optimization.m')
run('scripts/optimization/run_physics_constrained.m')
```

### 4) Validate with COMSOL (optional)

First ensure COMSOL server + model path are valid in your script/config.
Common entry points:

```matlab
run('scripts/comsol/test_comsol_connection.m')
run('scripts/comsol/run_comsol_test.m')
run('scripts/comsol/run_comsol_validation.m')
```

Or use the root-level scenario script:

```matlab
run('run_comsol_validation.m')
```

## Outputs

Generated artifacts are written to timestamped folders under `output/`, including:
- temperature profiles and convergence plots
- optimization logs/CSV summaries
- Pareto fronts and design ranking plots
- COMSOL candidate files and validation comparisons

## Known Caveats

- The repository is research-oriented and includes experimental scripts; not all scripts are equally maintained.
- Several analysis/diagnosis scripts reference legacy helper methods (for example `calculate_G`) that are not present in the current `TECGeometry` class; treat those scripts as investigational and update them before production use.
- Some COMSOL scripts include machine-specific absolute Windows paths; adjust model/server/LiveLink paths for your environment.
- The `output/` directory is large and historical; most workflows create new timestamped subfolders rather than overwriting older results.

## Verification

Basic solver verification script:

```matlab
run('src/tests/verify_solver.m')
```

## Notes and Documentation

Supporting theory and derivations are in:
- `notes/paper/Thermal_Network_For_Radial_TEC.tex`
- `notes/markdown/`

These notes describe the notation and modeling assumptions used by the codebase.
