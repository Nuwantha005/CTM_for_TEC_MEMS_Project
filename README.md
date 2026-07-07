# Radial Thermoelectric Cooler (TEC) Preliminary Optimization

This repository contains the Compact Thermal Model (CTM), optimization workflows, and COMSOL Multiphysics validation scripts for a **radial multi-stage thermoelectric cooler (TEC)** designed for MEMS and microchip cooling applications.

The project is structured around:
* **Analytical Thermal Solver**: A fast, iterative MATLAB-based solver for a radial multi-stage TEC thermal resistor network.
* **Optimization Framework**: Multiple optimization strategies (gradient-based local, global genetic algorithms, multi-objective Pareto optimization, and physics-constrained solvers).
* **High-Fidelity Validation**: Scripts to automate COMSOL Multiphysics simulations via LiveLink for MATLAB to validate compact model predictions.
* **SolidWorks CAD Integration**: Automated SolidWorks design tables and equations for parameter-driven geometry generation.

---

## Technical Overview

Given the chip size, target cooling heat flux, sink/coolant temperature, thermoelectric material properties (Seebeck coefficient, electrical resistivity, thermal conductivity), and TEC geometry parameters, this solver computes:
1. **Node Temperatures**: Complete temperature profile across the silicon chip and all TEC stages (junctions).
2. **Heat Flows**: Active heat pump rates ($Q_{\text{in}}$, $Q_{\text{out}}$) and cooling performance.
3. **Optimized Designs**: Parameters minimizing chip junction temperature or minimizing thickness/power while keeping the chip within safe limits.

Historical optimization runs and designs are summarized in:
* [OPTIMIZATION_RESULTS_SUMMARY.md](file:///run/media/nuwa/Work/Semester%207/ME4311%20-%20MicroNano%20Electro%20Mechanical%20Systems%20and%20Nanotechnology/Project/Preliminary%20Optimization/Algorithm/OPTIMIZATION_RESULTS_SUMMARY.md)

---

## Repository Structure

The cleaned directory structure is organized logically as follows:

```
├── main.m                      # Root setup script; initializes paths and lists workflows
├── README.md                   # Project documentation
├── OPTIMIZATION_RESULTS_SUMMARY.md # Summary of optimization results and target designs
├── data/                       # Reference datasets
│   ├── SW_equations.txt        # SolidWorks geometry equation exports
│   └── *.xlsx, *.csv           # SolidWorks design tables and template files
├── docs/                       # Technical paper & documentation
│   ├── Thermal_Network_For_Radial_TEC.tex # LaTeX documentation source
│   ├── Thermal_Network_For_Radial_TEC.pdf # Compiled paper PDF
│   └── images/                 # Image assets for the LaTeX document
├── src/                        # Core library implementation
│   ├── core/                   # Solver, optimizer, geometry, and material classes
│   ├── config/                 # Default parameters and optimization bounds (JSON/MATLAB)
│   ├── comsol/                 # COMSOL LiveLink wrapper and helper classes
│   ├── utils/                  # Result logging, file paths, and plotting helpers
│   └── tests/                  # Automated verification and regression test suite
└── scripts/                    # Runnable user workflows
    ├── sweeps/                 # Parametric and stage sweeps
    ├── optimization/           # Local, global, and multi-objective optimization pipelines
    ├── comsol/                 # COMSOL validation run automation (3-stage & 5-stage)
    ├── analysis/               # Diagnostic scripts for solver behaviors and physical limits
    ├── diagnosis/              # Thermal resistance breakdown and mismatch diagnostics
    └── solidworks/             # SolidWorks design table setup guides
```

---

## Core Classes (`src/core/`)

* **`RadialTECSolver.m`**: Handles the main iterative solve loop, convergence checks, and data management.
* **`ThermalNetwork.m`**: Assembles the system matrix representing the radial TEC equations and solves for junction temperatures.
* **`TECGeometry.m`**: Computes the geometry of each stage (inner/outer radii, volumes, cross-sectional areas) and thermoelectric leg resistances.
* **`MaterialProperties.m`**: Handles temperature-dependent thermoelectric properties (interpolated from reference data).
* **`TECOptimizer.m`**: Encapsulates optimization routines (gradient-based local optimization using `fmincon`).
* **`DesignOptimizer.m`**: Orchestrates design searches, current optimizations, and designs ranking.
* **`PhysicsConstrainedOptimizer.m`**: Runs optimizations constrained to physically realistic limits (e.g. non-negative temperatures, valid coefficient of performance).

---

## Getting Started

### Prerequisites
* **MATLAB** (R2021a or newer recommended).
* **Required Toolboxes**:
  * Optimization Toolbox (for `fmincon`)
  * Global Optimization Toolbox (for `ga`, `particleswarm`, and `gamultiobj`)
  * Parallel Computing Toolbox (optional; required only for parallel sweep scripts)

### Installation & Path Setup
Always launch the project by running `main.m` in the root folder to properly initialize the MATLAB paths:
```matlab
run('main.m')
```

This will set up all required paths and display a menu of available workflows.

---

## Typical Workflows

### 1. Run a Simple Simulation
To run the thermal network solver with the default design parameters:
```matlab
run('scripts/sweeps/run_solver.m')
```
*Results will be saved in a timestamped folder under `output/solverRuns/`.*

### 2. Run Parameter Sweeps
Explore the influence of individual geometric variables:
```matlab
run('scripts/sweeps/run_sweep.m')        % General geometry sweeps
run('scripts/sweeps/run_stage_sweep.m')  % Sweep number of stages (N)
run('scripts/sweeps/run_wedge_sweep.m')  % Sweep number of wedges (M)
```

### 3. Run Optimization Pipelines
Find the best design candidate for a target heat flux:
```matlab
run('scripts/optimization/run_optimization.m')             % Local gradient-based
run('scripts/optimization/run_global_optimization.m')      % Global genetic algorithm
run('scripts/optimization/run_multiobjective_optimization.m') % Multi-objective (Pareto)
```

### 4. High-Fidelity COMSOL Validation
To validate compact model predictions against 3D finite element simulations in COMSOL:
1. Open a terminal and start the COMSOL server:
   ```bash
   comsolmphserver -port 2036
   ```
2. Verify connection from MATLAB:
   ```matlab
   run('scripts/comsol/test_comsol_connection.m')
   ```
3. Run validation scripts:
   ```matlab
   run('scripts/comsol/run_comsol_validation.m')         % 3-stage validation
   run('scripts/comsol/run_comsol_validation_5stage.m')  % 5-stage validation
   ```

---

## Central Configuration

* **[default_params.json](file:///run/media/nuwa/Work/Semester%207/ME4311%20-%20MicroNano%20Electro%20Mechanical%20Systems%20and%20Nanotechnology/Project/Preliminary%20Optimization/Algorithm/src/config/default_params.json)**: Contains base dimensions, default boundary conditions, material characteristics, and solver tolerances.
* **[optimization_variables.m](file:///run/media/nuwa/Work/Semester%207/ME4311%20-%20MicroNano%20Electro%20Mechanical%20Systems%20and%20Nanotechnology/Project/Preliminary%20Optimization/Algorithm/src/config/optimization_variables.m)**: Configures optimization parameters, upper/lower bounds, integer constraints, and fixed bounds (e.g. stage bounds).

---

## Code Verification
Before committing or publishing changes, run the automated test suite to verify solver correctness:
```matlab
run('src/tests/verify_solver.m')
```
This script checks the solver initialization, executes a full test run, verifies vector dimensions, exports COMSOL properties, and performs a general energy balance sanity check.

---

## Documentation
The mathematical modeling, compact thermal network derivations, and notation details are fully documented in:
* **[Thermal_Network_For_Radial_TEC.tex](file:///run/media/nuwa/Work/Semester%207/ME4311%20-%20MicroNano%20Electro%20Mechanical%20Systems%20and%20Nanotechnology/Project/Preliminary%20Optimization/Algorithm/docs/Thermal_Network_For_Radial_TEC.tex)**: LaTeX document source with full derivations.
* **[Thermal_Network_For_Radial_TEC.pdf](file:///run/media/nuwa/Work/Semester%207/ME4311%20-%20MicroNano%20Electro%20Mechanical%20Systems%20and%20Nanotechnology/Project/Preliminary%20Optimization/Algorithm/docs/Thermal_Network_For_Radial_TEC.pdf)**: Compiled PDF of the paper.
* **[docs/images/](file:///run/media/nuwa/Work/Semester%207/ME4311%20-%20MicroNano%20Electro%20Mechanical%20Systems%20and%20Nanotechnology/Project/Preliminary%20Optimization/Algorithm/docs/images/)**: Diagnostic figures and geometry schematics used by the paper.

