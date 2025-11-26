# Proposal for Radial TEC Solver Enhancements

## 1. Codebase Understanding

### 1.1 Theory & Model
The project implements a **Finite Difference Method (FDM)** solver for a **Radial Multi-Stage Thermoelectric Cooler (TEC)**.
- **Geometry**: The TEC is modeled as concentric rings (stages). The geometry is defined by parameters like wedge angle, number of stages ($N$), and radial expansion factors.
- **Physics**: The model solves coupled thermal and electrical equations.
    - **Thermal**: Heat conduction, Peltier cooling/heating, Joule heating, and heat generation.
    - **Electrical**: Current flow through legs and interconnects.
- **System of Equations**: The problem is discretized into nodes (Silicon chip nodes and TEC junction nodes). This results in a linear system $[M][T] = [B]$ for a fixed temperature distribution.
- **Non-linearity**: Since material properties ($k, \rho, S$) depend on temperature, the system is non-linear and requires an iterative solution scheme.

### 1.2 Current Implementation
- **`RadialTECSolver.m`**: The main entry point. Currently, it performs a **single linear solve**, which is insufficient for the non-linear problem.
- **`ThermalNetwork.m`**: Assembles the system matrix $M$ and RHS vector $B$ based on a given temperature distribution `T_current`.
- **`TECGeometry.m`**: Handles the complex radial geometry calculations (geometric factors $G$, resistances $R$).
- **`MaterialProperties.m`**: Provides temperature-dependent material properties.

## 2. Proposed Enhancements

We will wrap the core solver with two new layers: a **Parametric Sweeper** and an **Optimizer**. Before that, the core solver needs to be robust.

### 2.1 Step 1: Robust Iterative Solver (Refactoring `RadialTECSolver`)
The current `run` method only solves once. We need to implement an iterative loop (Fixed-Point Iteration) to handle the temperature dependence.

**Plan:**
- Modify `RadialTECSolver` to include a `solve_iterative()` method.
- **Algorithm**:
    1. Initialize $T_{guess}$.
    2. Loop `k = 1` to `max_iter`:
        a. Update material properties using $T_{guess}$.
        b. Assemble and solve $[M][T_{new}] = [B]$.
        c. Check convergence: if $||T_{new} - T_{guess}|| < tol$, break.
        d. Update $T_{guess} = (1-\alpha)T_{guess} + \alpha T_{new}$ (Relaxation may be needed).
    3. Return converged $T$.

### 2.2 Step 2: Parametric Sweeper
A tool to systematically vary one or two parameters and observe the effect on performance (e.g., $T_{cold}$ vs. Current $I$).

**Plan:**
- Create `src/core/ParametricSweeper.m`.
- **Features**:
    - Input: `parameter_path` (e.g., "operating_conditions.I_current_A"), `range` (e.g., `0:0.1:5`).
    - Loop through the range, update the `solver.Config`, run `solve_iterative()`, and store results.
    - **Output**: A structure or table containing the sweep parameter and corresponding results ($T_{cold}$, $Q_c$, $COP$, etc.).
    - **Visualization**: Built-in plotting method to visualize the sweep.

### 2.3 Step 3: Optimizer
A tool to automatically find the optimal parameters to minimize/maximize an objective.

**Plan:**
- Create `src/core/TECOptimizer.m`.
- **Features**:
    - **Objective Function**: A wrapper that takes a vector of design variables $x$, updates the solver config, runs the simulation, and returns a scalar cost (e.g., $T_{cold}$).
    - **Variables**: Current ($I$), Geometry ($N_{stages}$, $L_1$, ratios).
    - **Algorithm**:
        - Use MATLAB's `fminsearch` (Nelder-Mead) for unconstrained/simple constrained problems.
        - Use `fmincon` if the Optimization Toolbox is available (better for constraints).
        - Alternatively, implement a simple Gradient Descent or Genetic Algorithm if toolboxes are restricted.
    - **Constraints**: Handle physical constraints (e.g., $T < T_{max}$, geometry limits) via penalty functions if using `fminsearch`.

## 3. Implementation Roadmap

1.  **Update `RadialTECSolver.m`**: Implement the convergence loop.
2.  **Create `ParametricSweeper.m`**: Implement single-parameter sweeping.
3.  **Create `TECOptimizer.m`**: Implement optimization wrapper.
4.  **Create `run_optimization.m`**: A script to demonstrate the usage.

## 4. File Structure Changes

```
src/
    core/
        RadialTECSolver.m   (Modified)
        ParametricSweeper.m (New)
        TECOptimizer.m      (New)
    ...
```
