# Tunable asymmetric swimming in biflagellate microswimmers

MATLAB codebase for analysing robot trajectory curvature data and simulating several swimmer models with tunable left/right actuation asymmetry.

## Overview

This repository contains:

- **Experimental data analysis** from supplied curvature data.
- **Minimal and intermediate ODE models** for qualitative phototaxis/curvature behaviour.
- **Two-link hydrodynamic model** (numerics + symbolic derivation support).
- **Shared geometry/utility functions** (circle fitting, curvature extraction, colormaps).

## Repository Structure

- `common/` — shared helper functions (`circleFit`, `curvature_from_traj`, colormap utility).
- `data_analysis/` — fit and plot curvature scaling law with frequency.
- `minimal_model/` — compact swimmer ODE model and figure scripts.
- `intermediate_model/` — intermediate asymmetry model and resonance.
- `two_link_model/numerics/` — full two-link numerical simulation, phototaxis controller studies, and video rendering.
- `two_link_model/symbolic/` — symbolic derivation scripts for linear systems/control matrices.

## Requirements

- MATLAB (R2021a+ recommended).
- Toolboxes used by scripts:
  - **Curve Fitting Toolbox** (`fit`, `fittype`) for curvature-law fitting scripts.
  - **Symbolic Math Toolbox** for `two_link_model/symbolic/` scripts.

## Quick Start

From MATLAB, set the current folder to the repository root, then run scripts by folder.

### 1 Data Analysis

`data_analysis/main_supplied_curvs.m`
- Uses supplied table `data_analysis/data/curvature_values_0731-1.csv`.
- Fits scaling law:

$$
\kappa = a\,\frac{k_L/k_R - 1}{k_L/k_R + 1}
$$

- Plots mean curvature with error bars across repeated multipliers.
- Exports `fit_to_robot_data.png` and `fit_to_robot_data.pdf`.

### 2 Minimal Model

`minimal_model/main.m`
- Solves a 3-state ODE for position and heading.
- Plots swimmer trajectory for chosen forcing/asymmetry functions.

`minimal_model/fig_curvavg_simplemodel.m`
- Sweeps frequency ratio and forcing functions.
- Produces curvature-scaling figure.
- Exports `fig_phototaxis_fit.pdf`.

`minimal_model/fig_phototaxis_simplemodel.m`
- Generates trajectory panels for different forcing profiles and frequency ratios.
- Exports `fig_phototaxis_example.pdf`.

### 3 Intermediate Model

`intermediate_model/main.m`
- Simulates an intermediate model with selectable asymmetry type:
  - `"amplitude"`
  - `"frequency"`
  - `"phase"`
- Plots the resulting trajectory.

`intermediate_model/resonant_freq_asymmetry.m`
- Computes and plots the resonance integral against frequency ratio.
- Exports `resonant_freq_asymmetry.pdf`.

### 4 Two-Link Model (Numerics)

`two_link_model/numerics/main.m`
- Runs the two-link swimmer simulation with prescribed gait functions.
- Calls `make_video(...)` to write an animation sequence (default output path `./ani`).

`two_link_model/numerics/phototaxis/main.m`
- Runs cycle-by-cycle phototaxis control simulation with a point light source.
- Produces trajectory/angle plots and optional animation output.

## Typical Workflow

1. Run `data_analysis/main_supplied_curvs.m` to recover the fitted experimental scaling law.
2. Run `minimal_model/fig_curvavg_simplemodel.m` to compare model scaling trends.
3. Use `intermediate_model/main.m` and/or `intermediate_model/resonant_freq_asymmetry.m` for asymmetry and resonance exploration.
4. Use `two_link_model/numerics/main.m` or `two_link_model/numerics/phototaxis/main.m` for higher-fidelity simulations.

## Notes

- Most scripts are written as top-level MATLAB scripts (not functions). Run them from their containing folder, or ensure relative paths are valid.
- Many scripts call `addpath(genpath('../common'))` or similar; keeping folder structure unchanged is recommended.
- Units are mixed by source (e.g., trajectory parser converts CSV positions from mm to m).

## Troubleshooting

- **`fit`/`fittype` not found**: install/enable Curve Fitting Toolbox.
- **Symbolic functions unavailable**: install Symbolic Math Toolbox for `two_link_model/symbolic`.
- **Missing files in `./data`**: verify script is launched from the expected folder and data files are present.
