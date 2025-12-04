# Energy System Optimization with Bernstein Polynomials

This repository contains Julia code for energy system optimization using Bernstein polynomial approximations. The code implements both continuous-time and discrete-time formulations for optimizing energy production, storage, and transmission in systems with thermal, hydropower, and wind generation.

## Overview

The optimization models minimize total system cost while satisfying energy demand across multiple scenarios. Key features:

- **Continuous model** (`run_continuous_model.jl`): Uses Bernstein polynomial approximation for continuous-time optimization
- **Discrete model** (`run_discrete_model.jl`): Traditional discrete-time optimization formulation
- Support for thermal, hydro, and wind power plants
- Multi-scenario stochastic optimization
- Reservoir storage and water value modeling
- Load shedding and power dumping capabilities

## Prerequisites

- **Julia** 1.6 or higher ([download here](https://julialang.org/downloads/))
- **CPLEX** optimization solver with valid license (academic or commercial)
  - Academic licenses available at [IBM Academic Initiative](https://www.ibm.com/academic/home)
  - Requires installation and configuration before running the code

## Installation

### 1. Install Julia

Download and install Julia from [julialang.org](https://julialang.org/downloads/). Verify installation:

```bash
julia --version
```

### 2. Install CPLEX

1. Download CPLEX from IBM (academic or commercial license)
2. Install CPLEX following IBM's instructions
3. Configure the CPLEX.jl Julia package to find your installation

For CPLEX.jl configuration help, see: https://github.com/jump-dev/CPLEX.jl

### 3. Install Julia Package Dependencies

Navigate to the repository directory and run:

```bash
cd vegardk-paper1
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This will install all required Julia packages listed in `Project.toml`.

## Directory Structure

```
vegardk-paper1/
├── input/                    # Input data files (required)
│   ├── gen.csv              # Generator specifications
│   ├── moduldata.xlsx       # Hydropower module data
│   ├── consumption.xlsx     # Load data
│   ├── DAY_AHEAD_wind.csv   # Wind forecast data
│   └── CAISO/NYISO data/    # Load time series
├── output/                   # Generated intermediate files (auto-created)
├── continuous_results/       # Continuous model outputs (auto-created)
├── discrete_results/         # Discrete model outputs (auto-created)
├── results/                  # Plot outputs (auto-created)
└── *.jl                     # Julia source files
```

Output directories are automatically created when you run the models.

## Usage

### Running the Discrete Model

```bash
julia --project=. run_discrete_model.jl
```

This runs the discrete-time optimization model and writes results to `discrete_results/`.

### Running the Continuous Model

```bash
julia --project=. run_continuous_model.jl
```

This runs the continuous-time optimization model with Bernstein polynomial approximation and writes results to `continuous_results/`.

### Key Parameters

To modify model parameters, edit the values in the main scripts:

**In `run_continuous_model.jl`:**
- `steps_per_hour = 4` - Time resolution (timesteps per hour)
- `input_data_scenarios = 10` - Number of scenarios to generate
- `simulation_scenarios = 9` - Number of scenarios to simulate
- `bernstein_degree = 4` - Degree of Bernstein polynomials
- `weights_calc_sampling_points = 1800` - Sampling resolution for weight calculation
- `res_sampling_points = 60` - Sampling resolution for results

**In `run_discrete_model.jl`:**
- `steps_per_hour = 4` - Time resolution
- `scenarios = 10` - Number of scenarios

## Output Files

Results are written to:

- **`discrete_results/results.xlsx`** - Discrete model results with sheets for:
  - Production by plant
  - First stage decisions
  - Transmission flows
  - Area-level results

- **`continuous_results/results.xlsx`** - Continuous model results with similar structure

- **`results/*.png`** - Plots and visualizations comparing models

- **`output/*.xlsx`** - Intermediate processed data files

## Code Structure

Main files:
- `run_continuous_model.jl` - Continuous model entry point
- `run_discrete_model.jl` - Discrete model entry point
- `continuous_model_simulator.jl` - Continuous model definition
- `discrete_model.jl` - Discrete model definition
- `process_input_data.jl` - Input data preprocessing
- `find_bernstein_weights.jl` - Bernstein polynomial calculations
- `plot_results.jl` - Result visualization
- `continuous_write_results.jl` - Output writing for continuous model
- `helper_functions.jl` - Utility functions
- `get_hydro_data.jl` - Hydropower data processing
- `elevation_matrix.jl` - Matrix operations for elevation/volume curves


## Citation

If you use this code in your research, please cite: https://ieeexplore.ieee.org/abstract/document/11050190
