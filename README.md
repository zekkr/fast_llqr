## Paper U11 implementation and reproduction

The formal U11 entry is `R/u11_functions.R`; package synchronization, fixed-source
simulation commands, archived summaries and environment evidence are documented
in [reproduction/README.md](reproduction/README.md). See
[release status](reproduction/RELEASE_STATUS.md) before using a release URL.

# Fast Local Linear Quantile Regression Research

A comprehensive simulation study comparing Local Linear Quantile Regression (LLQR) and Time-Varying Coefficient Quantile Regression (TVCQR).

## Project Structure

- `scripts/` - Main execution and setup scripts
- `R/` - Core function libraries
- `src/` - C++ and Fortran source code
- `data/` - Simulation results and data storage
- `results/` - Output figures and tables
- `docs/` - Documentation and papers

## Quick Start

1. **Open Project**: Directly open `fast_llqr.Rproj` to enter RStudio for simulation. All paths are relative, so there is no need to change any addresses—simply open the `.Rproj` file.

2. **Setup Environment**: Run `source("scripts/setup")` to load all required dependencies and packages.

3. **Core Functions**:
   - LLQR solver: Located in `R/llqr_functions.R`
   - TVCQR solver: Located in `R/tvcqr_functions.R`
   - Additional helper functions for simulation are also stored in the `R/` folder.

## Running Simulations

To execute the main simulation study:

1. Open `reproduce_results.R`
2. Set the desired parameters in the `sim_config` section
3. Run the script to start simulations

**Output Locations**:
- LLQR results will be saved in `data/llqr_results/`
- TVCQR results will be saved in `data/tvcqr_results/`

## Results Analysis

Use `run_performance_analysis.R` to analyze generated simulation results. This script outputs:

- **Running time statistics** for different algorithms
- **Average relative bias** of different algorithms (compared to results from the `quantreg` package without acceleration)
- Various performance metrics and comparative statistics

## Additional Scripts

- `scripts/simulation_llqr.R` - Dedicated LLQR simulation script
- `scripts/simulation_tvcqr.R` - Dedicated TVCQR simulation script

## Notes

- All file paths are relative to the project root
- The project is configured to run directly from the Rproj file
- Simulation configurations can be adjusted in `reproduce_results.R` to test different scenarios