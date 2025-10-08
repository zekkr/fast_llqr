# Fast Quantile Regression Research

A comprehensive simulation study comparing Local Linear Quantile Regression (LLQR) and Time-Varying Coefficient Quantile Regression (TVCQR).

## Project Structure

- `scripts/` - Main execution scripts
- `R/` - R function libraries  
- `src/` - C++ and Fortran source code
- `data/` - Data and simulation results
- `results/` - Output figures and tables
- `docs/` - Documentation and papers

## Quick Start

1. Open `fast_qr_research.Rproj` in RStudio
2. Run `scripts/00_setup.R` to load all dependencies
3. Execute specific simulation scripts:
   - LLQR: `scripts/01_simulation_llqr.R`
   - TVCQR: `scripts/02_simulation_tvcqr.R`

## Methods

- **LLQR**: Local Linear Quantile Regression
- **TVCQR**: Time-Varying Coefficient Quantile Regression