# Scheduling Estimation Methods Comparison

This folder contains a unified benchmarking script (`05_Methods_comparison.mlx`) designed to evaluate and compare the four quasi-LPV scheduling trajectory estimation methods supported by **CHRONOS**:

1. **Frozen** (`compute_scheduling_frozen`)
2. **LPV-RTI / Warm-Start** (`compute_scheduling_iterative_fast`)
3. **LPV-SQP / Iterative** (`compute_scheduling_iterative`)
4. **Taylor Recursive** (`compute_scheduling_recursive`)

## Overview
Choosing the right scheduling estimation method involves a trade-off between prediction accuracy and computational complexity. Using the nonlinear **Furuta pendulum** system as a baseline, this script simulates a tracking scenario and overlays the results of all four methods. 

By running this benchmark, you can directly analyze:
* **Closed-Loop Performance:** How well each method tracks the reference during nonlinear transients.
* **Computational Burden:** The online CPU time required per sample (compared against the sampling period $T_s$) and the number of solver iterations.

**Usage:** Simply open and run `05_Methods_comparison.mlx` in MATLAB to generate the comparison plots.