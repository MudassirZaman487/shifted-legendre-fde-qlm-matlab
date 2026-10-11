# Shifted Legendre Operational Matrix Methods for Nonlinear Fractional Differential Equations

**MATLAB research implementation for the MS Mathematics thesis of Mudassir Zaman (2026).**

**Thesis:** *A Legendre Polynomial-Based Operational Matrix Technique for Solving Nonlinear Fractional Differential Equations.*  
**University:** University of Education, Lahore (Attock Campus), Pakistan.  
**Supervisor:** Dr. Hammad Khalil. **Co-supervisor:** Dr. Abuzar Ghaffari.

>

## Overview

This research uses **shifted Legendre polynomial expansions and operational matrices** for fractional integration, Caputo differentiation, and multiplication by variable coefficients. Linear single-equation and coupled problems are reduced to finite-dimensional algebraic systems. Nonlinear problems are addressed by **quasi-linearization (QLM)**.

The four thesis computational examples cover:

1. First-order nonlinear Riccati-type fractional equation.
2. Second-order nonlinear single equation with variable coefficients and a lower-order derivative term.
3. First-order coupled nonlinear fractional system.
4. Second-order coupled nonlinear fractional system.

All original 22 MATLAB `.m` files are included in this repository.

## Requirements

- MATLAB with **Symbolic Math Toolbox** (`syms`, `symsum`, symbolic `int`, `vpa`, `subs` are used extensively).
- Graphical MATLAB environment for the original plotting scripts (`ezplot`, `figure`, `subplot`).
- These codes are run in MATLAB R2014b but you can use any other version.

## Getting started

Open MATLAB, set the current directory to the repository root, then run:

```matlab
setup_paths
run('examples/thesis/first_order_single_non_linear_complete_pipeline.m')
```

Other runs:

```matlab
run('examples/thesis/second_order_single_non_linear_complete_pipeline.m')
run('examples/thesis/first_order_coupled_system_non_linear_complete_pipeline.m')
run('examples/thesis/second_order_coupled_system_non_linear_complete_pipeline.m')
```

The scripts use symbolic computations, so they may be slow. Some contain deliberate `pause(1)` calls. These are legacy research scripts, not optimized fast demos. Do not run all examples at once until the first example has been tested locally.

For a basic **operator smoke test** :

```matlab
run('tests/test_basic_operators.m')
```

## Directory map

- `src/operators/` — shifted Legendre basis, spectral coefficients, fractional differentiation/integration matrices, variable-coefficient/product matrices.
- `src/linear_solvers/` — single/coupled, first/second-order matrix formulations.
- `src/qlm/` — nonlinear quasi-linearization routines.
- `examples/thesis/` — four full thesis numerical pipelines.
- `examples/matrix_demos/` — approximation, integration, differentiation, variable coefficient experiments.
- `results/figures/` — curated PNGs from the top-level graph directory of the submitted archive; existing visual outputs.
- `tests/` — basic mathematical smoke test for matrix operators.
