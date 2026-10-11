# Verification notes — code audit, October 2026

## What was checked

- Archive integrity, full filename inventory and all 22 original `.m` source files reviewed at code level.
- Function name dependencies between the 22 `.m` files were mapped; all project-specific functions called in the main pipelines are present in the ZIP.
- The four benchmark exact-solution/forcing relations were spot-checked independently with symbolic algebra **at integer-order parameter settings**.
- Basic mathematical projection identities for shifted Legendre polynomials (degrees 0–2 on `[0,1]`) were spot-checked independently.
- These independent mathematical checks **do not validate the MATLAB implementation or reproduce thesis error tables**.

## Runtime verification — not completed

Neither MATLAB nor a compatible Octave interpreter was available in the review environment. MATLAB scripts **have not been executed**. Symbolic Math Toolbox behavior, plotted outputs and numerical error values are not verified yet. Runtime status: **NOT VERIFIED**.

## Corrections made to this organized copy

1. `examples/thesis/first_order_single_non_linear_complete_pipeline.m`: replaced `subs(u_exact,t,dom_ex)` with `subs(u_exact,x,dom_ex)` in the final plotting section. `t` was never defined in this script.
2. `examples/thesis/second_order_single_non_linear_complete_pipeline.m`: corrected only figure comment `N = 10` to `N = 8` to match `N_fixed = 8`.
3. `examples/thesis/second_order_coupled_system_non_linear_complete_pipeline.m`: corrected only figure comment `N = 10` to `N = 8` to match `N_fixed = 8`.

No other algorithms, matrix formulas or experiment parameters were changed. Preserve and compare the user's original ZIP if modifying science code further.

## Known concerns before calling it reproduced research

- The Riccati exact solution `tanh(x)` vanishes at x=0. The current plot formula `abs(exact-approx)/abs(exact)` has an undefined value at that point. Some coupled exact benchmark functions also vanish at endpoints. Use explicit zero handling for relative-error plots; do not silently replace with an invented value.
- In several routines, solutions are rounded using `vpa(...,8)` at each QLM iteration, which may affect high-accuracy results. Do not claim full precision convergence without testing.
- `ezplot` is a legacy MATLAB plotting function and may be unavailable/deprecated in future MATLAB releases. Consider `fplot` only after verifying identical plots.
- Some scripts pause inside loops, and shifted Legendre operators use nested symbolic sums and integrals; performance can be slow for larger basis sizes.
- Figure directories are created but original thesis scripts do **not** consistently export figures automatically; PNGs in `results/figures/` are archived outputs, not demonstrated regenerated results.
- Output figures from an extra nested graph archive include filenames relating to third-order problems that are not included in these four thesis pipelines; they were not included in the clean release draft.
- No precise MATLAB version/toolbox configuration or runtime reproduction log was provided.

## Release and attribution checklist

- [x] Author states this is his own MS thesis MATLAB implementation and requests public sharing.
- [ ] Confirm no third-party or confidential material is included and check any institutional/manuscript restrictions that may apply.
- [ ] Execute smoke tests in the local MATLAB installation.
- [ ] Execute each of the four thesis scripts and compare figures/error tables to the MS thesis.
- [ ] Fix singular relative-error plots without altering scientific claims.
- [ ] Decide on licensing only after ownership and authorization are established.
- [ ] Publish the curated repository with its verification status stated clearly; link from the CV once public.
