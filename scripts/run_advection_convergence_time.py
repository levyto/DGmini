#!/usr/bin/env python3
"""
Run a 1D linear-advection *time-convergence* study for DGmini.

This script measures temporal accuracy on a fixed spatial discretization:
  1. Reads a base YAML config (typically examples/advection.yaml).
  2. Fixes the mesh size N for all runs.
  3. For each requested time integrator, uses an integrator-specific polynomial
     order p and starting CFL if provided, otherwise falls back to CLI defaults.
  4. Runs a sequence of CFL values CFL, CFL/2, CFL/4, ... for each integrator
     (equivalent to halving dt on a fixed mesh).
  5. Runs one additional *reference* solve for each integrator on the same mesh
     and with the same p, using a much smaller reference CFL chosen at the top
     of this file.
  6. Computes the L2 error of each run against the corresponding reference run.
  7. Prints convergence tables and writes summary CSV + log-log plots.
"""

from __future__ import annotations

import argparse, math
from pathlib import Path

from utils_convergence import (
  average_dt, 
  compute_l2_error_against_reference, 
  infer_repo_root_from_config, 
  load_yaml, 
  make_plot, 
  pairwise_eoc_from_x, 
  parse_solver_stdout, 
  prepare_results_dir, 
  print_time_table, 
  run_solver, 
  validate_time_settings, 
  write_case_config, 
  write_summary_csv
)

# ------------------------------------------------------------------------------
# User settings for the convergence study. Adjust as needed.
# ------------------------------------------------------------------------------

# List of integrators to test. This can be overridden by CLI args.
INTEGRATORS = ["forward_euler", "rk2", "rk3_ssp", "rk4"]

# Starting CFL for each integrator. Each next level uses CFL/2.
CFL_BY_INTEGRATOR = {"forward_euler": 0.10, "rk2": 0.60, "rk3_ssp": 0.80, "rk4": 0.80}

# Polynomial order p for each integrator. This is used for both the reference 
# and test runs of that integrator.
P_BY_INTEGRATOR = {"forward_euler": 1, "rk2": 2, "rk3_ssp": 3, "rk4": 4}

# Reference CFL used for the fine-time reference solution of each integrator.
# This should be smaller than the finest tested CFL for that integrator.
REF_CFL_BY_INTEGRATOR = {k: v / (2 ** 6) for k, v in CFL_BY_INTEGRATOR.items()}

# Base number of elements. The script refines this as N, 2N, 4N, ...
N_ELEM = 10

# Initial condition used both in YAML and in the exact solution after one period
INITIAL_CONDITION_EXPRESSION = "sin(2*pi*x)"

# ------------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
  """
  Parse command-line arguments for the convergence study.
  """
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--exe",         type=Path,  default=Path("./build/dgmini"),            help="Path to DGmini executable.")
  parser.add_argument("--config",      type=Path,  default=Path("./examples/advection.yaml"), help="Base YAML config.")
  parser.add_argument("--levels",      type=int,   default=4,      help="Number of temporal refinement levels.")
  parser.add_argument("--integrators", nargs="+",  default=INTEGRATORS, choices=INTEGRATORS,  help="Time integrators to test.")
  parser.add_argument("--n-elements",  type=int,   default=N_ELEM, help="Fixed number of elements for all runs.")
  parser.add_argument("--p-order",     type=int,   default=5,      help="Fallback polynomial order p if an integrator is not listed in P_BY_INTEGRATOR.")
  parser.add_argument("--base-cfl",    type=float, default=0.4,    help="Fallback starting CFL if an integrator is not listed in CFL_BY_INTEGRATOR.")
  parser.add_argument("--results-dir", type=Path,  default=Path("./output/advection_convergence_time"), help="Directory where generated configs, outputs, tables, and plots are stored.")
  return parser.parse_args()

def validate_base_config(cfg: dict) -> None:
  """
  Validate that the base YAML config has settings compatible with this 
  convergence study.
  """
  if cfg["pde"]["type"] != "linear_advection": 
    raise ValueError(f"This script expects pde.type = 'linear_advection', got {cfg['pde']['type']!r}.")
  if cfg["time_step_controller"]["type"] != "cfl": 
    raise ValueError("This script is intended for CFL-driven runs so that dt scales with h.")
  speed = float(cfg["pde"].get("advection_speed", 0.0))
  if abs(speed) <= 0.0: 
    raise ValueError("Advection speed must be nonzero for one-rotation test.")

def compute_final_time(cfg: dict) -> float:
  """
  Compute the final time corresponding to one full rotation of the wave around 
  the domain, based on the mesh and advection speed settings in the config. This 
  is the time at which we will measure the error against the reference solution.
  """
  return (float(cfg["mesh"]["x_right"]) - float(cfg["mesh"]["x_left"])) / abs(float(cfg["pde"]["advection_speed"]))

def get_integrator_settings(integrator: str, args: argparse.Namespace) -> tuple[float, int, float]:
  """
  Get the base CFL, polynomial order, and reference CFL for a given integrator,
  using the predefined dictionaries or falling back to CLI args if not listed.
  """
  base_cfl = float(CFL_BY_INTEGRATOR.get(integrator, args.base_cfl))
  p_order  = int(P_BY_INTEGRATOR.get(integrator, args.p_order))
  ref_cfl  = float(REF_CFL_BY_INTEGRATOR.get(integrator, base_cfl / (2 ** 6)))
  return base_cfl, p_order, ref_cfl

def main() -> int:
  """
  Main function to run the convergence study. It generates YAML configs for each
  case, runs the solver, computes errors against a fine-time reference solution,
  and produces tables and plots of the results.
  """
  args        = parse_args()
  exe         = args.exe.resolve() 
  config_path = args.config.resolve()
  results_dir = args.results_dir.resolve()
  repo_root   = infer_repo_root_from_config(config_path)

  if not exe.exists(): 
    raise FileNotFoundError(f"Executable not found: {exe}")
  if not config_path.exists(): 
    raise FileNotFoundError(f"Config not found: {config_path}")
  if args.n_elements <= 0: 
    raise ValueError("--n-elements must be positive.")
  if args.p_order < 0: 
    raise ValueError("--p-order must be nonnegative.")
  if args.levels < 2: 
    raise ValueError("--levels must be at least 2 to measure an observed order.")
  if args.base_cfl <= 0.0: 
    raise ValueError("--base-cfl must be positive.")

  validate_time_settings(args.integrators, CFL_BY_INTEGRATOR, P_BY_INTEGRATOR, REF_CFL_BY_INTEGRATOR)
  base_cfg                 = load_yaml(config_path)
  validate_base_config(base_cfg)
  final_time               = compute_final_time(base_cfg)
  configs_dir, outputs_dir = prepare_results_dir(results_dir)

  all_results = {}

  summary_rows = []
  for integrator in args.integrators:

    rows = []

    base_cfl_i, p_order_i, ref_cfl_i = get_integrator_settings(integrator, args)

    # Reference run
    finest_test_cfl = float(base_cfl_i / (2 ** (args.levels - 1)))

    if ref_cfl_i >= finest_test_cfl: 
      raise ValueError(f"For integrator {integrator!r}, REF_CFL_BY_INTEGRATOR = {ref_cfl_i} must be smaller than the finest tested CFL = {finest_test_cfl}.")

    ref_case_name  = f"{integrator}_p{p_order_i}_N{args.n_elements:04d}_reference";
    ref_output_dir = outputs_dir / ref_case_name
    ref_config     = configs_dir / f"{ref_case_name}.yaml"

    write_case_config(
      base_cfg=base_cfg, 
      integrator=integrator, 
      n_elements=args.n_elements,
      p_order=p_order_i, 
      cfl_number=ref_cfl_i, 
      final_time=final_time, 
      initial_condition_expression=INITIAL_CONDITION_EXPRESSION, 
      output_dir=ref_output_dir, 
      config_path=ref_config
    )
    
    ref_stdout, _ = run_solver(exe=exe, config_path=ref_config, cwd=repo_root)
    ref_n_steps, ref_t_final_reported = parse_solver_stdout(ref_stdout)
    ref_file      = ref_output_dir / "solution_final.dat"

    if not ref_file.exists(): 
      raise FileNotFoundError(f"Expected reference output file not found: {ref_file}")

    ref_dt_avg = average_dt(final_time, ref_n_steps, ref_t_final_reported)

    # Other runs
    for level in range(args.levels):

      cfl_number = float(base_cfl_i / (2 ** level))

      case_name       = f"{integrator}_p{p_order_i}_N{args.n_elements:04d}_cfl{level}"
      case_output_dir = outputs_dir / case_name
      case_config     = configs_dir / f"{case_name}.yaml"

      write_case_config(
        base_cfg=base_cfg, 
        integrator=integrator, 
        n_elements=args.n_elements, 
        p_order=p_order_i, 
        cfl_number=cfl_number, 
        final_time=final_time, 
        initial_condition_expression=INITIAL_CONDITION_EXPRESSION, 
        output_dir=case_output_dir, 
        config_path=case_config
      )

      stdout, _                 = run_solver(exe=exe, config_path=case_config, cwd=repo_root)
      n_steps, t_final_reported = parse_solver_stdout(stdout)
      final_file                = case_output_dir / "solution_final.dat"

      if not final_file.exists(): 
        raise FileNotFoundError(f"Expected final output file not found: {final_file}")

      error_l2, h, n_elements_in_file = compute_l2_error_against_reference(final_file, ref_file)

      if n_elements_in_file != args.n_elements: 
        raise RuntimeError(f"Mismatch between requested N={args.n_elements} and file Ne={n_elements_in_file} in {final_file}")

      rows.append(
        { "integrator": integrator, 
          "N": float(args.n_elements), 
          "p": float(p_order_i), 
          "h": h, 
          "cfl": cfl_number, 
          "dt_avg": average_dt(final_time, n_steps, t_final_reported), 
          "n_steps": float(n_steps) if n_steps is not None else math.nan, 
          "ref_cfl": ref_cfl_i, 
          "ref_dt_avg": ref_dt_avg, 
          "ref_n_steps": float(ref_n_steps) if ref_n_steps is not None else math.nan, 
          "error_l2": error_l2, 
          "t_final": float(t_final_reported) if t_final_reported is not None else math.nan
        }
      )
    # end for level

    eoc_t = pairwise_eoc_from_x([r["dt_avg"] for r in rows], [r["error_l2"] for r in rows])

    for row, eoc_i in zip(rows, eoc_t): 
      row["eoc_t"] = eoc_i
      summary_rows.append(row.copy())

    all_results[integrator] = rows;
    
    print_time_table("DGmini advection temporal convergence", rows)

  # end for integrator

  write_summary_csv(
    results_dir / "summary.csv", 
    summary_rows, 
    ["integrator", "p", "N", "cfl", "dt_avg", "n_steps", "h", "ref_cfl", "ref_dt_avg", "ref_n_steps", "error_l2", "eoc_t", "t_final"]
  )
  
  make_plot(
    results_dir / "error_vs_dt.png", 
    all_results, 
    x_builder=lambda rows: [r["dt_avg"] for r in rows], 
    y_key="error_l2", 
    xlabel="dt_avg", 
    ylabel=r"$L^2$ error vs reference", 
    title="DGmini advection temporal convergence after one period", 
    label_builder=lambda group, rows, slope: f"{group}, p={int(rows[0]['p'])}, N={int(rows[0]['N'])}, ref CFL={rows[0]['ref_cfl']:g} (slope ~ {slope:.2f})"
  )
  
  make_plot(
    results_dir / "error_vs_nsteps.png", 
    all_results, 
    x_builder=lambda rows: [r["n_steps"] for r in rows], 
    y_key="error_l2", 
    xlabel="number of time steps", 
    ylabel=r"$L^2$ error vs reference", 
    title="DGmini advection temporal convergence after one period", 
    label_builder=lambda group, rows, slope: f"{group}, p={int(rows[0]['p'])}, N={int(rows[0]['N'])}, ref CFL={rows[0]['ref_cfl']:g} (slope ~ {slope:.2f})"
  )
  
  print()
  print(f"Results written to: {results_dir}")
  print("  - summary.csv")
  print("  - error_vs_dt.png")
  print("  - error_vs_nsteps.png")
  print(f"  - generated YAML configs in {configs_dir}")
  print(f"  - solver outputs in {outputs_dir}")
  
  return 0

if __name__ == "__main__": 
  raise SystemExit(main())
