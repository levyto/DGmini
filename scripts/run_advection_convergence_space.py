#!/usr/bin/env python3
"""
Run a 1D linear-advection *space-convergence* study for DGmini.

This script measures spatial accuracy:
  1. Reads a base YAML config (typically examples/advection.yaml).
  2. Loops over the polynomial orders listed in P_ORDER.
  3. For each p, selects the associated time integrator through INTEGRATOR_BY_P
     and the corresponding CFL through CFL_BY_INTEGRATOR.
  4. Runs a sequence of meshes N, 2N, 4N, ... for each (p, integrator) pair.
  5. Computes the L2 error against the exact periodic solution after one full revolution.
  6. Prints convergence tables and writes summary CSV + log-log plots.

Notes:
  - To isolate spatial accuracy, choose the CFL small enough that temporal error
    stays below spatial error on the tested meshes.
  - The exact solution after one period is the initial condition again.
  - The initial condition expression is injected into the YAML from
    INITIAL_CONDITION_EXPRESSION below and also used for the exact solution.
"""

from __future__ import annotations

import argparse, math
from pathlib import Path

from utils_convergence import (
  average_dt, 
  compute_l2_error_against_exact, 
  infer_repo_root_from_config, 
  load_yaml, 
  make_advection_exact_solution, 
  make_plot, pairwise_eoc_from_x, 
  parse_solver_stdout, 
  prepare_results_dir, 
  print_space_table, 
  run_solver, 
  validate_space_settings, 
  write_case_config, 
  write_summary_csv
)

# ------------------------------------------------------------------------------
# User settings for the convergence study. Adjust as needed.
# ------------------------------------------------------------------------------

# Polynomial orders to test
P_ORDER = [1, 2, 3, 4]

# Time integrator assigned to each polynomial order
INTEGRATOR_BY_P = {
    1: "rk3_ssp",
    2: "rk3_ssp",
    3: "rk3_ssp",
    4: "rk4",
}

# CFL associated with each time integrator
CFL_BY_INTEGRATOR = {
    "forward_euler": 0.05,
    "rk2": 0.10,
    "rk3_ssp": 0.10,
    "rk4": 0.10,
}

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
  parser.add_argument("--exe",         type=Path, default=Path("./build/dgmini"),            help="Path to DGmini executable.")
  parser.add_argument("--config",      type=Path, default=Path("./examples/advection.yaml"), help="Base YAML config.")
  parser.add_argument("--levels",      type=int,  default=4,      help="Number of spatial refinement levels.")
  parser.add_argument("--n-elements",  type=int,  default=N_ELEM, help="Base number of elements. The study uses N, 2N, 4N, ...")
  parser.add_argument("--results-dir", type=Path, default=Path("./output/advection_convergence_space"), help="Directory where generated configs, outputs, tables, and plots are stored.")
  return parser.parse_args()

def validate_base_config(cfg: dict) -> None:
  """
  Validate that the base config is compatible with the intended advection 
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
  Compute the final time for the simulation, which should correspond to one full
  rotation of the wave around the domain. This is based on the domain length and
  the advection speed.
  """
  x_left = float(cfg["mesh"]["x_left"])
  x_right = float(cfg["mesh"]["x_right"])
  return (x_right - x_left) / abs(float(cfg["pde"]["advection_speed"]))

def main() -> int:
  """
  Main function to run the advection spatial convergence study. It will:
  - Parse command-line arguments.
  - Validate the base YAML config.
  - For each polynomial order and refinement level:
    - Generate a YAML config for that case.
    - Run the solver and capture output.
    - Compute the L2 error against the exact solution.
  - Compute observed orders of convergence.
  - Print tables and save summary CSV and plots.
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
  if args.levels < 2: 
    raise ValueError("--levels must be at least 2 to measure an observed order.")

  validate_space_settings(P_ORDER, INTEGRATOR_BY_P, CFL_BY_INTEGRATOR)
  base_cfg                 = load_yaml(config_path)
  validate_base_config(base_cfg)
  final_time               = compute_final_time(base_cfg)
  exact                    = make_advection_exact_solution(INITIAL_CONDITION_EXPRESSION)
  configs_dir, outputs_dir = prepare_results_dir(results_dir)

  all_results = {}

  summary_rows = []
  for p_order in P_ORDER:

    integrator = INTEGRATOR_BY_P[p_order]
    cfl        = float(CFL_BY_INTEGRATOR[integrator])

    rows = []
    for level in range(args.levels):

      n_elements      = int(args.n_elements * (2 ** level))
      case_name       = f"p{p_order}_{integrator}_N{n_elements:04d}"
      case_output_dir = outputs_dir / case_name
      case_config     = configs_dir / f"{case_name}.yaml"

      write_case_config(
        base_cfg=base_cfg, 
        integrator=integrator, 
        n_elements=n_elements, 
        p_order=p_order, 
        cfl_number=cfl, 
        final_time=final_time, 
        initial_condition_expression=INITIAL_CONDITION_EXPRESSION, 
        output_dir=case_output_dir, 
        config_path=case_config
      )
      
      stdout, _       = run_solver(exe=exe, config_path=case_config, cwd=repo_root)
      n_steps, t_final_reported = parse_solver_stdout(stdout)
      final_file      = case_output_dir / "solution_final.dat"

      if not final_file.exists(): 
        raise FileNotFoundError(f"Expected final output file not found: {final_file}")

      error_l2, h, n_elements_in_file = compute_l2_error_against_exact(final_file, exact)

      if n_elements_in_file != n_elements: 
        raise RuntimeError(f"Mismatch between requested N={n_elements} and file Ne={n_elements_in_file} in {final_file}")

      rows.append(
        { "integrator": integrator, 
          "N": float(n_elements), 
          "p": float(p_order), 
          "h": h, 
          "cfl": cfl, 
          "dt_avg": average_dt(final_time, n_steps, t_final_reported), 
          "n_steps": float(n_steps) if n_steps is not None else math.nan, 
          "error_l2": error_l2, 
          "t_final": float(t_final_reported) if t_final_reported is not None else math.nan
        }
      )
    # end for level

    eoc_h = pairwise_eoc_from_x([r["h"] for r in rows], [r["error_l2"] for r in rows])

    for row, eoc_i in zip(rows, eoc_h):
      row["eoc_h"] = eoc_i
      summary_rows.append(row.copy())

    all_results[f"p{p_order}"] = rows

    print_space_table("DGmini advection spatial convergence", rows, INITIAL_CONDITION_EXPRESSION)

  # end for p_order

  write_summary_csv(
    results_dir / "summary.csv", 
    summary_rows, 
    ["integrator", "p", "N", "h", "cfl", "dt_avg", "n_steps", "error_l2", "eoc_h", "t_final"]
  )
  
  make_plot(
    results_dir / "error_vs_h.png", 
    all_results, 
    x_builder=lambda rows: [r["h"] for r in rows], 
    y_key="error_l2", 
    xlabel="h", 
    ylabel=r"$L^2$ error vs exact", 
    title="DGmini advection spatial convergence after one period", 
    label_builder=lambda _, rows, slope: f"p={int(rows[0]['p'])}, {rows[0]['integrator']}, CFL={rows[0]['cfl']:g} (slope ~ {slope:.2f})"
  )
  
  make_plot(
    results_dir / "error_vs_N.png", 
    all_results, 
    x_builder=lambda rows: [r["N"] for r in rows], 
    y_key="error_l2", 
    xlabel="number of elements", 
    ylabel=r"$L^2$ error vs exact", 
    title="DGmini advection spatial convergence after one period", 
    label_builder=lambda _, rows, slope: f"p={int(rows[0]['p'])}, {rows[0]['integrator']}, CFL={rows[0]['cfl']:g} (slope ~ {slope:.2f})"
  
  )

  print()
  print(f"Results written to: {results_dir}")
  print("  - summary.csv")
  print("  - error_vs_h.png")
  print("  - error_vs_N.png")
  print(f"  - generated YAML configs in {configs_dir}")
  print(f"  - solver outputs in {outputs_dir}")

  return 0

if __name__ == "__main__":
  raise SystemExit(main())
