#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Callable, Iterable

import matplotlib.pyplot as plt
import numpy as np
import yaml

from utils_modal_solution_io import (
  read_modal_solution
)


# ------------------------------------------------------------------------------
# YAML / filesystem / solver IO
# ------------------------------------------------------------------------------

def load_yaml(path: Path) -> dict:
  """
  Load a YAML file and return its contents as a dictionary.
  """
  with path.open("r", encoding="utf-8") as f:
    return yaml.safe_load(f)


def dump_yaml(data: dict, path: Path) -> None:
  """
  Dump a dictionary to a YAML file.
  """
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", encoding="utf-8") as f:
    yaml.safe_dump(data, f, sort_keys=False)


def infer_repo_root_from_config(config_path: Path) -> Path:
  """
  Infer the DGmini repository root directory from the location of a config file.
  If the config file is located in a subdirectory named "examples", returns its 
  parent directory. Otherwise, returns the parent directory of the config file.
  """
  return config_path.parent.parent if config_path.parent.name == "examples" else config_path.parent


def prepare_results_dir(results_dir: Path) -> tuple[Path, Path]:
  """
  Prepare the results directory by creating subdirectories for configs and 
  solver outputs.
  
  If the results directory already exists, it will be removed and recreated.
  
  Returns the paths to the configs and outputs directories.
  """
  configs_dir = results_dir / "configs"
  outputs_dir = results_dir / "solver_output"

  if results_dir.exists():
    shutil.rmtree(results_dir)

  configs_dir.mkdir(parents=True, exist_ok=True)
  outputs_dir.mkdir(parents=True, exist_ok=True)

  return configs_dir, outputs_dir


def write_case_config(
  base_cfg: dict,
  *,
  integrator: str,
  n_elements: int,
  p_order: int,
  cfl_number: float,
  final_time: float,
  initial_condition_expression: str,
  output_dir: Path,
  config_path: Path,
) -> None:
  """
  Write a YAML configuration file for a single convergence test case, based on a
  base configuration and specific parameters for the case.
  """
  cfg = yaml.safe_load(yaml.safe_dump(base_cfg))

  cfg["mesh"]["n_elements"] = int(n_elements)
  cfg["fem"]["order"] = int(p_order)
  cfg["time_integrator"]["type"] = integrator
  cfg["run"]["final_time"] = float(final_time)
  cfg["time_step_controller"]["cfl"] = float(cfl_number)
  cfg["initial_condition"]["expression"] = initial_condition_expression

  cfg["output"]["directory"] = str(output_dir.resolve())
  cfg["output"]["prefix"] = "solution"
  cfg["output"]["output_dt"] = float(final_time + 1.0)
  cfg["output"]["write_initial"] = False
  cfg["output"]["write_final"] = True

  dump_yaml(cfg, config_path)


def run_solver(exe: Path, config_path: Path, cwd: Path) -> tuple[str, str]:
  """
  Run the DGmini solver with the specified executable and configuration file,
  in the given working directory. Returns a tuple of (stdout, stderr) from the
  solver process. Raises a RuntimeError if the solver returns a non-zero exit 
  code.
  """  
  cmd = [str(exe.resolve()), "--config", str(config_path.resolve())]
  proc = subprocess.run(
    cmd,
    cwd=str(cwd),
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
    check=False,
  )

  if proc.returncode != 0:
    raise RuntimeError(
      f"DGmini run failed for config {config_path}\n"
      f"returncode = {proc.returncode}\n"
      f"stdout:\n{proc.stdout}\n"
      f"stderr:\n{proc.stderr}"
    )

  return proc.stdout, proc.stderr


def parse_solver_stdout(stdout: str) -> tuple[int | None, float | None]:
  """
  Parse the solver's standard output to extract the number of time steps taken
  and the final time reached. Returns a tuple of (n_steps, final_time). If either 
  piece of information cannot be found, its value in the tuple will be None.
  """
  step_match = re.search(r"after\s+(\d+)\s+steps", stdout)
  time_match = re.search(r"Finished at t\s*=\s*([0-9eE+\-.]+)", stdout)

  n_steps = int(step_match.group(1)) if step_match else None
  final_time = float(time_match.group(1)) if time_match else None
  return n_steps, final_time


def average_dt(final_time: float, n_steps: int | None, reported_final_time: float | None) -> float:
  """
  Compute the average time step size given the final time, number of steps, and
  optionally the reported final time from the solver. If the number of steps is
  None or non-positive, returns NaN.
  """
  if n_steps is None or n_steps <= 0:
    return math.nan

  time_for_dt = reported_final_time if reported_final_time is not None else final_time
  return float(time_for_dt) / float(n_steps)


def validate_space_settings(
  p_order: list[int],
  integrator_by_p: dict[int, str],
  cfl_by_integrator: dict[str, float],
) -> None:
  """
  Validate the settings related to spatial convergence tests, including 
  polynomial orders, their associated time integrators, and the CFL numbers 
  for each integrator.  Raises ValueError or KeyError if any inconsistencies 
  or invalid values are found.
  """ 
  if not p_order:
    raise ValueError("P_ORDER must not be empty.")

  if len(set(p_order)) != len(p_order):
    raise ValueError(f"P_ORDER contains duplicates: {p_order}")

  valid_integrators = {"forward_euler", "rk2", "rk3_ssp", "rk4"}

  for p in p_order:
    if int(p) < 0:
      raise ValueError(f"Polynomial order must be nonnegative, got {p} in P_ORDER.")
    if p not in integrator_by_p:
      raise KeyError(f"Missing integrator mapping for p={p} in INTEGRATOR_BY_P.")
    integrator = integrator_by_p[p]
    if integrator not in cfl_by_integrator:
      raise KeyError(f"Missing CFL for integrator {integrator!r} referenced by p={p} in INTEGRATOR_BY_P.")

  unknown_integrators = sorted(set(cfl_by_integrator) - valid_integrators)
  if unknown_integrators:
    raise KeyError(f"Unknown integrator names in CFL_BY_INTEGRATOR: {unknown_integrators}")

  for name, value in cfl_by_integrator.items():
    if value <= 0.0:
      raise ValueError(f"CFL_BY_INTEGRATOR[{name!r}] must be positive, got {value}.")


def validate_time_settings(
  integrators: list[str],
  cfl_by_integrator: dict[str, float],
  p_by_integrator: dict[str, int],
  ref_cfl_by_integrator: dict[str, float],
) -> None:
  """
  Validate the settings related to temporal convergence tests, including time 
  integrators, their associated CFL numbers, polynomial orders, and reference 
  CFL numbers. Raises ValueError or KeyError if any inconsistencies or invalid 
  values are found.
  """
  valid_integrators = {"forward_euler", "rk2", "rk3_ssp", "rk4"}

  unknown_cfl = sorted(set(cfl_by_integrator) - valid_integrators)
  unknown_p = sorted(set(p_by_integrator) - valid_integrators)
  unknown_ref_cfl = sorted(set(ref_cfl_by_integrator) - valid_integrators)

  if unknown_cfl:
    raise KeyError(f"Unknown integrator names in CFL_BY_INTEGRATOR: {unknown_cfl}")
  if unknown_p:
    raise KeyError(f"Unknown integrator names in P_BY_INTEGRATOR: {unknown_p}")
  if unknown_ref_cfl:
    raise KeyError(f"Unknown integrator names in REF_CFL_BY_INTEGRATOR: {unknown_ref_cfl}")

  for integrator in integrators:
    if integrator not in valid_integrators:
      raise KeyError(f"Unknown integrator {integrator!r} in INTEGRATORS.")
    if integrator not in cfl_by_integrator:
      raise KeyError(f"Missing CFL for integrator {integrator!r}.")
    if integrator not in p_by_integrator:
      raise KeyError(f"Missing polynomial order for integrator {integrator!r}.")
    if integrator not in ref_cfl_by_integrator:
      raise KeyError(f"Missing reference CFL for integrator {integrator!r}.")

  for name, value in cfl_by_integrator.items():
    if value <= 0.0:
      raise ValueError(f"CFL_BY_INTEGRATOR[{name!r}] must be positive, got {value}.")
  for name, value in p_by_integrator.items():
    if int(value) < 0:
      raise ValueError(f"P_BY_INTEGRATOR[{name!r}] must be nonnegative, got {value}.")
  for name, value in ref_cfl_by_integrator.items():
    if value <= 0.0:
      raise ValueError(f"REF_CFL_BY_INTEGRATOR[{name!r}] must be positive, got {value}.")



# ------------------------------------------------------------------------------
# Exact solution / IC helpers
# ------------------------------------------------------------------------------

def make_initial_condition_function(expression: str) -> Callable[[np.ndarray], np.ndarray]:
  """
  Create a function that evaluates the given initial condition expression.
  The expression can use numpy functions and the variable 'x'.
  """
  allowed_names = {
    "np": np,
    "pi": np.pi,
    "x": None,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "exp": np.exp,
    "sqrt": np.sqrt,
    "log": np.log,
    "abs": np.abs,
    "tanh": np.tanh,
    "sinh": np.sinh,
    "cosh": np.cosh,
  }

  def func(x: np.ndarray) -> np.ndarray:
    """
    Evaluate the initial condition expression at the given x values.
    """
    local_env = dict(allowed_names)
    local_env["x"] = x
    values = eval(expression, {"__builtins__": {}}, local_env)
    return np.asarray(values, dtype=float)

  return func


def make_advection_exact_solution(expression: str) -> Callable[[np.ndarray], np.ndarray]:
  """
  Create a function that evaluates the exact solution for the advection equation.
  The exact solution is simply the initial condition, as one full period is 
  assumed.
  """
  return make_initial_condition_function(expression)


def make_burgers_exact_pre_shock_solution(
  expression: str,
  final_time: float,
  x_left: float,
  x_right: float,
) -> Callable[[np.ndarray], np.ndarray]:
  """
  Create a function that evaluates the exact solution for the inviscid Burgers' 
  equation before shock formation. This is done by using the method of 
  characteristics and performing a fixed-point iteration to find the 
  characteristic foot for each spatial point at the final time. The initial 
  condition is given by the expression, and the solution is assumed to be 
  periodic on the interval [x_left, x_right].
  """
  u0 = make_initial_condition_function(expression)
  length = x_right - x_left

  def exact(x: np.ndarray) -> np.ndarray:
    """
    Evaluate the exact solution at the given x values by performing a fixed-point
    iteration to find the characteristic foot for each point.
    """
    x = np.asarray(x, dtype=float)
    xi = x.copy()

    for _ in range(30):
      u = u0(xi)
      eps = 1.0e-8
      du_dxi = (u0(xi + eps) - u0(xi - eps)) / (2.0 * eps)

      F = xi + final_time * u - x
      dF = 1.0 + final_time * du_dxi

      delta = F / dF
      xi -= delta
      xi = x_left + np.mod(xi - x_left, length)

      if np.max(np.abs(delta)) < 1.0e-13:
        break

    return u0(xi)

  return exact


# ------------------------------------------------------------------------------
# Error evaluators
# ------------------------------------------------------------------------------

def infer_polynomial_order_from_coeffs(coeffs: np.ndarray) -> int:
  """
  Infer the polynomial order from the shape of the coefficient array. The array 
  is expected to have shape (n_elements, n_coeffs), where n_coeffs = p + 1 for a
  polynomial of order p.
  """
  if coeffs.ndim != 2 or coeffs.shape[1] < 1:
    raise ValueError("Coefficient array must have shape (n_elements, n_coeffs).")
  return int(coeffs.shape[1] - 1)


def compute_l2_error_against_exact(
  modal_file: Path,
  exact: Callable[[np.ndarray], np.ndarray],
) -> tuple[float, float, int]:
  """
  Compute the L2 error of the numerical solution against the exact solution,
  given a modal solution file and an exact solution function. Returns a tuple of
  (L2 error, average element size h, number of elements).
  """
  x_left, x_right, coeffs = read_modal_solution(str(modal_file))
  p = infer_polynomial_order_from_coeffs(coeffs)

  quad_order = max(2 * p + 8, 16)
  xi_q, w_q = np.polynomial.legendre.leggauss(quad_order)

  err_sq = 0.0
  h_values = x_right - x_left

  for e in range(len(x_left)):
    xl = x_left[e]
    xr = x_right[e]
    ce = coeffs[e]

    J = 0.5 * (xr - xl)
    x_q = 0.5 * (xl + xr) + J * xi_q
    uh_q = np.polynomial.legendre.legval(xi_q, ce)
    ue_q = exact(x_q)

    err_sq += np.sum(w_q * (uh_q - ue_q) ** 2) * J

  return math.sqrt(err_sq), float(np.mean(h_values)), int(len(x_left))


def compute_l2_error_against_reference(
  modal_file: Path,
  reference_file: Path,
) -> tuple[float, float, int]:
  """
  Compute the L2 error of the numerical solution against a reference solution,
  given modal solution files for both. Returns a tuple of (L2 error, average
  element size h, number of elements).
  """
  x_left, x_right, coeffs = read_modal_solution(str(modal_file))
  x_left_ref, x_right_ref, coeffs_ref = read_modal_solution(str(reference_file))

  if len(x_left) != len(x_left_ref):
    raise ValueError(
      f"Element-count mismatch between solution and reference:\n"
      f"  {modal_file}: {len(x_left)} elements\n"
      f"  {reference_file}: {len(x_left_ref)} elements"
    )

  if not np.allclose(x_left, x_left_ref) or not np.allclose(x_right, x_right_ref):
    raise ValueError(
      f"Mesh mismatch between solution and reference files:\n"
      f"  {modal_file}\n"
      f"  {reference_file}"
    )

  p = infer_polynomial_order_from_coeffs(coeffs)
  p_ref = infer_polynomial_order_from_coeffs(coeffs_ref)

  quad_order = max(2 * max(p, p_ref) + 8, 16)
  xi_q, w_q = np.polynomial.legendre.leggauss(quad_order)

  err_sq = 0.0
  h_values = x_right - x_left

  for e in range(len(x_left)):
    xl = x_left[e]
    xr = x_right[e]
    ce = coeffs[e]
    ce_ref = coeffs_ref[e]

    J = 0.5 * (xr - xl)
    uh_q = np.polynomial.legendre.legval(xi_q, ce)
    uref_q = np.polynomial.legendre.legval(xi_q, ce_ref)

    err_sq += np.sum(w_q * (uh_q - uref_q) ** 2) * J

  return math.sqrt(err_sq), float(np.mean(h_values)), int(len(x_left))


# ------------------------------------------------------------------------------
# Convergence math / reporting
# ------------------------------------------------------------------------------

def observed_order_from_fit(x: Iterable[float], err: Iterable[float]) -> float:
  """
  Compute the observed order of convergence from a fit of the error data.
  """
  x_arr = np.asarray(list(x), dtype=float)
  e_arr = np.asarray(list(err), dtype=float)
  slope, _ = np.polyfit(np.log(x_arr), np.log(e_arr), 1)
  return float(slope)


def pairwise_eoc_from_x(x_values: list[float], errors: list[float]) -> list[float | None]:
  """
  Compute the pairwise experimental order of convergence (EOC) from the error data.
  """
  eoc: list[float | None] = [None]
  for i in range(1, len(errors)):
    numerator = math.log(errors[i - 1] / errors[i])
    denominator = math.log(x_values[i - 1] / x_values[i])
    eoc.append(numerator / denominator)
  return eoc


def write_summary_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
  """
  Write a summary CSV file with the given rows and fieldnames.
  """
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in rows:
      writer.writerow(row)


def print_space_table(case_label: str, results: list[dict], ic_expression: str) -> None:
  """
  Print a table summarizing the spatial convergence results for a given case.
  """
  base_n = int(results[0]["N"]) if results else -1
  p = int(results[0]["p"])
  integrator = results[0]["integrator"]
  cfl = results[0]["cfl"]

  print()
  print(
    f"{case_label} | p = {p} | integrator = {integrator} | base N = {base_n} | CFL = {cfl:.6g} | IC = {ic_expression}"
  )
  print("-" * 120)
  print(f"{'N':>8} {'h':>14} {'dt_avg':>14} {'steps':>10} {'L2 error':>18} {'EOC_h':>10}")
  print("-" * 120)

  for row in results:
    eoc_val = "-" if row["eoc_h"] is None else f"{row['eoc_h']:.4f}"
    steps_val = "-" if math.isnan(row["n_steps"]) else f"{int(row['n_steps'])}"
    print(
      f"{int(row['N']):8d} {row['h']:14.6e} {row['dt_avg']:14.6e} {steps_val:>10} {row['error_l2']:18.10e} {eoc_val:>10}"
    )

  print("-" * 120)
  overall_order = observed_order_from_fit([r["h"] for r in results], [r["error_l2"] for r in results])
  print(f"overall fitted spatial order = {overall_order:.4f}")


def print_time_table(case_label: str, results: list[dict]) -> None:
  """
  Print a table summarizing the temporal convergence results for a given case.
  """
  p = int(results[0]["p"])
  integrator = results[0]["integrator"]
  n_elements = int(results[0]["N"])
  ref_cfl = results[0]["ref_cfl"]
  ref_dt_avg = results[0]["ref_dt_avg"]
  start_cfl = results[0]["cfl"]

  print()
  print(
    f"{case_label} | integrator = {integrator} | p = {p} | N = {n_elements} | start CFL = {start_cfl:.6g} | ref CFL = {ref_cfl:.6g} | ref dt_avg = {ref_dt_avg:.6e}"
  )
  print("-" * 120)
  print(f"{'CFL':>10} {'dt_avg':>14} {'steps':>10} {'h':>14} {'L2 error':>18} {'EOC_t':>10}")
  print("-" * 120)

  for row in results:
    eoc_val = "-" if row["eoc_t"] is None else f"{row['eoc_t']:.4f}"
    steps_val = "-" if math.isnan(row["n_steps"]) else f"{int(row['n_steps'])}"
    print(
      f"{row['cfl']:10.4e} {row['dt_avg']:14.6e} {steps_val:>10} {row['h']:14.6e} {row['error_l2']:18.10e} {eoc_val:>10}"
    )

  print("-" * 120)
  overall_order = observed_order_from_fit([r["dt_avg"] for r in results], [r["error_l2"] for r in results])
  print(f"overall fitted temporal order = {overall_order:.4f}")


def make_plot(
  output_path: Path,
  grouped_results: dict[str, list[dict]],
  *,
  x_builder: Callable[[list[dict]], list[float]],
  y_key: str,
  xlabel: str,
  ylabel: str,
  title: str,
  label_builder: Callable[[str, list[dict], float], str],
) -> None:
  """
  Create and save a log-log plot of the convergence results, with one line per 
  group.
  """
  fig = plt.figure()
  ax = fig.add_subplot(111)

  for group_name, rows in grouped_results.items():
    x = x_builder(rows)
    y = [r[y_key] for r in rows]
    slope = observed_order_from_fit(x, y)
    ax.loglog(x, y, marker="o", label=label_builder(group_name, rows, slope))

  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ax.set_title(title)
  ax.grid(True, which="both")
  ax.legend()
  fig.tight_layout()
  fig.savefig(output_path, dpi=200)
  plt.close(fig)
