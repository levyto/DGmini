from __future__ import annotations

import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = ROOT_DIR / "scripts"
SOLVER_EXE = ROOT_DIR / "build" / "dgmini"

EXPECTED_TIME_ORDER = {
  "forward_euler": 1.0,
  "rk2": 2.0,
  "rk3_ssp": 3.0,
  "rk4": 4.0,
}

ORDER_TOL = 0.2


def fitted_order(x_values: list[float], errors: list[float]) -> float:
  """
  Fit a line to log-log data and return the slope, which corresponds to the 
  observed order of convergence.
  """
  x = np.asarray(x_values, dtype=float)
  e = np.asarray(errors, dtype=float)
  slope, _ = np.polyfit(np.log(x), np.log(e), 1)
  return float(slope)


def run_convergence_script(
  *,
  script_path: Path,
  config_path: Path,
  results_dir: Path,
) -> Path:
  """
  Run the specified convergence script and return the path to the generated summary.csv.
  Raises an assertion error if the script fails or if summary.csv is not found.
  """
  assert script_path.exists(), f"Missing script: {script_path}"
  assert SOLVER_EXE.exists(), f"Missing solver executable: {SOLVER_EXE}"
  assert config_path.exists(), f"Missing config file: {config_path}"

  cmd = [
    sys.executable,
    str(script_path),
    "--exe",
    str(SOLVER_EXE),
    "--config",
    str(config_path),
    "--results-dir",
    str(results_dir),
  ]

  proc = subprocess.run(
    cmd,
    cwd=str(ROOT_DIR),
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
    check=False,
  )

  assert proc.returncode == 0, (
    f"Convergence script failed: {script_path.name}\\n\\n"
    f"STDOUT:\\n{proc.stdout}\\n\\n"
    f"STDERR:\\n{proc.stderr}"
  )

  summary_csv = results_dir / "summary.csv"
  assert summary_csv.exists(), f"Missing summary.csv: {summary_csv}"
  return summary_csv


def read_space_summary(summary_csv: Path) -> dict[int, list[dict[str, float]]]:
  """
  Read the space convergence summary CSV and return a dictionary grouped by polynomial order `p`.
  Each entry contains a list of dictionaries with keys: "integrator", "h", "error_l2", and "N".
  """
  grouped: dict[int, list[dict[str, float]]] = defaultdict(list)

  with summary_csv.open("r", encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
      p = int(float(row["p"]))
      grouped[p].append(
        {
          "integrator": row["integrator"],
          "h": float(row["h"]),
          "error_l2": float(row["error_l2"]),
          "N": float(row["N"]),
        }
      )

  for rows in grouped.values():
    rows.sort(key=lambda r: r["h"], reverse=True)

  return dict(grouped)


def read_time_summary(summary_csv: Path) -> dict[str, list[dict[str, float]]]:
  """
  Read the time convergence summary CSV and return a dictionary grouped by integrator.
  Each entry contains a list of dictionaries with keys: "dt_avg", "error_l2", "p", and "N".
  """
  grouped: dict[str, list[dict[str, float]]] = defaultdict(list)

  with summary_csv.open("r", encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
      grouped[row["integrator"]].append(
        {
          "dt_avg": float(row["dt_avg"]),
          "error_l2": float(row["error_l2"]),
          "p": float(row["p"]),
          "N": float(row["N"]),
        }
      )

  for rows in grouped.values():
    rows.sort(key=lambda r: r["dt_avg"], reverse=True)

  return dict(grouped)


def check_space_orders(summary_csv: Path, *, order_tol: float = ORDER_TOL) -> None:
  """
  Check the observed spatial convergence orders against the expected orders.
  Raises an assertion error if any measured order is outside the specified tolerance.
  """
  grouped = read_space_summary(summary_csv)
  failures: list[str] = []

  for p, rows in grouped.items():
    h = [r["h"] for r in rows]
    err = [r["error_l2"] for r in rows]
    integrator = rows[0]["integrator"]

    measured = fitted_order(h, err)
    expected = float(p + 1)
    o_min = expected - order_tol
    o_max = expected + order_tol

    if not (o_min <= measured <= o_max):
      failures.append(
        f"p={p} ({integrator}): measured order {measured:.4f}, "
        f"expected in [{o_min:.4f}, {o_max:.4f}]"
      )

  assert not failures, "Space verification failed:\\n" + "\\n".join(failures)


def check_time_orders(summary_csv: Path, *, order_tol: float = ORDER_TOL) -> None:
  """
  Check the observed time convergence orders against the expected orders.
  Raises an assertion error if any measured order is outside the specified tolerance.
  """
  grouped = read_time_summary(summary_csv)

  assert set(EXPECTED_TIME_ORDER).issubset(grouped.keys()), (
    "Not all expected integrators are present in summary.csv.\\n"
    f"Found: {sorted(grouped.keys())}"
  )

  failures: list[str] = []

  for integrator, q_expected in EXPECTED_TIME_ORDER.items():
    rows = grouped[integrator]
    dt = [r["dt_avg"] for r in rows]
    err = [r["error_l2"] for r in rows]

    q_measured = fitted_order(dt, err)
    q_min = q_expected - order_tol
    q_max = q_expected + order_tol

    if not (q_min <= q_measured <= q_max):
      failures.append(
        f"{integrator}: measured order {q_measured:.4f}, "
        f"expected in [{q_min:.4f}, {q_max:.4f}]"
      )

  assert not failures, "Time verification failed:\\n" + "\\n".join(failures)
  