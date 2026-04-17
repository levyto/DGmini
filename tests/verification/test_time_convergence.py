from __future__ import annotations

import csv
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = ROOT_DIR / "scripts"

TIME_SCRIPT = SCRIPTS_DIR / "run_advection_convergence_time.py"
SOLVER_EXE = ROOT_DIR / "build" / "dgmini"
CONFIG_FILE = ROOT_DIR / "examples" / "advection.yaml"

EXPECTED_TIME_ORDER = {
    "forward_euler": 1.0,
    "rk2": 2.0,
    "rk3_ssp": 3.0,
    "rk4": 4.0,
}

ORDER_TOL = 0.2


def fitted_order(x_values: list[float], errors: list[float]) -> float:
    x = np.asarray(x_values, dtype=float)
    e = np.asarray(errors, dtype=float)
    slope, _ = np.polyfit(np.log(x), np.log(e), 1)
    return float(slope)


def read_summary(summary_csv: Path) -> dict[str, list[dict[str, float]]]:
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


def test_time_convergence_orders(tmp_path: Path) -> None:
    assert TIME_SCRIPT.exists(), f"Missing script: {TIME_SCRIPT}"
    assert SOLVER_EXE.exists(), f"Missing solver executable: {SOLVER_EXE}"
    assert CONFIG_FILE.exists(), f"Missing config file: {CONFIG_FILE}"

    results_dir = tmp_path / "time_verification"

    cmd = [
        sys.executable,
        str(TIME_SCRIPT),
        "--exe",
        str(SOLVER_EXE),
        "--config",
        str(CONFIG_FILE),
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
        "Time convergence script failed.\n\n"
        f"STDOUT:\n{proc.stdout}\n\n"
        f"STDERR:\n{proc.stderr}"
    )

    summary_csv = results_dir / "summary.csv"
    assert summary_csv.exists(), f"Missing summary.csv: {summary_csv}"

    grouped = read_summary(summary_csv)

    assert set(EXPECTED_TIME_ORDER).issubset(grouped.keys()), (
        "Not all expected integrators are present in summary.csv.\n"
        f"Found: {sorted(grouped.keys())}"
    )

    failures: list[str] = []

    for integrator, q_expected in EXPECTED_TIME_ORDER.items():
        rows = grouped[integrator]
        dt = [r["dt_avg"] for r in rows]
        err = [r["error_l2"] for r in rows]

        q_measured = fitted_order(dt, err)
        q_min = q_expected - ORDER_TOL
        q_max = q_expected + ORDER_TOL

        if not (q_min <= q_measured <= q_max):
            failures.append(
                f"{integrator}: measured order {q_measured:.4f}, expected in [{q_min:.4f}, {q_max:.4f}]"
            )

    assert not failures, "Time verification failed:\n" + "\n".join(failures)
