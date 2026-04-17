from __future__ import annotations

import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = ROOT_DIR / "scripts"

SPACE_SCRIPT = SCRIPTS_DIR / "run_burgers_convergence_space.py"
SOLVER_EXE = ROOT_DIR / "build" / "dgmini"
CONFIG_FILE = ROOT_DIR / "examples" / "burgers.yaml"

ORDER_TOL = 0.2


def fitted_order(x_values: list[float], errors: list[float]) -> float:
    x = np.asarray(x_values, dtype=float)
    e = np.asarray(errors, dtype=float)
    slope, _ = np.polyfit(np.log(x), np.log(e), 1)
    return float(slope)


def read_summary(summary_csv: Path) -> dict[int, list[dict[str, float]]]:
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


def test_space_convergence_orders(tmp_path: Path) -> None:
    assert SPACE_SCRIPT.exists(), f"Missing script: {SPACE_SCRIPT}"
    assert SOLVER_EXE.exists(), f"Missing solver executable: {SOLVER_EXE}"
    assert CONFIG_FILE.exists(), f"Missing config file: {CONFIG_FILE}"

    results_dir = tmp_path / "space_verification"

    cmd = [
        sys.executable,
        str(SPACE_SCRIPT),
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
        "Space convergence script failed.\n\n"
        f"STDOUT:\n{proc.stdout}\n\n"
        f"STDERR:\n{proc.stderr}"
    )

    summary_csv = results_dir / "summary.csv"
    assert summary_csv.exists(), f"Missing summary.csv: {summary_csv}"

    grouped = read_summary(summary_csv)
    failures: list[str] = []

    for p, rows in grouped.items():
        h = [r["h"] for r in rows]
        err = [r["error_l2"] for r in rows]
        integrator = rows[0]["integrator"]

        measured = fitted_order(h, err)
        expected = float(p + 1)
        o_min = expected - ORDER_TOL
        o_max = expected + ORDER_TOL

        if not (o_min <= measured <= o_max):
            failures.append(
                f"p={p} ({integrator}): measured order {measured:.4f}, "
                f"expected in [{o_min:.4f}, {o_max:.4f}]"
            )

    assert not failures, "Space verification failed:\n" + "\n".join(failures)
