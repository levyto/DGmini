from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT_DIR = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = ROOT_DIR / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from modal_solution_io import read_modal_solution


SOLVER_EXE = ROOT_DIR / "build" / "dgmini"


def infer_polynomial_order_from_coeffs(coeffs: np.ndarray) -> int:
    if coeffs.ndim != 2 or coeffs.shape[1] < 1:
        raise ValueError("Coefficient array must have shape (n_elements, n_coeffs).")
    return int(coeffs.shape[1] - 1)


def compute_l2_error_between_modal_files(
    solution_file: Path,
    reference_file: Path,
) -> float:
    x_left, x_right, coeffs = read_modal_solution(str(solution_file))
    x_left_ref, x_right_ref, coeffs_ref = read_modal_solution(str(reference_file))

    if len(x_left) != len(x_left_ref):
        raise ValueError(
            f"Element-count mismatch:\n"
            f"  solution:  {len(x_left)}\n"
            f"  reference: {len(x_left_ref)}"
        )

    if not np.allclose(x_left, x_left_ref) or not np.allclose(x_right, x_right_ref):
        raise ValueError("Mesh mismatch between solution and reference.")

    p = infer_polynomial_order_from_coeffs(coeffs)
    p_ref = infer_polynomial_order_from_coeffs(coeffs_ref)

    quad_order = max(2 * max(p, p_ref) + 8, 16)
    xi_q, w_q = np.polynomial.legendre.leggauss(quad_order)

    err_sq = 0.0

    for e in range(len(x_left)):
        xl = x_left[e]
        xr = x_right[e]
        ce = coeffs[e]
        ce_ref = coeffs_ref[e]

        J = 0.5 * (xr - xl)
        uh_q = np.polynomial.legendre.legval(xi_q, ce)
        uref_q = np.polynomial.legendre.legval(xi_q, ce_ref)

        err_sq += np.sum(w_q * (uh_q - uref_q) ** 2) * J

    return math.sqrt(err_sq)


def get_solution_file_from_config(config_file: Path) -> Path:
    with config_file.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    output_dir = Path(cfg["output"]["directory"])
    prefix = cfg["output"].get("prefix", "solution")

    return ROOT_DIR / output_dir / f"{prefix}_final.dat"


def run_solver_from_root(config_file: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(SOLVER_EXE), "--config", str(config_file.resolve())],
        cwd=str(ROOT_DIR),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )


def run_regression_case(case_dir: Path, l2_tol: float = 1.0e-10) -> None:
    config_file = case_dir / "config.yaml"
    reference_file = case_dir / "reference.dat"

    assert SOLVER_EXE.exists(), f"Solver executable not found: {SOLVER_EXE}"
    assert config_file.exists(), f"Missing config file: {config_file}"
    assert reference_file.exists(), f"Missing reference file: {reference_file}"

    solution_file = get_solution_file_from_config(config_file)

    if solution_file.exists():
        solution_file.unlink()

    result = run_solver_from_root(config_file)

    assert result.returncode == 0, (
        f"Solver failed for case: {case_dir.name}\n\n"
        f"STDOUT:\n{result.stdout}\n\n"
        f"STDERR:\n{result.stderr}"
    )

    assert solution_file.exists(), f"Expected output file not found: {solution_file}"

    error_l2 = compute_l2_error_between_modal_files(solution_file, reference_file)

    assert error_l2 < l2_tol, (
        f"Regression mismatch detected for case: {case_dir.name}\n"
        f"L2 error = {error_l2:.16e}\n"
        f"Tolerance = {l2_tol:.16e}"
    )