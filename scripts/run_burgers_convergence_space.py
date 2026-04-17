#!/usr/bin/env python3
"""Run a 1D Burgers' *space-convergence* study for DGmini.

This script measures spatial accuracy for the inviscid Burgers equation:
  1. Reads a base YAML config (typically examples/burgers.yaml).
  2. Loops over the polynomial orders listed in P_ORDER.
  3. For each p, selects the associated time integrator through INTEGRATOR_BY_P
     and the corresponding CFL through CFL_BY_INTEGRATOR.
  4. Runs a sequence of meshes N, 2N, 4N, ... for each (p, integrator) pair.
  5. Computes the L2 error against the exact *pre-shock* solution at the final time.
  6. Prints convergence tables and writes summary CSV + log-log plots.

Notes:
  - This script is intended only for times before shock formation.
  - For the default initial condition u0(x) = sin(2*pi*x), the shock forms at
        t_s = - 1/min(u'_0(x)) = 1 / (2*pi) ≈ 0.159155.
    Therefore the default final time T = 0.15 is still in the smooth regime.
  - The exact solution is evaluated through the method of characteristics:
        u(x,t) = u0(xi),    where    x = xi + t*u0(xi).
    The footpoint xi is found numerically for each evaluation point.
  - To isolate spatial accuracy, choose the CFL small enough that temporal error
    remains below spatial error on the tested meshes.
  - The initial condition expression is injected into the YAML from
    INITIAL_CONDITION_EXPRESSION below and is also used to define the exact solution.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import yaml

from modal_solution_io import read_modal_solution


# ------------------------------------------------------------------------------
# Default settings
# ------------------------------------------------------------------------------

# Polynomial orders to test.
P_ORDER = [1, 2, 3, 4]

# Time integrator assigned to each polynomial order.
INTEGRATOR_BY_P = {
    1: "rk3_ssp",
    2: "rk3_ssp",
    3: "rk3_ssp",
    4: "rk4",
}

# CFL associated with each time integrator.
CFL_BY_INTEGRATOR = {
    "forward_euler": 0.05,
    "rk2": 0.10,
    "rk3_ssp": 0.1,
    "rk4": 0.1,
}

# Base number of elements. The script refines this as N, 2N, 4N, ...
N_ELEM = 50

# Initial condition used both in YAML and in the exact solution before shock formation.
INITIAL_CONDITION_EXPRESSION = "sin(2*pi*x)"

# Final time for the convergence test. Should be before shock formation for the chosen IC.
FINAL_TIME = 0.1


# ------------------------------------------------------------------------------
# Argument parsing and config helpers
# ------------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--exe",
        type=Path,
        default=Path("./build/dgmini"),
        help="Path to DGmini executable.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("./examples/burgers.yaml"),
        help="Base YAML config.",
    )
    parser.add_argument(
        "--levels",
        type=int,
        default=4,
        help="Number of spatial refinement levels.",
    )
    parser.add_argument(
        "--n-elements",
        type=int,
        default=N_ELEM,
        help="Base number of elements. The study uses N, 2N, 4N, ...",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("./output/burgers_convergence_space"),
        help="Directory where generated configs, outputs, tables, and plots are stored.",
    )
    return parser.parse_args()


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def dump_yaml(data: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def validate_base_config(cfg: dict) -> None:
    pde_type = cfg["pde"]["type"]
    if pde_type != "burgers":
        raise ValueError(f"This script expects pde.type = 'burgers', got {pde_type!r}.")

    tsc_type = cfg["time_step_controller"]["type"]
    if tsc_type != "cfl":
        raise ValueError(
            "This script is intended for CFL-driven runs so that dt scales with h. "
            f"Got time_step_controller.type = {tsc_type!r}."
        )


def validate_settings() -> None:
    if not P_ORDER:
        raise ValueError("P_ORDER must not be empty.")

    if len(set(P_ORDER)) != len(P_ORDER):
        raise ValueError(f"P_ORDER contains duplicates: {P_ORDER}")

    for p in P_ORDER:
        if int(p) < 0:
            raise ValueError(f"Polynomial order must be nonnegative, got {p} in P_ORDER.")

        if p not in INTEGRATOR_BY_P:
            raise KeyError(f"Missing integrator mapping for p={p} in INTEGRATOR_BY_P.")

        integrator = INTEGRATOR_BY_P[p]
        if integrator not in CFL_BY_INTEGRATOR:
            raise KeyError(
                f"Missing CFL for integrator {integrator!r} referenced by p={p} in INTEGRATOR_BY_P."
            )

    unknown_integrators = sorted(set(CFL_BY_INTEGRATOR) - {"forward_euler", "rk2", "rk3_ssp", "rk4"})
    if unknown_integrators:
        raise KeyError(f"Unknown integrator names in CFL_BY_INTEGRATOR: {unknown_integrators}")

    for name, value in CFL_BY_INTEGRATOR.items():
        if value <= 0.0:
            raise ValueError(f"CFL_BY_INTEGRATOR[{name!r}] must be positive, got {value}.")


def get_case_settings(p_order: int) -> Tuple[str, float]:
    integrator = INTEGRATOR_BY_P[p_order]
    cfl = float(CFL_BY_INTEGRATOR[integrator])
    return integrator, cfl


# ------------------------------------------------------------------------------
# Exact solution and error evaluation
# ------------------------------------------------------------------------------

def infer_polynomial_order_from_coeffs(coeffs: np.ndarray) -> int:
    if coeffs.ndim != 2 or coeffs.shape[1] < 1:
        raise ValueError("Coefficient array must have shape (n_elements, n_coeffs).")
    return int(coeffs.shape[1] - 1)


def compute_final_time(cfg: dict) -> float:
    return FINAL_TIME


def make_initial_condition_function(expression: str) -> Callable[[np.ndarray], np.ndarray]:
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
        local_env = dict(allowed_names)
        local_env["x"] = x
        values = eval(expression, {"__builtins__": {}}, local_env)
        return np.asarray(values, dtype=float)

    return func


def make_burgers_exact_pre_shock_solution(
    expression: str,
    final_time: float,
    x_left: float,
    x_right: float,
) -> Callable[[np.ndarray], np.ndarray]:
    u0 = make_initial_condition_function(expression)
    length = x_right - x_left

    def exact(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        xi = x.copy()

        # Newton iteration for:
        # x = xi + t * u0(xi)
        # F(xi) = xi + t*u0(xi) - x = 0
        for _ in range(30):
            u = u0(xi)

            # Numerical derivative of u0(xi)
            eps = 1.0e-8
            du_dxi = (u0(xi + eps) - u0(xi - eps)) / (2.0 * eps)

            F = xi + final_time * u - x
            dF = 1.0 + final_time * du_dxi

            delta = F / dF
            xi -= delta

            # periodic wrap-back to domain
            xi = x_left + np.mod(xi - x_left, length)

            if np.max(np.abs(delta)) < 1.0e-13:
                break

        return u0(xi)

    return exact


def compute_l2_error_against_exact(
    modal_file: Path,
    exact: Callable[[np.ndarray], np.ndarray],
) -> Tuple[float, float, int]:
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


# ------------------------------------------------------------------------------
# Solver IO
# ------------------------------------------------------------------------------

def write_case_config(
    base_cfg: dict,
    integrator: str,
    n_elements: int,
    p_order: int,
    cfl_number: float,
    output_dir: Path,
    config_path: Path,
    final_time: float,
) -> None:
    cfg = yaml.safe_load(yaml.safe_dump(base_cfg))

    cfg["mesh"]["n_elements"] = int(n_elements)
    cfg["fem"]["order"] = int(p_order)
    cfg["time_integrator"]["type"] = integrator
    cfg["run"]["final_time"] = float(final_time)
    cfg["time_step_controller"]["cfl"] = float(cfl_number)

    cfg["initial_condition"]["expression"] = INITIAL_CONDITION_EXPRESSION

    cfg["output"]["directory"] = str(output_dir.resolve())
    cfg["output"]["prefix"] = "solution"
    cfg["output"]["output_dt"] = float(final_time + 1.0)
    cfg["output"]["write_initial"] = False
    cfg["output"]["write_final"] = True

    dump_yaml(cfg, config_path)


def parse_solver_stdout(stdout: str) -> Tuple[int | None, float | None]:
    step_match = re.search(r"after\s+(\d+)\s+steps", stdout)
    time_match = re.search(r"Finished at t\s*=\s*([0-9eE+\-.]+)", stdout)
    n_steps = int(step_match.group(1)) if step_match else None
    final_time = float(time_match.group(1)) if time_match else None
    return n_steps, final_time


def run_solver(exe: Path, config_path: Path, cwd: Path) -> Tuple[str, str]:
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


# ------------------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------------------

def observed_order_from_fit(x: Iterable[float], err: Iterable[float]) -> float:
    x_arr = np.asarray(list(x), dtype=float)
    e_arr = np.asarray(list(err), dtype=float)
    slope, _ = np.polyfit(np.log(x_arr), np.log(e_arr), 1)
    return float(slope)


def pairwise_eoc_from_x(x_values: List[float], errors: List[float]) -> List[float | None]:
    eoc: List[float | None] = [None]
    for i in range(1, len(errors)):
        numerator = math.log(errors[i - 1] / errors[i])
        denominator = math.log(x_values[i - 1] / x_values[i])
        eoc.append(numerator / denominator)
    return eoc


def print_table(
    integrator: str,
    p: int,
    cfl: float,
    results: List[Dict[str, float]],
) -> None:
    base_n = int(results[0]["N"]) if results else -1

    print()
    print(
        f"p = {p} | integrator = {integrator} | base N = {base_n} | CFL = {cfl:.6g} | "
        f"IC = {INITIAL_CONDITION_EXPRESSION}"
    )
    print("-" * 120)
    print(f"{'N':>8} {'h':>14} {'dt_avg':>14} {'steps':>10} {'L2 error':>18} {'EOC_h':>10}")
    print("-" * 120)
    for row in results:
        eoc_val = "-" if row["eoc_h"] is None else f"{row['eoc_h']:.4f}"
        steps_val = "-" if math.isnan(row["n_steps"]) else f"{int(row['n_steps'])}"
        print(
            f"{int(row['N']):8d} {row['h']:14.6e} {row['dt_avg']:14.6e} {steps_val:>10} "
            f"{row['error_l2']:18.10e} {eoc_val:>10}"
        )
    print("-" * 120)
    overall_order = observed_order_from_fit(
        [r["h"] for r in results],
        [r["error_l2"] for r in results],
    )
    print(f"overall fitted spatial order = {overall_order:.4f}")


def write_summary_csv(path: Path, all_results: Dict[str, List[Dict[str, float]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "integrator",
                "p",
                "N",
                "h",
                "cfl",
                "dt_avg",
                "n_steps",
                "error_l2",
                "eoc_h",
            ]
        )
        for _, rows in all_results.items():
            for row in rows:
                writer.writerow(
                    [
                        row["integrator"],
                        int(row["p"]),
                        int(row["N"]),
                        row["h"],
                        row["cfl"],
                        row["dt_avg"],
                        row["n_steps"],
                        row["error_l2"],
                        row["eoc_h"],
                    ]
                )


def make_plots(results_dir: Path, all_results: Dict[str, List[Dict[str, float]]]) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for _, rows in all_results.items():
        h = [r["h"] for r in rows]
        err = [r["error_l2"] for r in rows]
        order = observed_order_from_fit(h, err)
        p_i = int(rows[0]["p"])
        integrator_i = rows[0]["integrator"]
        cfl_i = rows[0]["cfl"]
        ax.loglog(
            h,
            err,
            marker="o",
            label=f"p={p_i}, {integrator_i}, CFL={cfl_i:g} (slope ~ {order:.2f})",
        )
    ax.set_xlabel("h")
    ax.set_ylabel(r"$L^2$ error vs exact")
    ax.set_title("DGmini Burgers' spatial convergence after one period")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(results_dir / "error_vs_h.png", dpi=200)
    plt.close(fig)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for _, rows in all_results.items():
        n_elem = [r["N"] for r in rows]
        err = [r["error_l2"] for r in rows]
        order = observed_order_from_fit([1.0 / n for n in n_elem], err)
        p_i = int(rows[0]["p"])
        integrator_i = rows[0]["integrator"]
        cfl_i = rows[0]["cfl"]
        ax.loglog(
            n_elem,
            err,
            marker="o",
            label=f"p={p_i}, {integrator_i}, CFL={cfl_i:g} (slope ~ {order:.2f})",
        )
    ax.set_xlabel("number of elements")
    ax.set_ylabel(r"$L^2$ error vs exact")
    ax.set_title("DGmini Burgers' spatial convergence after one period")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(results_dir / "error_vs_N.png", dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------------------
# Main driver
# ------------------------------------------------------------------------------

def main() -> int:
    args = parse_args()

    exe = args.exe.resolve()
    config_path = args.config.resolve()
    results_dir = args.results_dir.resolve()
    repo_root = config_path.parent.parent if config_path.parent.name == "examples" else config_path.parent

    if not exe.exists():
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")
    if args.n_elements <= 0:
        raise ValueError("--n-elements must be positive.")
    if args.levels < 2:
        raise ValueError("--levels must be at least 2 to measure an observed order.")

    validate_settings()

    base_cfg = load_yaml(config_path)
    validate_base_config(base_cfg)
    final_time = compute_final_time(base_cfg)
    
    x_left = float(base_cfg["mesh"]["x_left"])
    x_right = float(base_cfg["mesh"]["x_right"])

    exact = make_burgers_exact_pre_shock_solution(
        INITIAL_CONDITION_EXPRESSION,
        final_time,
        x_left,
        x_right,
    )

    configs_dir = results_dir / "configs"
    outputs_dir = results_dir / "solver_output"
    if results_dir.exists():
        shutil.rmtree(results_dir)
    configs_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    all_results: Dict[str, List[Dict[str, float]]] = {}

    for p_order in P_ORDER:
        integrator, cfl = get_case_settings(p_order)
        rows: List[Dict[str, float]] = []

        for level in range(args.levels):
            n_elements = int(args.n_elements * (2 ** level))
            case_name = f"p{p_order}_{integrator}_N{n_elements:04d}"
            case_output_dir = outputs_dir / case_name
            case_config = configs_dir / f"{case_name}.yaml"

            write_case_config(
                base_cfg=base_cfg,
                integrator=integrator,
                n_elements=n_elements,
                p_order=p_order,
                cfl_number=cfl,
                output_dir=case_output_dir,
                config_path=case_config,
                final_time=final_time,
            )

            stdout, _ = run_solver(exe=exe, config_path=case_config, cwd=repo_root)
            n_steps, t_final_reported = parse_solver_stdout(stdout)

            final_file = case_output_dir / "solution_final.dat"
            if not final_file.exists():
                raise FileNotFoundError(f"Expected final output file not found: {final_file}")

            error_l2, h, n_elements_in_file = compute_l2_error_against_exact(final_file, exact)
            if n_elements_in_file != n_elements:
                raise RuntimeError(
                    f"Mismatch between requested N={n_elements} and file Ne={n_elements_in_file} "
                    f"in {final_file}"
                )

            if n_steps is None or n_steps <= 0:
                dt_avg = math.nan
            else:
                final_time_for_dt = t_final_reported if t_final_reported is not None else final_time
                dt_avg = float(final_time_for_dt) / float(n_steps)

            rows.append(
                {
                    "integrator": integrator,
                    "N": float(n_elements),
                    "p": float(p_order),
                    "h": h,
                    "cfl": cfl,
                    "dt_avg": dt_avg,
                    "n_steps": float(n_steps) if n_steps is not None else math.nan,
                    "error_l2": error_l2,
                    "t_final": float(t_final_reported) if t_final_reported is not None else math.nan,
                }
            )

        h_values = [r["h"] for r in rows]
        errors = [r["error_l2"] for r in rows]
        eoc_h = pairwise_eoc_from_x(h_values, errors)
        for row, eoc_i in zip(rows, eoc_h):
            row["eoc_h"] = eoc_i

        all_results[f"p{p_order}"] = rows
        print_table(integrator, p_order, cfl, rows)

    write_summary_csv(results_dir / "summary.csv", all_results)
    make_plots(results_dir, all_results)

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
