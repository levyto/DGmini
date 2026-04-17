#!/usr/bin/env python3
"""Run a 1D Burgers' *time-convergence* study for DGmini.

This script measures temporal accuracy for the inviscid Burgers equation on a fixed spatial discretization:
  1. Reads a base YAML config (typically examples/burgers.yaml).
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

Notes:
  - This script measures the temporal order against a fine-time numerical reference
    on the same spatial discretization, rather than against the analytical solution.
  - This is the preferred setup for Burgers' equation, since even in the smooth regime
    the exact solution is implicit, and after shock formation a smooth exact solution
    no longer exists.
  - For the default initial condition u0(x) = sin(2*pi*x), the shock forms at
        t_s = - 1/min(u'_0(x)) = 1 / (2*pi) ≈ 0.159155.
    Therefore the default final time T = 0.15 is still in the pre-shock regime.
  - The reference solution uses the same mesh and polynomial order as the tested run,
    but a much smaller CFL, so the measured error is intended to isolate the temporal
    discretization error.
  - The printed dt is the effective average time step dt_avg = T / n_steps.
    If the solver clips the last step to hit final_time exactly, the last actual step
    may be smaller.
  - The initial condition expression is injected into the YAML from
    INITIAL_CONDITION_EXPRESSION below.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import yaml

from modal_solution_io import read_modal_solution

# ------------------------------------------------------------------------------
# Default settings for each integrator, which can be overridden by CLI args.
# ------------------------------------------------------------------------------

# List of integrators to test. This can be overridden by CLI args.
INTEGRATORS = ["forward_euler", "rk2", "rk3_ssp", "rk4"]

# Starting CFL for each integrator. Each next level uses CFL/2.
CFL_BY_INTEGRATOR = {
    "forward_euler":    0.10,
    "rk2":              0.60,
    "rk3_ssp":          0.80,
    "rk4":              0.80,
}

# Polynomial order p for each integrator. This is used for both the reference 
# and test runs of that integrator.
P_BY_INTEGRATOR = {
    "forward_euler":    1,
    "rk2":              2,
    "rk3_ssp":          3,
    "rk4":              4,
}

# Reference CFL used for the fine-time reference solution of each integrator.
# This should be smaller than the finest tested CFL for that integrator.
REF_CFL_BY_INTEGRATOR = {
    "forward_euler": CFL_BY_INTEGRATOR.get("forward_euler", 0.10) / (2 ** 6),
    "rk2":           CFL_BY_INTEGRATOR.get("rk2",           0.60) / (2 ** 6),
    "rk3_ssp":       CFL_BY_INTEGRATOR.get("rk3_ssp",       0.80) / (2 ** 6),
    "rk4":           CFL_BY_INTEGRATOR.get("rk4",           0.80) / (2 ** 6),
}

# Base number of elements. The script refines this as N, 2N, 4N, ...
N_ELEM = 100

# Initial condition used both in YAML and in the exact solution before shock formation.
INITIAL_CONDITION_EXPRESSION = "sin(2*pi*x)"

# Final time for the convergence test. Should be before shock formation for the chosen IC.
FINAL_TIME = 0.15


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
        help="Number of temporal refinement levels.",
    )
    parser.add_argument(
        "--integrators",
        nargs="+",
        default=INTEGRATORS,
        choices=INTEGRATORS,
        help="Time integrators to test.",
    )
    parser.add_argument(
        "--n-elements",
        type=int,
        default=N_ELEM,
        help="Fixed number of elements for all runs.",
    )
    parser.add_argument(
        "--p-order",
        type=int,
        default=5,
        help="Fallback polynomial order p if an integrator is not listed in P_BY_INTEGRATOR.",
    )
    parser.add_argument(
        "--base-cfl",
        type=float,
        default=0.4,
        help="Fallback starting CFL if an integrator is not listed in CFL_BY_INTEGRATOR. Each next level uses CFL/2.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("./output/burgers_convergence_time"),
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


def validate_override_mappings() -> None:
    unknown_cfl = sorted(set(CFL_BY_INTEGRATOR) - set(INTEGRATORS))
    unknown_p = sorted(set(P_BY_INTEGRATOR) - set(INTEGRATORS))
    unknown_ref_cfl = sorted(set(REF_CFL_BY_INTEGRATOR) - set(INTEGRATORS))

    if unknown_cfl:
        raise KeyError(f"Unknown integrator names in CFL_BY_INTEGRATOR: {unknown_cfl}")
    if unknown_p:
        raise KeyError(f"Unknown integrator names in P_BY_INTEGRATOR: {unknown_p}")
    if unknown_ref_cfl:
        raise KeyError(f"Unknown integrator names in REF_CFL_BY_INTEGRATOR: {unknown_ref_cfl}")

    for name, value in CFL_BY_INTEGRATOR.items():
        if value <= 0.0:
            raise ValueError(f"CFL_BY_INTEGRATOR[{name!r}] must be positive, got {value}.")

    for name, value in P_BY_INTEGRATOR.items():
        if int(value) < 0:
            raise ValueError(f"P_BY_INTEGRATOR[{name!r}] must be nonnegative, got {value}.")

    for name, value in REF_CFL_BY_INTEGRATOR.items():
        if value <= 0.0:
            raise ValueError(f"REF_CFL_BY_INTEGRATOR[{name!r}] must be positive, got {value}.")


def get_integrator_settings(integrator: str, args: argparse.Namespace) -> Tuple[float, int, float]:
    base_cfl = float(CFL_BY_INTEGRATOR.get(integrator, args.base_cfl))
    p_order = int(P_BY_INTEGRATOR.get(integrator, args.p_order))
    ref_cfl = float(REF_CFL_BY_INTEGRATOR.get(integrator, base_cfl / (2 ** 6)))
    return base_cfl, p_order, ref_cfl


def infer_polynomial_order_from_coeffs(coeffs: np.ndarray) -> int:
    if coeffs.ndim != 2 or coeffs.shape[1] < 1:
        raise ValueError("Coefficient array must have shape (n_elements, n_coeffs).")
    return int(coeffs.shape[1] - 1)


def compute_final_time(cfg: dict) -> float:
    return FINAL_TIME


def compute_l2_error_against_reference(
    modal_file: Path,
    reference_file: Path,
) -> Tuple[float, float, int]:
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
    n_elements: int,
    ref_cfl: float,
    ref_dt_avg: float,
    results: List[Dict[str, float]],
) -> None:
    start_cfl = results[0]["cfl"] if results else float("nan")

    print()
    print(
        f"Integrator: {integrator} | p = {p} | N = {n_elements} | "
        f"start CFL = {start_cfl:.6g} | ref CFL = {ref_cfl:.6g} | ref dt_avg = {ref_dt_avg:.6e}"
    )
    print("-" * 120)
    print(f"{'CFL':>10} {'dt_avg':>14} {'steps':>10} {'h':>14} {'L2 error':>18} {'EOC_t':>10}")
    print("-" * 120)
    for row in results:
        eoc_val = "-" if row["eoc_t"] is None else f"{row['eoc_t']:.4f}"
        steps_val = "-" if math.isnan(row["n_steps"]) else f"{int(row['n_steps'])}"
        print(
            f"{row['cfl']:10.4e} {row['dt_avg']:14.6e} {steps_val:>10} "
            f"{row['h']:14.6e} {row['error_l2']:18.10e} {eoc_val:>10}"
        )
    print("-" * 120)
    overall_order = observed_order_from_fit(
        [r["dt_avg"] for r in results],
        [r["error_l2"] for r in results],
    )
    print(f"overall fitted temporal order = {overall_order:.4f}")


def write_summary_csv(path: Path, all_results: Dict[str, List[Dict[str, float]]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "integrator",
                "p",
                "N",
                "cfl",
                "dt_avg",
                "n_steps",
                "h",
                "ref_cfl",
                "ref_dt_avg",
                "ref_n_steps",
                "error_l2",
                "eoc_t",
            ]
        )
        for integrator, rows in all_results.items():
            for row in rows:
                writer.writerow(
                    [
                        integrator,
                        int(row["p"]),
                        int(row["N"]),
                        row["cfl"],
                        row["dt_avg"],
                        row["n_steps"],
                        row["h"],
                        row["ref_cfl"],
                        row["ref_dt_avg"],
                        row["ref_n_steps"],
                        row["error_l2"],
                        row["eoc_t"],
                    ]
                )


def make_plots(results_dir: Path, all_results: Dict[str, List[Dict[str, float]]]) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for integrator, rows in all_results.items():
        dt = [r["dt_avg"] for r in rows]
        err = [r["error_l2"] for r in rows]
        order = observed_order_from_fit(dt, err)
        p_i = rows[0]["p"]
        n_i = rows[0]["N"]
        ref_cfl = rows[0]["ref_cfl"]
        ax.loglog(
            dt,
            err,
            marker="o",
            label=f"{integrator}, p={p_i}, N={n_i}, ref CFL={ref_cfl:g} (slope ~ {order:.2f})",
        )
    ax.set_xlabel("dt_avg")
    ax.set_ylabel(r"$L^2$ error vs reference")
    ax.set_title("DGmini Burgers' temporal convergence after one period")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(results_dir / "error_vs_dt.png", dpi=200)
    plt.close(fig)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for integrator, rows in all_results.items():
        n_steps = [r["n_steps"] for r in rows]
        err = [r["error_l2"] for r in rows]
        order = observed_order_from_fit([1.0 / s for s in n_steps], err)
        p_i = rows[0]["p"]
        n_i = rows[0]["N"]
        ref_cfl = rows[0]["ref_cfl"]
        ax.loglog(
            n_steps,
            err,
            marker="o",
            label=f"{integrator}, p={p_i}, N={n_i}, ref CFL={ref_cfl:g} (slope ~ {order:.2f})",
        )
    ax.set_xlabel("number of time steps")
    ax.set_ylabel(r"$L^2$ error vs reference")
    ax.set_title("DGmini Burgers' temporal convergence after one period")
    ax.grid(True, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(results_dir / "error_vs_nsteps.png", dpi=200)
    plt.close(fig)


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
    if args.p_order < 0:
        raise ValueError("--p-order must be nonnegative.")
    if args.levels < 2:
        raise ValueError("--levels must be at least 2 to measure an observed order.")
    if args.base_cfl <= 0.0:
        raise ValueError("--base-cfl must be positive.")

    validate_override_mappings()

    base_cfg = load_yaml(config_path)
    validate_base_config(base_cfg)
    final_time = compute_final_time(base_cfg)

    configs_dir = results_dir / "configs"
    outputs_dir = results_dir / "solver_output"
    if results_dir.exists():
        shutil.rmtree(results_dir)
    configs_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    all_results: Dict[str, List[Dict[str, float]]] = {}

    for integrator in args.integrators:
        rows: List[Dict[str, float]] = []
        base_cfl_i, p_order_i, ref_cfl_i = get_integrator_settings(integrator, args)

        finest_test_cfl = float(base_cfl_i / (2 ** (args.levels - 1)))
        if ref_cfl_i >= finest_test_cfl:
            raise ValueError(
                f"For integrator {integrator!r}, REF_CFL_BY_INTEGRATOR = {ref_cfl_i} "
                f"must be smaller than the finest tested CFL = {finest_test_cfl}."
            )

        ref_case_name = f"{integrator}_p{p_order_i}_N{args.n_elements:04d}_reference"
        ref_output_dir = outputs_dir / ref_case_name
        ref_config = configs_dir / f"{ref_case_name}.yaml"

        write_case_config(
            base_cfg=base_cfg,
            integrator=integrator,
            n_elements=args.n_elements,
            p_order=p_order_i,
            cfl_number=ref_cfl_i,
            output_dir=ref_output_dir,
            config_path=ref_config,
            final_time=final_time,
        )

        ref_stdout, _ = run_solver(exe=exe, config_path=ref_config, cwd=repo_root)
        ref_n_steps, ref_t_final_reported = parse_solver_stdout(ref_stdout)

        ref_file = ref_output_dir / "solution_final.dat"
        if not ref_file.exists():
            raise FileNotFoundError(f"Expected reference output file not found: {ref_file}")

        if ref_n_steps is None or ref_n_steps <= 0:
            ref_dt_avg = math.nan
        else:
            ref_time_for_dt = ref_t_final_reported if ref_t_final_reported is not None else final_time
            ref_dt_avg = float(ref_time_for_dt) / float(ref_n_steps)

        for level in range(args.levels):
            cfl_number = float(base_cfl_i / (2 ** level))
            case_name = f"{integrator}_p{p_order_i}_N{args.n_elements:04d}_cfl{level}"
            case_output_dir = outputs_dir / case_name
            case_config = configs_dir / f"{case_name}.yaml"

            write_case_config(
                base_cfg=base_cfg,
                integrator=integrator,
                n_elements=args.n_elements,
                p_order=p_order_i,
                cfl_number=cfl_number,
                output_dir=case_output_dir,
                config_path=case_config,
                final_time=final_time,
            )

            stdout, _ = run_solver(exe=exe, config_path=case_config, cwd=repo_root)
            n_steps, t_final_reported = parse_solver_stdout(stdout)

            final_file = case_output_dir / "solution_final.dat"
            if not final_file.exists():
                raise FileNotFoundError(f"Expected final output file not found: {final_file}")

            error_l2, h, n_elements_in_file = compute_l2_error_against_reference(final_file, ref_file)
            if n_elements_in_file != args.n_elements:
                raise RuntimeError(
                    f"Mismatch between requested N={args.n_elements} and file Ne={n_elements_in_file} "
                    f"in {final_file}"
                )

            if n_steps is None or n_steps <= 0:
                dt_avg = math.nan
            else:
                final_time_for_dt = t_final_reported if t_final_reported is not None else final_time
                dt_avg = float(final_time_for_dt) / float(n_steps)

            rows.append(
                {
                    "N": float(args.n_elements),
                    "p": float(p_order_i),
                    "h": h,
                    "cfl": cfl_number,
                    "dt_avg": dt_avg,
                    "n_steps": float(n_steps) if n_steps is not None else math.nan,
                    "ref_cfl": ref_cfl_i,
                    "ref_dt_avg": ref_dt_avg,
                    "ref_n_steps": float(ref_n_steps) if ref_n_steps is not None else math.nan,
                    "error_l2": error_l2,
                    "t_final": float(t_final_reported) if t_final_reported is not None else math.nan,
                }
            )

        dt_values = [r["dt_avg"] for r in rows]
        errors = [r["error_l2"] for r in rows]
        eoc_t = pairwise_eoc_from_x(dt_values, errors)
        for row, eoc_i in zip(rows, eoc_t):
            row["eoc_t"] = eoc_i

        all_results[integrator] = rows
        print_table(integrator, p_order_i, args.n_elements, ref_cfl_i, ref_dt_avg, rows)

    write_summary_csv(results_dir / "summary.csv", all_results)
    make_plots(results_dir, all_results)

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
