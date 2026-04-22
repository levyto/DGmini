#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import pytest

from utils_verification import (
  ROOT_DIR,
  SCRIPTS_DIR,
  check_time_orders,
  run_convergence_script,
)


TIME_CASES = [
  (
    "advection",
    SCRIPTS_DIR / "run_advection_convergence_time.py",
    ROOT_DIR / "examples" / "advection.yaml",
  ),
  (
    "burgers",
    SCRIPTS_DIR / "run_burgers_convergence_time.py",
    ROOT_DIR / "examples" / "burgers.yaml",
  ),
]


@pytest.mark.parametrize(
  ("case_name", "script_path", "config_path"),
  TIME_CASES,
  ids=[case[0] for case in TIME_CASES],
)
def test_time_convergence_orders(
  tmp_path: Path,
  case_name: str,
  script_path: Path,
  config_path: Path,
) -> None:
  results_dir = tmp_path / f"{case_name}_time_verification"
  summary_csv = run_convergence_script(
    script_path=script_path,
    config_path=config_path,
    results_dir=results_dir,
  )
  check_time_orders(summary_csv)