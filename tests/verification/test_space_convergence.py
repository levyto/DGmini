#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import pytest

from utils_verification import (
  ROOT_DIR,
  SCRIPTS_DIR,
  check_space_orders,
  run_convergence_script,
)


SPACE_CASES = [
  (
    "advection",
    SCRIPTS_DIR / "run_advection_convergence_space.py",
    ROOT_DIR / "examples" / "advection.yaml",
  ),
  (
    "burgers",
    SCRIPTS_DIR / "run_burgers_convergence_space.py",
    ROOT_DIR / "examples" / "burgers.yaml",
  ),
]

@pytest.mark.parametrize(
  ("case_name", "script_path", "config_path"),
  SPACE_CASES,
  ids=[case[0] for case in SPACE_CASES],
)
def test_space_convergence_orders(
  tmp_path: Path,
  case_name: str,
  script_path: Path,
  config_path: Path,
) -> None:
  results_dir = tmp_path / f"{case_name}_space_verification"
  summary_csv = run_convergence_script(
    script_path=script_path,
    config_path=config_path,
    results_dir=results_dir,
  )
  check_space_orders(summary_csv)