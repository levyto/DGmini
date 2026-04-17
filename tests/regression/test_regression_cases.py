from pathlib import Path

import pytest

from regression_utils import run_regression_case


REGRESSION_DIR = Path(__file__).resolve().parent

CASE_DIRS = [
    d for d in REGRESSION_DIR.iterdir()
    if d.is_dir() and (d / "config.yaml").exists() and (d / "reference.dat").exists()
]


@pytest.mark.parametrize("case_dir", CASE_DIRS, ids=[d.name for d in CASE_DIRS])
def test_regression_case(case_dir: Path) -> None:
    run_regression_case(case_dir)