"""Shared inputs and result locations for the paper workflows."""
from pathlib import Path

from pas.config import DATA_PATHS


LMARENA_CSV = Path(DATA_PATHS["lmarena"]["cleaned"]).with_name("clean_summary_v2.csv")
LMARENA_RESULTS = Path("tmp/lmarena_ci_results")
