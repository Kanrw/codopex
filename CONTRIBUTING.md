# Contributing to codopex

Thanks for considering a contribution.  This project is research software:
correctness of the enumeration and the regression against the original
workflow matter more than feature count.

## Development setup

```bash
git clone https://github.com/Kanrw/codopex.git
cd codopex
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Checks

Run all of these before opening a pull request:

```bash
pytest                      # tests (includes the SAGAR regression fixtures)
ruff check .                # lint
ruff format --check .       # formatting
mypy                        # type checking
```

CI runs the same commands on Python 3.10-3.13.

## What to keep in mind

- **Structure invariant.**  Every structure returned by the pipeline keeps
  the pristine site order; substitutions only change species.  Do not add
  code that reorders or moves sites in returned structures.
- **Determinism.**  Orbit representatives, shell grouping and criterion
  selection must stay deterministic (ties broken by site/config index).
- **Regression fixtures.**  `tests/data/ppo/*` is reference output of the
  original SAGAR-based workflow.  Do not edit the fixtures or the expected
  CSVs to make a failing test pass; explain and document the deviation in
  `docs/regression.md` instead (see the "intentional deviations" section).
- **Public API.**  Anything exported from `codopex/__init__.py` is part of
  the public API; new parameters should have sensible defaults so existing
  calls keep working.
- **Docs.**  Update the README and the changelog (`CHANGELOG.md`,
  "Unreleased" section) when you change user-facing behavior.

## Reporting issues

Please include the codopex version, the Python/pymatgen/spglib versions, a
minimal input structure (or its composition and space group) and the full
traceback.  For symmetry-related issues, state the `symprec` values used.
