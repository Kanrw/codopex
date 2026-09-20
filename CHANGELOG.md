# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2026-09-20

Complexity-audit cleanup: dead code, unused API and redundant layers were
removed.  No numerical behavior changed and the SAGAR regression fixtures are
untouched; result files saved by 0.2.0 still load (unknown keys are ignored).

### Removed

- `codopex.io.to_dataframe()` — use `pandas.DataFrame(cp.pair_rows(...))` or
  `cp.triple_rows(...)`; pandas is now only a test dependency.
- `PairsResult.complex_for()` — use `pairs.complexes[combo]`.
- `AllTriplesResult.pair_labels` — use `list(batch.results)`.
- The `center` parameter of `generate_pairs()`; the supercell center is fixed
  at (0.5, 0.5, 0.5) as in the original workflow.
- The `type_c` field of `TripleCandidate` — read it from the parent
  `TripleComplex.type_c`; the `config_index` field of `SingleDefectBase` and
  the `shell` field of `PairShell`.
- The `codopex.cli` `__main__` guard; use the `codopex` console script.

## [0.2.0] - 2026-09-20

### Added

- Command-line interface: `codopex pairs`, `codopex triples` and
  `codopex classify` (`codopex --help`).
- `generate_all_triples()` batch API over every parent pair (or a subset),
  with `AllTriplesResult.skipped` for pairs lacking the requested shell.
- Result persistence: `codopex.io.save()` / `codopex.io.load()` for pair,
  triple and batch results (JSON, gzip when the path ends in `.gz`), so the
  two phases can run in separate sessions or machines.
- GitHub Actions CI (Python 3.10-3.14, ruff, mypy), `py.typed`, coverage
  configuration, `CITATION.cff`, `CONTRIBUTING.md` and this changelog.
- `codopex.__version__` is now read from a single source of truth.

### Changed

- `classify_defects()` pairs sites with a globally optimal one-to-one
  assignment (maximum number of matches below the distance threshold,
  minimum total displacement) computed from a single distance matrix;
  ambiguous matchings may differ from the previous per-site greedy pick.
- `classify_defects()` no longer labels interstitial atoms: a defect site
  without a pristine counterpart raises `ValueError`, because codopex only
  generates substitutional complexes.

## [0.1.0] - 2026-09-09

### Added

- Initial release: pure pymatgen/spglib engine that replaces the SAGAR-based
  enumerator of the original notebooks.
- Defect-type discovery from pristine site symmetry (Schoenflies + Wyckoff),
  single-defect base selection (dopant closest to the supercell center).
- Symmetry-inequivalent defect pairs, nearest (`nn`) and next-nearest
  (`nnn`) shells, all `A+B` combos with repetition.
- Triple complexes expanded from a parent pair with the `nearA_farB`,
  `nearB_farA` and `both_nn` placement criteria.
- `classify_defects()` for labelling vacancies and substitutions in
  externally produced structures.
- POSCAR tree and CSV manifest export helpers.
- Regression suite against the original SAGAR workflow on the
  Pb3(PO4)2 Cu-S case (`tests/data/ppo`, see `docs/regression.md`).

[Unreleased]: https://github.com/Kanrw/codopex/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/Kanrw/codopex/releases/tag/v0.3.0
[0.2.0]: https://github.com/Kanrw/codopex/releases/tag/v0.2.0
[0.1.0]: https://github.com/Kanrw/codopex/releases/tag/v0.1.0
