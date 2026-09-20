# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/Kanrw/codopex/commits/main
[0.1.0]: https://github.com/Kanrw/codopex/releases
