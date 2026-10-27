# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `TriplesResult.shells` records the tier ("nn"/"nnn") the placement
  criteria were evaluated in, and is saved/loaded with the result.  Saved
  files without the key load as `"nn"`.
- `codopex.triples.generate_triples(..., dist_tol=...)` (also
  `generate_all_triples`) overrides the shell tolerance used by the placement
  criteria; `None` (default) keeps `pairs_result.dist_tol`.
- `TriplesResult.stats["missing"]` counts the criteria that had no candidate
  in the requested shell; `codopex triples` also prints a warning for them.

### Changed

- Triple export layout is now `{combo}/{shells}/{criterion}/POSCAR...`
  (e.g. `.../nn/nearA_farB/POSCAR`).  The criteria keep their names in both
  tiers, so an `nn` and an `nnn` run into the same `--out` used to overwrite
  each other's representative `POSCAR`; they now write side by side.
- `--write-all` output is counted correctly: the CLI reports the number of
  POSCAR files actually on disk (bases, representatives and `work_*`).
- `codopex triples --dist-tol` now takes effect with `--load-pairs`
  (previously it was silently ignored).  `--type-order` and `--site-symprec`
  only apply when phase 1 is rebuilt and now raise instead of being ignored.
- `export_pairs` / `export_triples` also return the paths written under
  `work_{combo}/` when `write_all=True`.
- `io.save()` and `io.write_rows_csv()` create missing parent directories, so
  the README quickstart (`cp.io.save(pairs, "runs/pairs.json.gz")`) works.

### Removed

- Dead "empty combo" bookkeeping in the pair and triple pipelines
  (`stats["empty_combos"]`, `stats["empty"]` were always empty lists).

### Fixed

- A typo in `shells` / `pair_shell` now fails before any enumeration instead
  of after it, with the same message for all call sites.

## [0.3.1] - 2026-10-16

### Added

- `skills/codopex-workflow/` — an agent skill that walks a user (or an AI
  assistant) through the physically guided co-doping workflow of the paper:
  system setup with doped, pair prioritization by chemical scoring, DFT
  binding energies, expansion of bound complexes, and the equilibrium
  Fermi-level check.  Ships the procedure (`SKILL.md`), step references, the
  stdlib-only helper script `scripts/codopex_tools.py` (binding-energy curves,
  ΔE_b, many-body departure, tiers, abundance, Fermi-level shift) and fill-in
  templates.  The directory is excluded from ruff: the helper script is a
  standalone asset, not library code.

### Changed

- README now matches the second-round revision of the paper: the workflow
  figure is replaced by the pair-prioritization / evaluation-and-expansion
  version (tiers, complex–complex expansion, abundance and E_F^eq check), and
  the caption and introduction describe the workflow in the same terms.

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

[Unreleased]: https://github.com/Kanrw/codopex/compare/v0.3.1...HEAD
[0.3.1]: https://github.com/Kanrw/codopex/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/Kanrw/codopex/releases/tag/v0.3.0
[0.2.0]: https://github.com/Kanrw/codopex/releases/tag/v0.2.0
[0.1.0]: https://github.com/Kanrw/codopex/releases/tag/v0.1.0
