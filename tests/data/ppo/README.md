# PPO regression fixtures

These files reproduce the Pb3(PO4)2 co-doping case of the original workflow
that codopex was extracted from (R-3m `Pb3(PO4)2`, 78-atom supercell, Cu-S
codoping).  They serve as the reference for the regression test
`tests/test_ppo_regression.py`.

Provenance:

* `bulk.vasp` — pristine supercell used as the input of the original
  notebooks (DFT-relaxed CONTCAR of the defect-free bulk).
* `single_defect/{type}/POSCAR` — the five base structures (one per defect
  type) as produced by the original SAGAR-based `DefectMaker` pipeline.
* `pairs/{combo}/POSCAR` — the fifteen nearest-neighbour pair structures
  (original pipeline output).
* `pairs_expected.csv` — shell distances produced by codopex, the SAGAR
  nearest distance, and whether the codopex representative structure is
  symmetry-equivalent to the SAGAR file (StructureMatcher).
* `triples_expected.csv` — triple representatives produced by codopex with
  the SAGAR reference distances and a per-row verification mode:
  - `verified`: codopex value equals the SAGAR value exactly (up to
    1e-3 A), possibly with the A/B orientation swapped for pairs whose two
    dopants are symmetry-equivalent;
  - `shell`: the values lie in the same distance shell (0.05 A) of the
    SAGAR reference; representative choices inside a shell are not
    reproducible when several inequivalent configurations share the shell
    distance (the original pipeline picked one arbitrarily; see
    `docs/regression.md`).

The structures are the authors' own DFT inputs/outputs and are distributed
under the repository license for testing purposes.
