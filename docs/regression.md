# Regression against the original SAGAR workflow

codopex replaces the enumerator of the original co-doping notebooks
(pyvaspflow's `DefectMaker`, built on SAGAR) with a pure
pymatgen/spglib engine.  This document states what the engine reproduces
exactly, what it does *not* reproduce by design, and how the regression
(`tests/test_ppo_regression.py`) verifies each level.

## What is verified exactly

On the Pb3(PO4)2 case (fixtures in `tests/data/ppo`, references produced by
the original SAGAR-based workflow):

1. **Defect-type discovery** — same five types, same labels, same canonical
   order: `Cu_Pb_C3v, Cu_Pb_D3d, S_O_Cs, S_O_C3v, S_P_C3v`.
2. **Base structures** — the five single-defect base structures are
   symmetry-equivalent to the SAGAR files.  codopex selects, for each type,
   the site closest to the supercell center among *all* sites of the
   type's pristine site class.
3. **Pair distances** — all fifteen nearest distances equal the SAGAR
   values to machine precision on the fixtures (the original files carry
   ~1e-4 A writing noise; the committed tolerance is 1e-3 A).
4. **Pair representatives** — fourteen of fifteen nearest structures are
   symmetry-equivalent to the SAGAR files (StructureMatcher).
5. **Triple criteria** — 98 of 105 (combo, criterion) rows match the SAGAR
   values exactly or with the documented A/B orientation swap; the
   remaining 7 rows lie in the same distance shell.

## Intentional deviations (representative choice)

The notebook picks one representative *inside* each distance shell.  When
several symmetry-inequivalent configurations share the same shell distance
(exact distance ties) the choice depends on SAGAR's internal enumeration
order, which is arbitrary and not reproducible.  codopex picks
deterministically (lowest configuration index); a different but equally
valid representative may result.  Consequences:

* `Cu_Pb_D3d+Cu_Pb_D3d` (nn) — three inequivalent configurations at the
  same distance: the codopex representative is not the SAGAR file, while
  the distance is identical.  Triples grown from this pair are
  correspondingly affected (marked `shell` in the expectation table).
* A few triple criteria on the `Cu_Pb_C3v+Cu_Pb_C3v` pair select different
  in-shell maxima because the SAGAR-written coordinates carry ~1e-4 A
  noise that shifts near-degenerate candidates across the 0.05 A window;
  all reported distances are within one shell of the reference.

A/B convention for the criteria: codopex always labels the pair's
base-type dopant as A and the added dopant as B.  For pairs whose two
dopants are symmetry-equivalent (same element, same site symmetry) the
original notebook fell back to file order, which may swap A and B; the
regression therefore accepts the swapped orientation for those pairs
(`nearA_farB` compared against `nearB_farA` with d_AC/d_BC exchanged).

## Shell definitions

* Pair shell: distances between the two dopants, clustered with a
  tolerance of 0.05 A ("nn" = first cluster, "nnn" = second).  A new shell
  starts when a distance exceeds the previous shell start by more than the
  tolerance (the original notebook's algorithm).
* Triple criteria: `nearA_farB` = third dopant in the chosen shell of
  d_AC, maximizing d_BC; `nearB_farA` symmetric; `both_nn` = in the chosen
  shell of both distances, minimizing the mean (falling back to the overall
  mean-distance minimum).

## Reproducing the full dev-time parity check

The exhaustive per-job comparison against the SAGAR artifacts (configuration
counts and degeneracies per enumeration job, structure matching of all
representatives) requires the original output directories, which live
outside this repository:

```
/Users/mlwang/Documents/Research/cslp/gen_classify_codoping/
    codope_260812/            # 15 pair folders + work_* jobs + deg.txt
    codope_260812_3def/       # 35 triple folders + CSVs
```

The per-job enumeration (configuration counts and degeneracy multisets)
matched SAGAR exactly for all eleven pair jobs and all triple jobs.
