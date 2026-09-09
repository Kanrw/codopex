# codopex

Generate symmetry-inequivalent co-doping defect-pair and higher-order
complexes in a host supercell.

**codopex** automates the configurational enumeration behind co-doping
studies: given a pristine supercell and a set of substitution reactions
(e.g. Cu at Pb sites, S at O sites), it

1. discovers every inequivalent *single-defect type* (a dopant on one
   pristine site-symmetry class, labelled `{dopant}_{host}_{symmetry}`,
   e.g. `Cu_Pb_C3v`),
2. enumerates all symmetry-inequivalent *defect pairs* and returns the
   nearest-neighbour (default) or next-nearest-neighbour complexes of every
   pair combo (with repetition, including `A+A`),
3. expands any pair by a third dopant into *triple complexes* placed near
   A, near B, or near both (`nearA_farB`, `nearB_farA`, `both_nn`), at the
   nearest or next-nearest shell.

Everything is computed in memory with pymatgen/spglib — no external
enumerator, no file side effects.  Structures keep the pristine site order;
export helpers write POSCAR trees and manifests when you need VASP inputs.

```text
codopex 为您生成共掺杂缺陷复合体的不等价构型：
第一阶段给定超胞与取代反应，生成全部缺陷对（最近邻，可选加次近邻）；
第二阶段指定缺陷对并添加第三种掺杂，按 near A / near B / near both
三种放置判据生成三缺陷复合体。纯 pymatgen 实现，无外部枚举器。
```

## Installation

```bash
pip install -e .          # from the repository root
```

Requires Python ≥ 3.10, `pymatgen>=2024` and `numpy`.

## Quick start

```python
import codopex as cp

# Phase 1: nearest-neighbour pairs among Cu@Pb, S@O and S@P
pairs = cp.generate_pairs(
    "CONTCAR",                          # pristine supercell (Structure or path)
    [("Pb", "Cu"), ("O", "S"), ("P", "S")],
)

pairs.type_labels
# ['Cu_Pb_C3v', 'Cu_Pb_D3d', 'S_O_Cs', 'S_O_C3v', 'S_P_C3v']

for row in cp.pair_rows(pairs):
    print(row["combo"], row["shell"], row["dist_AB"])

# include the next-nearest shell as well
pairs = cp.generate_pairs("CONTCAR", [("Pb", "Cu"), ("O", "S"), ("P", "S")],
                          shells="nnn")

# Phase 2: expand one pair by an additional S@O dopant
triples = cp.generate_triples(
    pairs,
    pair="Cu_Pb_C3v+S_O_C3v",
    third=[("O", "S")],     # optional; default: all reactions of the pairs
)
for row in cp.triple_rows(triples):
    print(row["combo"], row["criterion"], row["d_AC"], row["d_BC"])
```

### What you get back

`PairsResult` / `TriplesResult` hold the structures themselves together
with full metadata:

| attribute | content |
|---|---|
| `.types` | discovered `DefectType`s in canonical order (reactions in the order given, then site order within a reaction) |
| `.bases` | one base structure per type: the configuration whose dopant sits closest to the supercell center |
| `.complexes` | per combo: `.candidates` (every inequivalent configuration with distances, degeneracy/orbit size and structure) and `.shells` / `.selections` (the representatives) |
| `.stats` | per-job configuration counts, skipped combos |

Every candidate records its **degeneracy** (orbit size of the substituted
site), which is the weight needed for thermodynamic averages.

### Exporting for VASP runs

```python
cp.io.export_pairs(pairs, "out_pairs", write_all=False)
cp.io.export_triples(triples, "out_triples")
cp.io.write_rows_csv(cp.pair_rows(pairs), "pairs.csv")
```

The export mirrors the layout of the original notebooks:
`single_defect/{type}/POSCAR`, `{combo}/{nn|nnn}/POSCAR`,
`{combo}/{criterion}/POSCAR`, with an unindexed `POSCAR` copy per
representative and all candidates under `work_{combo}/` when
`write_all=True`.  POSCAR files are written species-grouped by atomic
number (VASP convention).

## How it works

* **Enumeration engine** — a single substitution on a base structure is
  symmetry-inequivalent up to the space group of that base; the
  inequivalent configurations are exactly the orbits of the candidate host
  sites under spglib.  Orbit representatives (smallest site index) give a
  deterministic order.  This reproduces the SAGAR enumeration used in the
  original workflow (see `docs/regression.md`).
* **Type discovery** — host-site classes are defined by the pristine site
  symmetry (Schoenflies) plus a Wyckoff tag; classes sharing a symmetry
  label stay distinct types.  Labels follow the original convention
  `{dopant}_{host}_{symmetry}`.
* **Shells** — distances are clustered with a 0.05 A tolerance; "nn" is
  the first cluster, "nnn" the second (a new shell starts past
  start + 0.05 A).  Within a shell, several inequivalent configurations may
  exist; codopex keeps one deterministic representative per (combo, shell)
  and exposes all candidates.
* **Pair combos** — all `combinations_with_replacement` over the types
  (so `A+A` pairs are included and require at least two sites of the type).
* **Triples** — a triple whose sorted types are `(i, j, c)` is generated
  from the representative of its *canonical parent pair* `(i, j)`.  The
  placement criteria of the original workflow select one representative per
  criterion; all configurations are kept in `.candidates` with d_AC/d_BC.

A/B convention: A is the pair's base-type dopant, B the added one; C is
the third dopant.  Site order is never changed: site `i` of any returned
structure is pristine site `i` with (possibly) a different element.

## Testing

```bash
pip install -e ".[dev]"
python -m pytest tests/
```

* `tests/test_synthetic.py` — mechanics on small crystals (shells, base
  rule, criteria semantics, exports, edge cases).
* `tests/test_ppo_regression.py` — regression against the original
  SAGAR-based workflow on the Pb3(PO4)2 Cu–S case (`tests/data/ppo`).
  See `docs/regression.md` for what is exact, what differs by design, and
  why.

## Publishing

The repository is set up to publish as a public GitHub project (MIT).
The `gh` CLI is not installed on the development machine, so publishing is
manual:

```bash
# create an empty public repository named codopex on github.com, then
git remote add origin git@github.com:<you>/codopex.git
git push -u origin main
```

## Limitations

* Only *substitutional* complexes are enumerated (vacancies and
  interstitials are labelled by `classify_defects` but are not part of the
  generation pipeline yet).
* The input is a single pristine supercell; codopex does not build or
  enlarge supercells and does not check that the cell is large enough for
  the defects of interest (it warns only when a site class is exhausted).
* Representative *selection inside* an equal-distance shell is
  deterministic but not unique; a different (equally valid) configuration
  may be picked than an earlier SAGAR-based run (see `docs/regression.md`).
