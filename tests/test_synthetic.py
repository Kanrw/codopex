"""Unit tests on small synthetic crystals.

These verify the mechanics of the pipeline (enumeration, shells, criteria,
labels, exports) on structures whose symmetry is easy to reason about; the
physics-level equivalence with the original SAGAR workflow is covered by
``test_ppo_regression.py``.
"""

import os

import pytest
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, Structure

import codopex as cp
from codopex.pipeline import shell_groups

matcher = StructureMatcher(primitive_cell=False, attempt_supercell=False)


def cubic_2x2x2():
    """Simple cubic lattice (a=4 A), 8 host atoms, one symmetry class."""
    return Structure(Lattice.cubic(4.0), ["Na"] * 8,
                     [[x / 2, y / 2, z / 2]
                      for x in (0, 1) for y in (0, 1) for z in (0, 1)])


def rocksalt_2x2x2():
    """Rocksalt supercell: 64 atoms, Na and Cl each one symmetry class."""
    base = Structure(Lattice.cubic(4.0),
                     ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    return base * [2, 2, 2]


def rocksalt_1x1x1():
    return Structure(Lattice.cubic(4.0),
                     ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])


# ---------------------------------------------------------------------------
# shells
# ---------------------------------------------------------------------------


def test_shell_groups():
    # group boundary: strictly beyond start + tol opens a new shell
    assert shell_groups([1.0, 1.049, 1.051, 1.06], 0.05) == [0, 0, 1, 1]
    assert shell_groups([1.0, 1.051, 1.10, 1.20], 0.05) == [0, 1, 1, 2]
    assert shell_groups([1.0, 1.04, 1.049], 0.05) == [0, 0, 0]
    assert shell_groups([], 0.05) == []


# ---------------------------------------------------------------------------
# single-defect types and bases
# ---------------------------------------------------------------------------


def test_type_discovery_and_base_rule():
    res = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")])
    # one site class -> one type
    assert res.type_labels == ["Li_Na_Oh"]
    # base = dopant closest to the cell center (0.5,0.5,0.5)
    base = res.bases["Li_Na_Oh"]
    assert base.dopant_site == 7  # the (0.5,0.5,0.5) site
    assert base.distance_to_center == pytest.approx(0.0, abs=1e-9)
    assert base.degeneracy == 8
    assert base.structure[7].species_string == "Li"


def test_reaction_dedup():
    a = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li"), ("Na", "Li")])
    b = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")])
    assert a.reactions == b.reactions == [("Na", "Li")]


def test_type_label_reports_pristine_symmetry():
    res = cp.generate_pairs(rocksalt_2x2x2(), [("Na", "Li"), ("Cl", "F")])
    # rock-salt sites are octahedral (m-3m -> Oh)
    assert res.type_labels == ["Li_Na_Oh", "F_Cl_Oh"]


# ---------------------------------------------------------------------------
# pairs on the cubic cell: shells are exactly lattice-neighbor shells
# ---------------------------------------------------------------------------


def test_pair_shells_are_neighbor_shells():
    res = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    comp = res.complexes["Li_Na_Oh+Li_Na_Oh"]
    # first shell = true nearest-neighbour distance (half the lattice constant)
    assert comp.shells[0].label == "nn"
    assert comp.shells[0].distance == pytest.approx(2.0, abs=1e-6)
    # shells appear in increasing distance order with the documented labels
    labels = [s.label for s in comp.shells]
    assert labels[0] == "nn"
    assert labels == ["nn"] or labels == ["nn", "nnn"]
    assert all(
        comp.shells[i].distance <= comp.shells[i + 1].distance
        for i in range(len(comp.shells) - 1)
    )
    # shell assignments are consistent with the shell-grouping algorithm
    by_d = sorted(comp.candidates, key=lambda c: c.distance)
    assert [c.shell for c in by_d] == shell_groups(
        [c.distance for c in by_d], 0.05
    )
    assert len(comp.candidates) >= 2


def test_nn_only_by_default():
    res = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")])  # default shells="nn"
    comp = res.complexes["Li_Na_Oh+Li_Na_Oh"]
    assert [s.label for s in comp.shells] == ["nn"]


def test_positions_preserved_and_single_substitution():
    bulk = cubic_2x2x2()
    res = cp.generate_pairs(bulk, [("Na", "Li")])
    comp = res.complexes["Li_Na_Oh+Li_Na_Oh"]
    for cand in comp.candidates:
        assert len([s for s in cand.structure if s.species_string == "Li"]) == 2
        for i in range(len(bulk)):
            assert (cand.structure[i].frac_coords - bulk[i].frac_coords).max() < 1e-9


# ---------------------------------------------------------------------------
# degenerate combos: not enough sites
# ---------------------------------------------------------------------------


def test_self_pair_needs_two_sites():
    res = cp.generate_pairs(rocksalt_1x1x1(), [("Na", "Li")])
    # only one Na site: the Li@Na + Li@Na combo cannot exist
    assert "Li_Na_Oh+Li_Na_Oh" not in res.complexes
    assert res.stats["single_defect"]["Li_Na_Oh"]["n_sites"] == 1


# ---------------------------------------------------------------------------
# triples: criteria semantics on the cubic cell
# ---------------------------------------------------------------------------


def test_triple_criteria_semantics():
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    tres = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh", shells="nn")
    tc = tres.complexes["Li_Na_Oh+Li_Na_Oh+Li_Na_Oh"]
    assert len(tc.candidates) >= 1
    cands = tc.candidates
    d_ac = sorted(c.d_ac for c in cands)
    d_bc = sorted(c.d_bc for c in cands)

    near_a = tc.selection("nearA_farB")
    near_b = tc.selection("nearB_farA")
    both = tc.selection("both_nn")

    # nearA: third dopant within the first d_AC shell, maximizing d_BC
    assert near_a.d_ac <= d_ac[0] + 0.05
    assert near_a.d_bc == max(c.d_bc for c in cands if c.d_ac <= d_ac[0] + 0.05)
    # nearB: mirror statement
    assert near_b.d_bc <= d_bc[0] + 0.05
    assert near_b.d_ac == max(c.d_ac for c in cands if c.d_bc <= d_bc[0] + 0.05)
    # both_nn: minimizes the mean distance over candidates NN to both
    both_batch = [c for c in cands if c.d_ac <= d_ac[0] + 0.05 and c.d_bc <= d_bc[0] + 0.05]
    pick = both_batch if both_batch else cands
    assert both.d_ac + both.d_bc == min(c.d_ac + c.d_bc for c in pick)


def test_triple_nnn_shells():
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    tres = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh", shells="nnn")
    tc = tres.complexes["Li_Na_Oh+Li_Na_Oh+Li_Na_Oh"]
    near_a = tc.selection("nearA_farB")
    d_ac = sorted(c.d_ac for c in tc.candidates)
    lo = d_ac[0] + 0.05  # boundary of first shell
    # with the pair at the cell center + a nearest neighbour, remaining sites
    # sit in distinct shells; the nnn criterion requires the second shell
    assert near_a is not None
    assert near_a.d_ac >= lo or near_a.d_ac > d_ac[0] + 1e-6


def test_pair_must_exist():
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")])
    with pytest.raises(ValueError, match="not among generated combos"):
        cp.generate_triples(pairs, "K_Na_Oh+K_Na_Oh")
    pairs2 = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nn")
    with pytest.raises(ValueError, match="has no 'nnn' representative"):
        cp.generate_triples(pairs2, "Li_Na_Oh+Li_Na_Oh", pair_shell="nnn")


# ---------------------------------------------------------------------------
# deterministic output
# ---------------------------------------------------------------------------


def test_deterministic():
    r1 = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    r2 = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    assert [c.distance for c in r1.complexes["Li_Na_Oh+Li_Na_Oh"].candidates] == \
           [c.distance for c in r2.complexes["Li_Na_Oh+Li_Na_Oh"].candidates]


# ---------------------------------------------------------------------------
# classify_defects utility
# ---------------------------------------------------------------------------


def test_classify_substitution_and_vacancy():
    bulk = rocksalt_2x2x2()
    sub = bulk.copy()
    sub.replace(3, "Li")
    labels = [d.defect_name for d in cp.classify_defects(bulk, sub)]
    assert labels == ["Li_Na_Oh"]

    vac = Structure(
        bulk.lattice,
        [s.species_string for i, s in enumerate(bulk) if i != 5],
        [s.frac_coords for i, s in enumerate(bulk) if i != 5],
        validate_proximity=False,
    )
    labels = [d.defect_name for d in cp.classify_defects(bulk, vac)]
    assert labels == ["v_Na_Oh"]


# ---------------------------------------------------------------------------
# export helpers round-trip
# ---------------------------------------------------------------------------


def test_export_round_trip(tmp_path):
    res = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    comp = res.complexes["Li_Na_Oh+Li_Na_Oh"]
    written = cp.io.export_pairs(res, str(tmp_path), write_all=True)
    assert set(written) == {"Li_Na_Oh+Li_Na_Oh"}
    # representative files: nn & nnn (2 copies each)
    for shell in ("nn", "nnn"):
        folder = os.path.join(str(tmp_path), "Li_Na_Oh+Li_Na_Oh", shell)
        assert os.path.exists(os.path.join(folder, "POSCAR"))
    # work dir holds every candidate
    work = os.path.join(str(tmp_path), "work_Li_Na_Oh+Li_Na_Oh")
    n = len([f for f in os.listdir(work) if f.startswith("POSCAR")])
    assert n == len(comp.candidates)
    # re-read matches in-memory structures
    comp = res.complexes["Li_Na_Oh+Li_Na_Oh"]
    for shell in comp.shells:
        path = os.path.join(str(tmp_path), "Li_Na_Oh+Li_Na_Oh", shell.label, "POSCAR")
        assert matcher.fit(shell.candidate.structure, Structure.from_file(path))


def test_csv_export(tmp_path):
    res = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")])
    path = os.path.join(str(tmp_path), "manifest.csv")
    cp.io.write_rows_csv(cp.pair_rows(res), path)
    import csv

    with open(path, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    assert rows[0]["combo"] == "Li_Na_Oh+Li_Na_Oh"
    assert rows[0]["shell"] == "nn"
