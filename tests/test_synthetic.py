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
from helpers import cubic_2x2x2, rocksalt_1x1x1, rocksalt_2x2x2

matcher = StructureMatcher(primitive_cell=False, attempt_supercell=False)


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
        comp.shells[i].distance <= comp.shells[i + 1].distance for i in range(len(comp.shells) - 1)
    )
    # shell assignments are consistent with the shell-grouping algorithm
    by_d = sorted(comp.candidates, key=lambda c: c.distance)
    assert [c.shell for c in by_d] == shell_groups([c.distance for c in by_d], 0.05)
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


def test_shells_validated_before_enumeration():
    # a typo must fail fast, not after the whole spglib run
    with pytest.raises(ValueError, match="shells must be one of"):
        cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="n")
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")])
    with pytest.raises(ValueError, match="pair_shell must be one of"):
        cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh", pair_shell="next")


def test_missing_shell_representatives_are_recorded():
    """A criterion with no candidate in the requested shell is reported."""
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    comp = pairs.complexes["Li_Na_Oh+Li_Na_Oh"]
    cand = comp.shells[0].candidate
    # one candidate only -> a single d_AC shell, so the second shell is empty
    tc = cp.TripleComplex(
        combo="Li_Na_Oh+Li_Na_Oh+Li_Na_Oh",
        pair=comp.combo,
        type_c=pairs.types[0],
        candidates=[
            cp.TripleCandidate(
                config_index=0,
                site_index=cand.site_index,
                structure=cand.structure,
                d_ac=2.0,
                d_bc=2.1,
                degeneracy=1,
            )
        ],
    )
    from codopex.triples import _select_criteria

    missing = _select_criteria(tc, tier=1, dist_tol=0.05)
    assert missing == ["nearA_farB", "nearB_farA"]
    assert tc.selection("nearA_farB") is None
    assert tc.selection("both_nn") is not None  # documented fallback


def test_third_restricts_the_added_dopant():
    bulk = rocksalt_2x2x2()
    pairs = cp.generate_pairs(bulk, [("Na", "Li"), ("Cl", "F")])
    all_third = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh")
    only_f = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh", third=[("Cl", "F")])
    assert len(only_f.complexes) < len(all_third.complexes)
    assert all(combo.endswith("F_Cl_Oh") for combo in only_f.complexes)


def test_dist_tol_overrides_the_pair_tolerance():
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    assert cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh").stats["dist_tol"] == pairs.dist_tol
    overridden = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh", dist_tol=0.2)
    assert overridden.stats["dist_tol"] == 0.2
    assert overridden.shells == "nn"


# ---------------------------------------------------------------------------
# deterministic output
# ---------------------------------------------------------------------------


def test_deterministic():
    r1 = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    r2 = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    assert [c.distance for c in r1.complexes["Li_Na_Oh+Li_Na_Oh"].candidates] == [
        c.distance for c in r2.complexes["Li_Na_Oh+Li_Na_Oh"].candidates
    ]


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


def test_classify_uses_global_assignment():
    """A greedy match would leave b2 and d2 unpaired (now an error); the
    optimal assignment pairs both Na sites with the two defect Na sites."""
    bulk = Structure(
        Lattice.cubic(1.0),
        ["Na", "Na", "Cl"],
        [[0.51, 0.5, 0.5], [0.50, 0.51, 0.5], [0.0, 0.0, 0.0]],
        validate_proximity=False,
    )
    defect = Structure(
        Lattice.cubic(1.0),
        ["Na", "Na", "Cl"],
        [[0.50, 0.5, 0.5], [0.60, 0.5, 0.5], [0.0, 0.0, 0.0]],
        validate_proximity=False,
    )
    # threshold 0.14 * (1/3)^(1/3) ~ 0.097 A: the direct greedy pick
    # (b1-d1 = 0.01 A) leaves b2 0.1005 A from d2, just beyond the
    # threshold, while the global pairing b1-d2 (0.09 A) + b2-d1 (0.01 A)
    # keeps both matches
    assert cp.classify_defects(bulk, defect, stol=0.14) == []


def test_classify_rejects_interstitial_atoms():
    """A site without a pristine counterpart is an error, not a defect label:
    codopex-generated structures only ever contain substitutions."""
    bulk = cubic_2x2x2()
    defect = bulk.copy()
    defect.append("Li", [0.25, 0.25, 0.25], validate_proximity=False)
    with pytest.raises(ValueError, match="no pristine counterpart"):
        cp.classify_defects(bulk, defect)


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


def test_export_triples_separates_shell_tiers(tmp_path):
    """nn and nnn runs must not share a folder: the criteria keep their names."""
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    combo = "Li_Na_Oh+Li_Na_Oh+Li_Na_Oh"
    nn = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh", shells="nn")
    nnn = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh", shells="nnn")
    cp.io.export_triples(nn, str(tmp_path))
    cp.io.export_triples(nnn, str(tmp_path))
    for shells in ("nn", "nnn"):
        folder = os.path.join(str(tmp_path), combo, shells, "nearA_farB")
        assert os.path.exists(os.path.join(folder, "POSCAR")), shells
    # and each tier keeps its own representative
    assert matcher.fit(
        nn.complexes[combo].selection("nearA_farB").structure,
        Structure.from_file(os.path.join(str(tmp_path), combo, "nn", "nearA_farB", "POSCAR")),
    )


def test_export_triples_write_all_counts_candidates(tmp_path):
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    tres = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh")
    written = cp.io.export_triples(tres, str(tmp_path), write_all=True)
    combo = "Li_Na_Oh+Li_Na_Oh+Li_Na_Oh"
    work = os.path.join(str(tmp_path), f"work_{combo}")
    assert len(os.listdir(work)) == len(tres.complexes[combo].candidates)
    # every written path is reported, so callers can count files
    assert str(tmp_path) in written[combo][0]
    assert len(written[combo]) == 6 + len(tres.complexes[combo].candidates)


def test_csv_export(tmp_path):
    res = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")])
    path = os.path.join(str(tmp_path), "nested", "manifest.csv")  # parent dir missing
    cp.io.write_rows_csv(cp.pair_rows(res), path)
    import csv

    with open(path, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    assert rows[0]["combo"] == "Li_Na_Oh+Li_Na_Oh"
    assert rows[0]["shell"] == "nn"
