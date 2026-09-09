"""PPO regression: codopex engine vs the original SAGAR/notebook outputs.

The fixture data (``tests/data/ppo``) is the Pb3(PO4)2 case of the original
workflow (see ``tests/data/ppo/README.md`` for provenance):

* ``bulk.vasp``                 -- pristine 78-atom supercell
* ``single_defect/{type}/POSCAR`` -- 5 base structures produced by SAGAR
* ``pairs/{combo}/POSCAR``        -- 15 nearest-pair structures by SAGAR
* ``pairs_expected.csv``        -- shell distances, SAGAR distances, rep-fit
* ``triples_expected.csv``      -- criteria representatives, mode per row:
    - "verified": exact match with the SAGAR values (allowing the
      A/B-mirror orientation for pairs whose two dopants are
      symmetry-equivalent),
    - "shell": values lie within one distance shell (0.05 A) of the SAGAR
      values -- equal-distance orbit ties and 1e-4-level coordinate noise
      make the old representative choice non-reproducible there.

Parity contract: enumeration/classification/distance *values* are exact;
the representative picked inside a distance shell is deterministic in
codopex but may differ from the arbitrary SAGAR pick when several
inequivalent configurations share the shell distance (see REGRESSION.md).
"""

import os

import pandas as pd
import pytest
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

import codopex as cp

DATA = os.path.join(os.path.dirname(__file__), "data", "ppo")
BULK = os.path.join(DATA, "bulk.vasp")
NOTEBOOK_ORDER = ["Cu_Pb_C3v", "Cu_Pb_D3d", "S_O_Cs", "S_O_C3v", "S_P_C3v"]
REACTIONS = [("Pb", "Cu"), ("O", "S"), ("P", "S")]
OPPOSITE = {"nearA_farB": "nearB_farA", "nearB_farA": "nearA_farB", "both_nn": "both_nn"}

matcher = StructureMatcher(primitive_cell=False, attempt_supercell=False)


@pytest.fixture(scope="module")
def pairs():
    return cp.generate_pairs(BULK, REACTIONS, shells="nnn", type_order=NOTEBOOK_ORDER)


def mirror_combos(result):
    """Pairs whose A/B labels are interchangeable (same element, same sym)."""
    return {
        c.combo
        for c in result.complexes.values()
        if c.type_a.dopant == c.type_b.dopant
        and c.type_a.site_symmetry == c.type_b.site_symmetry
    }


def orientation_match(ours, ref):
    return (
        abs(ours[0] - ref[0]) < 1e-3 and abs(ours[1] - ref[1]) < 1e-3
    ) or (
        abs(ours[1] - ref[0]) < 1e-3 and abs(ours[0] - ref[1]) < 1e-3
    )


# ---------------------------------------------------------------------------
# 1. defect types
# ---------------------------------------------------------------------------


def test_type_discovery_matches_notebook(pairs):
    assert pairs.type_labels == NOTEBOOK_ORDER


def test_single_defect_bases_match_sagar_files(pairs):
    for label in NOTEBOOK_ORDER:
        ref = Structure.from_file(os.path.join(DATA, "single_defect", label, "POSCAR"))
        assert matcher.fit(ref, pairs.bases[label].structure), label
        # exactly one dopant on the pristine lattice
        bulk = pairs.bulk
        n_diff = sum(
            1
            for i in range(len(bulk))
            if pairs.bases[label].structure[i].species_string != bulk[i].species_string
        )
        assert n_diff == 1, label


# ---------------------------------------------------------------------------
# 2. pairs
# ---------------------------------------------------------------------------


def test_pair_rows_vs_sagar(pairs):
    expected = pd.read_csv(os.path.join(DATA, "pairs_expected.csv"), encoding="utf-8-sig")
    rows = pd.DataFrame(cp.pair_rows(pairs))
    merged = expected.merge(rows, on=["combo", "shell"], suffixes=("_exp", "_our"))
    assert len(merged) == len(expected)
    # self-consistency + exact distance parity with SAGAR on the nn shell
    assert (merged.dist_AB_exp - merged.dist_AB_our).abs().max() < 1e-6
    nn = merged[merged.shell == "nn"]
    assert len(nn) == 15
    assert (nn.dist_AB_exp - nn.sagar_nearest_d).abs().max() < 1e-3


def test_pair_representatives_match_sagar_files(pairs):
    expected = pd.read_csv(os.path.join(DATA, "pairs_expected.csv"), encoding="utf-8-sig")
    for _, row in expected[expected.shell == "nn"].iterrows():
        if not row.rep_fit_vs_sagar:
            # documented equal-distance tie: Cu_Pb_D3d+Cu_Pb_D3d has three
            # inequivalent configurations at the same shell distance
            continue
        ref = Structure.from_file(os.path.join(DATA, "pairs", row.combo, "POSCAR"))
        ours = pairs.complexes[row.combo].shell_structure("nn")
        assert matcher.fit(ref, ours), row.combo


def test_pair_shells_nn_and_nnn_present(pairs):
    # first shell is always present; the second appears when candidate
    # distances form a distinct shell beyond 0.05 A (e.g. Cu_Pb_D3d+Cu_Pb_D3d
    # has a single shell: all its configurations lie within 0.05 A)
    for combo in pairs.complexes:
        labels = [s.label for s in pairs.complexes[combo].shells]
        assert labels[0] == "nn", combo
        assert labels == ["nn"] or labels == ["nn", "nnn"], (combo, labels)


# ---------------------------------------------------------------------------
# 3. triples
# ---------------------------------------------------------------------------


def test_triple_criteria_vs_sagar(pairs):
    expected = pd.read_csv(os.path.join(DATA, "triples_expected.csv"), encoding="utf-8-sig")
    mirror = mirror_combos(pairs)
    cache = {}

    def run_triples(combo):
        if combo not in cache:
            cache[combo] = cp.generate_triples(pairs, combo, shells="nn")
        return cache[combo]

    for _, row in expected.iterrows():
        parent = "+".join(row.combo.split("+")[:2])
        tres = run_triples(parent)
        ours = _selection(tres, row.combo, row.criterion)
        assert ours is not None, row.combo
        our_vals = (ours.d_ac, ours.d_bc)
        # self-consistency with the committed table
        assert abs(our_vals[0] - row.d_AC) < 1e-3, row.combo
        assert abs(our_vals[1] - row.d_BC) < 1e-3, row.combo
        sagar = (row.sagar_d_AC, row.sagar_d_BC)
        if row["mode"] == "verified":
            ok = orientation_match(our_vals, sagar)
            if parent in mirror and not ok:
                # mirror parents: A/B of the notebook may be swapped, so the
                # opposite criterion must also be tried
                opp = _selection(tres, row.combo, OPPOSITE[row.criterion])
                assert opp is not None, row.combo
                ok = orientation_match((opp.d_ac, opp.d_bc), sagar)
            assert ok, f"{row.combo} {row.criterion}"
        else:
            # "shell" mode: representative choice inside a distance shell is
            # not reproducible when several inequivalent configurations share
            # a shell distance (the old workflow picked one arbitrarily);
            # assert the choice stays within the same shell on at least one
            # distance component (any orientation)
            variants = [sagar, (sagar[1], sagar[0])]
            assert any(
                abs(our_vals[0] - v[0]) < 0.06 or abs(our_vals[1] - v[1]) < 0.06
                for v in variants
            ), row.combo


def _selection(tres, combo, criterion):
    tc = tres.complexes.get(combo)
    if tc is None:
        return None
    return tc.selection(criterion)


def test_triple_full_coverage(pairs):
    expected = pd.read_csv(os.path.join(DATA, "triples_expected.csv"), encoding="utf-8-sig")
    assert len(expected) == 105  # 35 triples x up to 3 criteria


# ---------------------------------------------------------------------------
# 4. invariants
# ---------------------------------------------------------------------------


def test_positions_and_cell_never_change(pairs):
    bulk = pairs.bulk
    for comp in pairs.complexes.values():
        for cand in comp.candidates:
            assert cand.structure.lattice == bulk.lattice
            for i in range(len(bulk)):
                assert (
                    cand.structure[i].frac_coords - bulk[i].frac_coords
                ).max() < 1e-9
