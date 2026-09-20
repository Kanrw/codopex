"""Persistence (save/load) and batch triple generation.

The mechanics are tested on a small cubic crystal whose symmetry is easy to
reason about; structure-level fidelity of a save/load round trip is checked
with StructureMatcher.
"""

import pytest
from pymatgen.analysis.structure_matcher import StructureMatcher

import codopex as cp
from helpers import cubic_2x2x2

matcher = StructureMatcher(primitive_cell=False, attempt_supercell=False)


# ---------------------------------------------------------------------------
# save / load round trips
# ---------------------------------------------------------------------------


def test_pairs_save_load_round_trip(tmp_path):
    result = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    path = tmp_path / "pairs.json.gz"
    assert cp.io.save(result, str(path)) == str(path)

    loaded = cp.io.load(str(path))
    assert isinstance(loaded, cp.PairsResult)
    assert loaded.type_labels == result.type_labels
    assert loaded.reactions == result.reactions
    assert loaded.dist_tol == result.dist_tol
    assert loaded.stats == result.stats
    assert loaded.site_type_maps[0][0] is loaded.types[0]

    for label, base in result.bases.items():
        assert loaded.bases[label].dopant_site == base.dopant_site
        assert loaded.bases[label].degeneracy == base.degeneracy
        assert matcher.fit(base.structure, loaded.bases[label].structure)

    for combo, comp in result.complexes.items():
        other = loaded.complexes[combo]
        assert other.type_a is loaded.types[comp.type_a.index]
        assert other.type_b is loaded.types[comp.type_b.index]
        assert other.base_dopant_site == comp.base_dopant_site
        assert [s.label for s in other.shells] == [s.label for s in comp.shells]
        assert [c.distance for c in other.candidates] == [c.distance for c in comp.candidates]
        assert [c.shell for c in other.candidates] == [c.shell for c in comp.candidates]
        assert [c.degeneracy for c in other.candidates] == [c.degeneracy for c in comp.candidates]
        for ours, theirs in zip(other.shells, comp.shells, strict=True):
            assert theirs.distance == pytest.approx(ours.distance)
            assert matcher.fit(theirs.candidate.structure, ours.candidate.structure)


def test_triples_save_load_round_trip(tmp_path):
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    result = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh", shells="nnn")
    path = tmp_path / "triples.json"
    cp.io.save(result, str(path))

    loaded = cp.io.load(str(path))
    assert isinstance(loaded, cp.TriplesResult)
    assert loaded.pair == result.pair
    assert loaded.pair_shell == result.pair_shell
    assert loaded.stats == result.stats
    assert loaded.pairs.type_labels == pairs.type_labels

    for combo, tc in result.complexes.items():
        other = loaded.complexes[combo]
        assert other.type_c is loaded.pairs.types[tc.type_c.index]
        assert [s.criterion for s in other.selections] == [s.criterion for s in tc.selections]
        assert [c.d_ac for c in other.candidates] == [c.d_ac for c in tc.candidates]
        for ours, theirs in zip(other.selections, tc.selections, strict=True):
            assert ours.candidate is not None
            assert theirs.candidate is not None
            assert theirs.candidate.d_ac == pytest.approx(ours.candidate.d_ac)
            assert theirs.candidate.d_bc == pytest.approx(ours.candidate.d_bc)
            assert theirs.candidate.degeneracy == ours.candidate.degeneracy
            assert matcher.fit(theirs.candidate.structure, ours.candidate.structure)


def test_batch_save_load_round_trip(tmp_path):
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    batch = cp.generate_all_triples(pairs)
    path = tmp_path / "batch.json.gz"
    cp.io.save(batch, str(path))

    loaded = cp.io.load(str(path))
    assert isinstance(loaded, cp.AllTriplesResult)
    assert list(loaded.results) == list(batch.results)
    assert loaded.skipped == batch.skipped
    for parent in batch:
        assert set(loaded[parent].complexes) == set(batch[parent].complexes)
        assert loaded[parent].pair == parent


def test_load_rejects_foreign_file(tmp_path):
    foreign = tmp_path / "foreign.json"
    foreign.write_text("{}")
    with pytest.raises(ValueError, match="not a codopex result file"):
        cp.io.load(str(foreign))


def test_save_rejects_other_objects(tmp_path):
    with pytest.raises(TypeError, match="PairsResult"):
        cp.io.save("not a result", str(tmp_path / "x.json"))


# ---------------------------------------------------------------------------
# batch triples
# ---------------------------------------------------------------------------


def test_generate_all_triples_matches_single_calls():
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    batch = cp.generate_all_triples(pairs)
    assert list(batch.results) == sorted(pairs.complexes)
    single = cp.generate_triples(pairs, "Li_Na_Oh+Li_Na_Oh")
    assert set(batch["Li_Na_Oh+Li_Na_Oh"].complexes) == set(single.complexes)
    assert batch.skipped == []


def test_generate_all_triples_skips_missing_shell():
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")])  # nn only
    batch = cp.generate_all_triples(pairs, pair_shell="nnn")
    assert len(batch) == 0
    assert batch.skipped == ["Li_Na_Oh+Li_Na_Oh"]


def test_generate_all_triples_rejects_unknown_pair():
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")])
    with pytest.raises(ValueError, match="not among generated combos"):
        cp.generate_all_triples(pairs, pairs=["K_Na_Oh+K_Na_Oh"])
