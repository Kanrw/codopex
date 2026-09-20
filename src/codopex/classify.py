"""Bulk/defect comparison: label every defect in a defective structure.

This is a structure-level (in-memory) port of the ``classify_all_defects``
helper from the original notebooks: match sites between the pristine bulk and
a defective structure, then name each mismatch with the pristine site
symmetry attached:

* ``v_{host}_{sym}``         vacancy
* ``{dop}_{host}_{sym}``     substitution

Sites are paired by a globally optimal one-to-one assignment (as many
matches as possible below the distance threshold, then the smallest total
displacement), computed from a single periodic distance matrix.  A site of
the defective structure without a pristine counterpart (an interstitial
atom) raises ``ValueError``: every site of a generated structure must
correspond to a pristine site, so codopex only produces substitutions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pymatgen.core import Structure
from scipy.optimize import linear_sum_assignment

from codopex.symmetry import site_table

# Cost of a pair that must not be matched (wrong species or beyond the
# threshold).  Larger than any possible real total displacement, so the
# assignment first maximizes the number of real matches, then minimizes
# their total distance.
_NO_MATCH_COST = 1e6


def _stol_threshold(structure: Structure, stol: float) -> float:
    return stol * (structure.volume / len(structure)) ** (1 / 3)


def _distances(bulk: Structure, defect: Structure) -> np.ndarray:
    """Periodic distance matrix between the sites of two structures."""
    return bulk.lattice.get_all_distances(bulk.frac_coords, defect.frac_coords)


def _match_pairs(
    distances: np.ndarray, allowed: np.ndarray, threshold: float
) -> list[tuple[int, int]]:
    """Best one-to-one matching of allowed site pairs below ``threshold``.

    Pairs that are not allowed or exceed the threshold are dropped after the
    assignment, leaving the maximum possible number of matches.
    """
    cost = np.where(allowed & (distances < threshold), distances, _NO_MATCH_COST)
    rows, cols = linear_sum_assignment(cost)
    return [
        (int(i), int(j))
        for i, j in zip(rows, cols, strict=True)
        if allowed[i, j] and distances[i, j] < threshold
    ]


def _match_same_species(
    bulk: Structure, defect: Structure, distances: np.ndarray, threshold: float
) -> tuple[dict[int, int | None], set[int]]:
    """Match sites of equal species (defect-free sites of the bulk)."""
    bulk_species = np.array([site.species_string for site in bulk])
    defect_species = np.array([site.species_string for site in defect])
    allowed = bulk_species[:, None] == defect_species[None, :]

    bulk_to_defect: dict[int, int | None] = {i: None for i in range(len(bulk))}
    defect_matched: set[int] = set()
    for i, j in _match_pairs(distances, allowed, threshold):
        bulk_to_defect[i] = j
        defect_matched.add(j)
    return bulk_to_defect, defect_matched


def _match_cross_species(
    unmatched_bulk: list[int],
    unmatched_defect: list[int],
    distances: np.ndarray,
    threshold: float,
) -> tuple[list[tuple[int, int]], list[int], list[int]]:
    """Match the remaining sites against each other (substitutions)."""
    if not unmatched_bulk or not unmatched_defect:
        return [], unmatched_bulk, unmatched_defect

    sub_distances = distances[np.ix_(unmatched_bulk, unmatched_defect)]
    pairs = _match_pairs(sub_distances, np.ones_like(sub_distances, dtype=bool), threshold)

    sub_pairs = [(unmatched_bulk[i], unmatched_defect[j]) for i, j in pairs]
    matched_bulk = {i for i, _ in pairs}
    matched_defect = {j for _, j in pairs}
    remaining_bulk = [site for k, site in enumerate(unmatched_bulk) if k not in matched_bulk]
    remaining_defect = [site for k, site in enumerate(unmatched_defect) if k not in matched_defect]
    return sub_pairs, remaining_bulk, remaining_defect


@dataclass
class DefectRecord:
    """One identified defect in a defective structure."""

    defect_type: str  # "vacancy" | "substitution"
    defect_name: str  # canonical label, e.g. "Cu_Pb_C3v"


def classify_defects(
    bulk: Structure,
    defect: Structure,
    symprec: float = 0.01,
    stol: float = 0.5,
) -> list[DefectRecord]:
    """Identify all defects in ``defect`` relative to the pristine ``bulk``.

    Both structures must share the same lattice.  The reported names use the
    *pristine* site symmetry of the host position (from ``bulk``), matching
    the convention ``{dop}_{host}_{sym}`` used throughout the co-doping
    workflow.  This routine is mainly useful for validation and for labeling
    externally generated structures; the codopex pipeline itself labels
    configurations directly from the pristine site table.

    Raises
    ------
    ValueError
        If ``defect`` has a site that does not correspond to any pristine
        site (an interstitial atom), which a codopex-generated structure
        never does.
    """
    bulk_info = site_table(bulk, symprec=symprec)
    threshold = _stol_threshold(bulk, stol)
    distances = _distances(bulk, defect)

    bulk_to_defect, defect_matched = _match_same_species(bulk, defect, distances, threshold)
    unmatched_bulk = [i for i, j in bulk_to_defect.items() if j is None]
    unmatched_defect = [j for j in range(len(defect)) if j not in defect_matched]
    sub_pairs, vacancy_bulk, extra_defect_sites = _match_cross_species(
        unmatched_bulk, unmatched_defect, distances, threshold
    )
    if extra_defect_sites:
        raise ValueError(
            f"{len(extra_defect_sites)} site(s) of the defective structure have no "
            "pristine counterpart; interstitial atoms are not supported "
            "(codopex only generates substitutional complexes)"
        )

    records: list[DefectRecord] = []

    for i in vacancy_bulk:
        info = bulk_info[i]
        records.append(DefectRecord("vacancy", f"v_{info.element}_{info.site_symmetry}"))

    for bulk_idx, defect_idx in sub_pairs:
        info = bulk_info[bulk_idx]
        records.append(
            DefectRecord(
                "substitution",
                f"{defect[defect_idx].species_string}_{info.element}_{info.site_symmetry}",
            )
        )

    return records
