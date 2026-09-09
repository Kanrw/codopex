"""Bulk/defect comparison: label every defect in a defective structure.

This is a structure-level (in-memory) port of the ``classify_all_defects``
helper from the original notebooks: match sites between the pristine bulk and
a defective structure, then name each mismatch as a vacancy, a substitution
or an interstitial, with the pristine site symmetry attached:

* ``v_{host}_{sym}``         vacancy
* ``{dop}_{host}_{sym}``     substitution
* ``{el}_i_{sym}``           interstitial
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from pymatgen.core import Structure

from codopex.symmetry import format_site_symmetry, site_table


def _stol_threshold(structure: Structure, stol: float) -> float:
    return stol * (structure.volume / len(structure)) ** (1 / 3)


def _match_same_species(
    bulk: Structure, defect: Structure, threshold: float
) -> Tuple[Dict[int, Optional[int]], set]:
    lattice = bulk.lattice
    bulk_to_defect: Dict[int, Optional[int]] = {}
    defect_matched: set = set()

    for i, b_site in enumerate(bulk):
        best_j, min_dist = None, threshold + 1.0
        for j, d_site in enumerate(defect):
            if j in defect_matched:
                continue
            if d_site.species_string != b_site.species_string:
                continue
            dist = lattice.get_all_distances([b_site.frac_coords], [d_site.frac_coords])[0, 0]
            if dist < min_dist:
                min_dist, best_j = dist, j
        if best_j is not None and min_dist < threshold:
            bulk_to_defect[i] = best_j
            defect_matched.add(best_j)
        else:
            bulk_to_defect[i] = None
    return bulk_to_defect, defect_matched


def _match_cross_species(
    unmatched_bulk: List[int],
    unmatched_defect: List[int],
    bulk: Structure,
    defect: Structure,
    threshold: float,
) -> Tuple[List[Tuple[int, int]], List[int], List[int]]:
    lattice = bulk.lattice
    used_defect: set = set()
    sub_pairs: List[Tuple[int, int]] = []

    for i in unmatched_bulk:
        best_j, min_dist = None, threshold + 1.0
        for j in unmatched_defect:
            if j in used_defect:
                continue
            dist = lattice.get_all_distances([bulk[i].frac_coords], [defect[j].frac_coords])[0, 0]
            if dist < min_dist:
                min_dist, best_j = dist, j
        if best_j is not None and min_dist < threshold:
            sub_pairs.append((i, best_j))
            used_defect.add(best_j)

    matched_bulk = {i for i, _ in sub_pairs}
    matched_defect = {j for _, j in sub_pairs}
    remaining_bulk = [i for i in unmatched_bulk if i not in matched_bulk]
    remaining_defect = [j for j in unmatched_defect if j not in matched_defect]
    return sub_pairs, remaining_bulk, remaining_defect


def _interstitial_site_symmetry(
    defect: Structure, frac_coords, symprec: float = 0.01
) -> str:
    temp = defect.copy()
    temp.append("He", frac_coords, validate_proximity=False)
    try:
        sym = site_table(temp, symprec=symprec)[-1].site_symmetry
    except Exception:  # noqa: BLE001 - symmetry may fail on weird geometries
        sym = "C1"
    return sym


@dataclass
class DefectRecord:
    """One identified defect in a defective structure."""

    defect_type: str  # "vacancy" | "substitution" | "interstitial"
    defect_name: str  # canonical label, e.g. "Cu_Pb_C3v"


def classify_defects(
    bulk: Structure,
    defect: Structure,
    symprec: float = 0.01,
    stol: float = 0.5,
) -> List[DefectRecord]:
    """Identify all defects in ``defect`` relative to the pristine ``bulk``.

    Both structures must share the same lattice.  The reported names use the
    *pristine* site symmetry of the host position (from ``bulk``), matching
    the convention ``{dop}_{host}_{sym}`` used throughout the co-doping
    workflow.  This routine is mainly useful for validation and for labeling
    externally generated structures; the codopex pipeline itself labels
    configurations directly from the pristine site table.
    """
    bulk_info = site_table(bulk, symprec=symprec)
    threshold = _stol_threshold(bulk, stol)

    bulk_to_defect, defect_matched = _match_same_species(bulk, defect, threshold)
    unmatched_bulk = [i for i, j in bulk_to_defect.items() if j is None]
    unmatched_defect = [j for j in range(len(defect)) if j not in defect_matched]
    sub_pairs, vacancy_bulk, interstitial_defect = _match_cross_species(
        unmatched_bulk, unmatched_defect, bulk, defect, threshold
    )

    records: List[DefectRecord] = []

    for i in vacancy_bulk:
        info = bulk_info[i]
        records.append(
            DefectRecord("vacancy", f"v_{info.element}_{info.site_symmetry}")
        )

    for bulk_idx, defect_idx in sub_pairs:
        info = bulk_info[bulk_idx]
        records.append(
            DefectRecord(
                "substitution",
                f"{defect[defect_idx].species_string}_{info.element}_{info.site_symmetry}",
            )
        )

    for j in interstitial_defect:
        site = defect[j]
        sym = _interstitial_site_symmetry(defect, site.frac_coords, symprec=symprec)
        records.append(
            DefectRecord("interstitial", f"{site.species_string}_i_{sym}")
        )

    return records
