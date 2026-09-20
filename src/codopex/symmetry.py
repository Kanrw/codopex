"""Site-symmetry analysis helpers (pure pymatgen/spglib).

Every function takes and returns pymatgen objects; no file I/O happens here.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Hermann-Mauguin -> Schoenflies site-symmetry labels (copy of the mapping
# used in the original classification notebooks).
_HM_TO_SCHOENFLIES = {
    "1": "C1",
    "-1": "Ci",
    "2": "C2",
    "m": "Cs",
    "2/m": "C2h",
    "222": "D2",
    "mm2": "C2v",
    "mmm": "D2h",
    "4": "C4",
    "-4": "S4",
    "4/m": "C4h",
    "422": "D4",
    "4mm": "C4v",
    "-42m": "D2d",
    "4/mmm": "D4h",
    "3": "C3",
    "-3": "C3i",
    "32": "D3",
    "3m": "C3v",
    "-3m": "D3d",
    "6": "C6",
    "-6": "C3h",
    "6/m": "C6h",
    "622": "D6",
    "6mm": "C6v",
    "-6m2": "D3h",
    "6/mmm": "D6h",
    "23": "T",
    "m-3": "Th",
    "432": "O",
    "-43m": "Td",
    "m-3m": "Oh",
}


def format_site_symmetry(raw_symbol: str) -> str:
    """Normalize an spglib site-symmetry symbol to Schoenflies notation."""
    hm = raw_symbol.strip().replace(" ", "").replace(".", "")
    return _HM_TO_SCHOENFLIES.get(hm, hm)


@dataclass
class SiteInfo:
    """Pristine-site properties of one atom."""

    element: str
    wyckoff: str  # e.g. "6c"
    site_symmetry: str  # Schoenflies label, e.g. "C3v"


def _dataset(structure: Structure, symprec: float):
    """spglib symmetry dataset, tolerant of pymatgen API changes."""
    sga = SpacegroupAnalyzer(structure, symprec=symprec)
    data = sga.get_symmetry_dataset()
    if isinstance(data, dict):
        return data
    # newer pymatgen returns a SpglibDataset-like object with attributes
    return {
        "wyckoffs": data.wyckoffs,
        "site_symmetry_symbols": data.site_symmetry_symbols,
        "equivalent_atoms": data.equivalent_atoms,
    }


def site_table(structure: Structure, symprec: float = 0.01) -> list[SiteInfo]:
    """Per-site symmetry information for ``structure``.

    ``wyckoff`` combines the multiplicity of the site's orbit (within this
    cell) with the spglib Wyckoff letter.  Site symmetry is converted from
    Hermann-Mauguin to Schoenflies notation.
    """
    data = _dataset(structure, symprec=symprec)
    wyckoffs = data["wyckoffs"]
    sym_symbols = data["site_symmetry_symbols"]
    equiv = data["equivalent_atoms"]

    n = len(structure)
    if len(wyckoffs) != n or len(sym_symbols) != n or len(equiv) != n:
        raise ValueError(
            "spglib dataset inconsistent with structure size; "
            "check the structure is a valid crystal with PBC"
        )

    multiplicity: dict[int, int] = {}
    for rep in equiv:
        multiplicity[rep] = multiplicity.get(rep, 0) + 1

    info: list[SiteInfo] = []
    for i, site in enumerate(structure):
        rep = int(equiv[i])
        info.append(
            SiteInfo(
                element=site.species_string,
                wyckoff=f"{multiplicity[rep]}{wyckoffs[i]}",
                site_symmetry=format_site_symmetry(sym_symbols[i]),
            )
        )
    return info


def equivalent_atoms(structure: Structure, symprec: float = 1e-3) -> list[int]:
    """Orbit representative index per site (spglib ``equivalent_atoms``)."""
    return [int(i) for i in _dataset(structure, symprec=symprec)["equivalent_atoms"]]


def orbits_of(
    structure: Structure,
    species: Sequence[str] | None = None,
    symprec: float = 1e-3,
) -> list[list[int]]:
    """Partition site indices into symmetry orbits.

    Only sites whose species are in ``species`` (or all sites if None) are
    kept.  Orbits are returned sorted by their smallest site index, which
    defines the deterministic configuration order of the enumeration engine.
    """
    equiv = equivalent_atoms(structure, symprec=symprec)
    grouped: dict[int, list[int]] = {}
    for i, rep in enumerate(equiv):
        if species is not None and structure[i].species_string not in species:
            continue
        grouped.setdefault(int(rep), []).append(i)
    orbits = [sorted(members) for members in grouped.values()]
    orbits.sort(key=lambda o: o[0])
    return orbits


def fractional_distance_to(frac: np.ndarray, point: Sequence[float]) -> float:
    """Minimum-image distance between fractional coordinates (mod 1 metric).

    Used for the "defect closest to the supercell center" rule of the
    original notebooks; the metric is deliberately the same (norm of
    periodic fractional differences), not a Cartesian distance.
    """
    diff = np.asarray(frac, dtype=float) - np.asarray(point, dtype=float)
    return float(np.linalg.norm((diff + 0.5) % 1.0 - 0.5))
