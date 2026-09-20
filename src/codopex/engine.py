"""Symmetry-inequivalent substitution enumeration.

The engine replaces the SAGAR/pyvaspflow ``DefectMaker`` of the original
notebooks with a pure pymatgen/spglib implementation.

For a *base* structure (pristine bulk, or bulk already carrying some
dopants) and a substitution reaction ``host -> dopant``, the inequivalent
single-substitution configurations are exactly the orbits of the host atoms
under the space group of the base structure: two host sites are equivalent
iff a symmetry operation of the base maps one onto the other while leaving
the existing decoration intact.  Orbit representatives are the members of
each orbit with the smallest site index, giving a deterministic order.

Site order is never changed by the engine: ``base`` and all returned
structures have identical site lists, only species differ at the substituted
indices.  No files are written and no directory is changed.
"""

from __future__ import annotations

from collections.abc import Sequence

from pymatgen.core import Structure

from codopex.symmetry import orbits_of
from codopex.types import SubstitutionConfig


def enumerate_substitutions(
    base: Structure,
    host: str,
    dopant: str,
    symprec: float = 1e-3,
    equiv: Sequence[int] | None = None,
) -> list[SubstitutionConfig]:
    """Enumerate all inequivalent ways to substitute one ``host`` by ``dopant``.

    Parameters
    ----------
    base : Structure
        The starting structure (pristine bulk or already doped).
    host : str
        Element of the sites to be substituted.
    dopant : str
        Substituting element.
    symprec : float
        Symmetry tolerance (Angstrom) for the spglib space-group search.
        The default 1e-3 matches the SAGAR tolerance used by the original
        DefectMaker workflow.
    equiv : optional sequence of int
        Precomputed orbit representatives of ``base`` (``equivalent_atoms``),
        to reuse one spglib search across several reactions on the same base.

    Returns
    -------
    List[SubstitutionConfig]
        One entry per inequivalent configuration, in deterministic orbit
        order (by smallest host-site index).
    """
    if host not in {site.species_string for site in base}:
        raise ValueError(f"host element {host!r} not present in base structure")

    candidates: list[SubstitutionConfig] = []
    for orbit in orbits_of(base, species=[host], symprec=symprec, equiv=equiv):
        rep = orbit[0]  # smallest index = deterministic representative
        substituted = base.copy()
        substituted.replace(rep, dopant)
        candidates.append(
            SubstitutionConfig(
                site_index=rep,
                structure=substituted,
                degeneracy=len(orbit),
                config_index=len(candidates),
            )
        )
    return candidates


def periodic_distance(structure: Structure, i: int, j: int) -> float:
    """Minimum periodic distance (Angstrom) between two sites of a structure."""
    return float(
        structure.lattice.get_distance_and_image(
            structure[i].frac_coords, structure[j].frac_coords
        )[0]
    )
