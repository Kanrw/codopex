"""Shared helpers of the pair/triple pipelines.

* pristine-site bookkeeping: every site of the bulk keeps a stable index;
  substitutions only change species, never order, so the pristine identity
  of any site in any later structure is its original bulk index;
* defect-type construction: a defect type is a reaction (host, dopant)
  applied to one symmetry class of host sites, labelled with the pristine
  site symmetry (e.g. ``Cu_Pb_C3v``);
* distance-shell clustering (the notebook's shell algorithm).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from pymatgen.core import Structure

from codopex.symmetry import site_table
from codopex.types import DefectType, Reaction

SHELL_LABELS = {0: "nn", 1: "nnn"}  # first and second distance shell

SiteKey = Tuple[str, str, str]  # (host element, site symmetry, wyckoff tag)


def load_bulk(bulk) -> Structure:
    """Accept a Structure or a structure file path."""
    if isinstance(bulk, Structure):
        return bulk
    return Structure.from_file(bulk)


def dedupe_reactions(reactions: List[Reaction]) -> List[Reaction]:
    """Remove duplicate (host, dopant) entries, keeping first order."""
    seen: set = set()
    out: List[Reaction] = []
    for r in reactions:
        if r not in seen:
            seen.add(r)
            out.append(r)
    return out


def build_types(
    bulk: Structure,
    reactions: List[Reaction],
    site_symprec: float = 0.01,
    type_order: Optional[List[str]] = None,
) -> Tuple[List[DefectType], List[Dict[int, Optional[DefectType]]]]:
    """Discover defect types from the pristine bulk.

    A "type" is one (host, dopant) reaction applied to one symmetry class of
    host sites.  Host-site classes are keyed by (site symmetry, wyckoff), so
    two inequivalent classes that share a site-symmetry label stay distinct
    types (their labels then carry a wyckoff disambiguator).

    Canonical order: reactions in the order given; within a reaction, types
    in order of first occurrence along the bulk site list.  ``type_order``
    overrides the final ordering by type label (e.g. to reproduce the hand
    ordering of an earlier study).

    Returns
    -------
    types : list of DefectType in canonical order
    type_by_reaction : per reaction index, a site-index -> DefectType map
        (None for sites that are not host sites of that reaction).
    """
    reactions = dedupe_reactions(reactions)
    info = site_table(bulk, symprec=site_symprec)
    host_set = {host for host, _ in reactions}

    # per reaction: site key -> first bulk index where it occurs
    first_site: List[Dict[SiteKey, int]] = []
    for host, _dop in reactions:
        seen: Dict[SiteKey, int] = {}
        for s, it in enumerate(info):
            if it.element == host:
                key = (host, it.site_symmetry, it.wyckoff)
                if key not in seen:
                    seen[key] = s
        first_site.append(seen)

    # raw entries in canonical order (reaction order, then first-site order)
    raw: List[Tuple[int, int, SiteKey, str]] = []  # (r, site, key, base label)
    for r, (host, dop) in enumerate(reactions):
        for key, s in first_site[r].items():  # dict is insertion-ordered
            raw.append((r, s, key, f"{dop}_{host}_{key[1]}"))

    # disambiguate duplicated plain labels (same host+dopant+symmetry on two
    # different Wyckoff classes) by appending the wyckoff tag
    base_count: Dict[Tuple[int, str], int] = {}
    for r, _s, _key, label in raw:
        base_count[(r, label)] = base_count.get((r, label), 0) + 1
    labels: List[str] = []
    for r, _s, key, label in raw:
        if base_count[(r, label)] > 1:
            labels.append(f"{label}_{key[2]}")
        else:
            labels.append(label)

    types = [
        DefectType(
            host=reactions[r][0],
            dopant=reactions[r][1],
            site_symmetry=key[1],
            site_tag=key[2],
            index=i,
            label=labels[i],
        )
        for i, (r, _s, key, _label) in enumerate(raw)
    ]
    if type_order is not None:
        types = _apply_type_order(types, type_order)

    # per-reaction site -> type map (None for non-host sites)
    type_by_reaction: List[Dict[int, Optional[DefectType]]] = []
    for r, (host, dop) in enumerate(reactions):
        matching = [t for t in types if t.host == host and t.dopant == dop]
        mapping: Dict[int, Optional[DefectType]] = {}
        for s, it in enumerate(info):
            if it.element != host:
                mapping[s] = None
                continue
            key = (host, it.site_symmetry, it.wyckoff)
            hit = [t for t in matching if t.site_symmetry == key[1] and t.site_tag == key[2]]
            mapping[s] = hit[0] if len(hit) == 1 else None
        type_by_reaction.append(mapping)
    return types, type_by_reaction


def _apply_type_order(types: List[DefectType], type_order: List[str]) -> List[DefectType]:
    by_label = {t.label: t for t in types}
    if set(type_order) != set(by_label):
        raise ValueError(
            f"type_order must list every discovered type exactly once; "
            f"got {type_order}, expected {sorted(by_label)}"
        )
    return [
        DefectType(
            host=t.host, dopant=t.dopant, site_symmetry=t.site_symmetry,
            site_tag=t.site_tag, index=i, label=t.label,
        )
        for i, t in enumerate(by_label[l] for l in type_order)
    ]


def shell_groups(distances: Sequence[float], dist_tol: float) -> List[int]:
    """Assign 0-based shell-group indices to *sorted* distances.

    Reproduces the notebook algorithm: a new shell starts when a distance
    exceeds the start of the previous shell by more than ``dist_tol``
    Angstrom.
    """
    groups: List[int] = []
    starts: List[float] = []
    for d in distances:
        if not starts or d - starts[-1] > dist_tol:
            starts.append(d)
        groups.append(len(starts) - 1)
    return groups
