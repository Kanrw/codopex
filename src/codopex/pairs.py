"""Phase-1 pipeline: single-defect types and symmetry-inequivalent defect pairs.

Faithful, structure-level port of the ``gen_codope_pairs`` notebook logic:

1. discover defect types per (host, dopant) reaction (pristine site-symmetry
   classes);
2. build one base structure per type: the inequivalent single-dopant
   configuration whose dopant sits closest to the supercell center;
3. for every pair combo (with repetition) enumerate all inequivalent
   configurations of the second dopant around the base structure of the
   first type, group them into distance shells (tolerance ``dist_tol``), and
   keep one representative per included shell: "nn" (first shell, default)
   or "nn" + "nnn" (first two shells).

All structures share the site order of the pristine bulk; no files are
written here.
"""

from __future__ import annotations

from collections.abc import Sequence

from codopex import engine
from codopex.pipeline import (
    SHELL_LABELS,
    build_types,
    dedupe_reactions,
    load_bulk,
    shell_groups,
)
from codopex.symmetry import fractional_distance_to, site_table
from codopex.types import (
    PairCandidate,
    PairComplex,
    PairShell,
    PairsResult,
    Reaction,
    SingleDefectBase,
)

CENTER = (0.5, 0.5, 0.5)


def generate_pairs(
    bulk,
    reactions: list[Reaction],
    shells: str = "nn",
    dist_tol: float = 0.05,
    enum_symprec: float = 1e-3,
    site_symprec: float = 0.01,
    center: Sequence[float] = CENTER,
    type_order: list[str] | None = None,
) -> PairsResult:
    """Generate nearest (and optionally next-nearest) defect-pair complexes.

    Parameters
    ----------
    bulk : Structure or str
        Pristine defect-free supercell (structure object or path to a
        structure file readable by pymatgen).
    reactions : list of (host, dopant)
        Substitution reactions in canonical order, e.g.
        ``[("Pb", "Cu"), ("O", "S"), ("P", "S")]``.
    shells : "nn" or "nnn"
        Which distance shells to keep as representatives per combo:
        "nn" = first shell only (default), "nnn" = first and second shell.
    dist_tol : float
        Distance-shell tolerance in Angstrom.
    enum_symprec : float
        spglib symmetry tolerance for configuration enumeration
        (1e-3 matches the SAGAR engine of the original workflow).
    site_symprec : float
        Symmetry tolerance for pristine site-symmetry labels.
    center : sequence of 3 floats
        Supercell center used by the base-selection rule.
    type_order : optional list of type labels
        Override the canonical type ordering (order of the returned
        ``types`` and of combo components/parents).

    Returns
    -------
    PairsResult
    """
    bulk = load_bulk(bulk)
    reactions = dedupe_reactions(list(reactions))
    site_info = site_table(bulk, symprec=site_symprec)
    types, type_by_reaction = build_types(
        bulk,
        reactions,
        site_symprec=site_symprec,
        type_order=type_order,
        site_info=site_info,
    )
    label_index = {t.label: t.index for t in types}

    stats: dict = {"single_defect": {}, "jobs": [], "empty_combos": [], "dropped": 0}

    # ---- step 1: single-defect bases -------------------------------------
    # One base structure per defect type: the configuration whose dopant sits
    # closest to the supercell center, chosen over *all* sites of the type's
    # pristine site class (not just orbit representatives at the strict
    # enumeration tolerance).  This reproduces the base selection of the
    # original SAGAR workflow (verified against its saved base files) and
    # gives the most isolated defect position.
    bases: dict[str, SingleDefectBase] = {}
    for t in types:
        members = [
            i
            for i, s in enumerate(site_info)
            if s.element == t.host
            and s.site_symmetry == t.site_symmetry
            and s.wyckoff == t.site_tag
        ]
        if not members:
            stats["single_defect"][t.label] = {"n_sites": 0}
            continue
        site = min(members, key=lambda i: fractional_distance_to(bulk[i].frac_coords, center))
        d = fractional_distance_to(bulk[site].frac_coords, center)
        structure = bulk.copy()
        structure.replace(site, t.dopant)
        bases[t.label] = SingleDefectBase(
            structure=structure,
            dopant_site=site,
            distance_to_center=d,
            degeneracy=len(members),
            config_index=members.index(site),  # ordinal within the class
        )
        stats["single_defect"][t.label] = {
            "n_sites": len(members),
            "center_distance": round(d, 4),
        }

    # ---- step 2+3: pair combos -------------------------------------------
    complexes: dict[str, PairComplex] = {}
    n_types = len(types)
    for i, ti in enumerate(types):
        base_rec = bases[ti.label]
        seen_reactions: set = set()
        for j in range(i, n_types):
            tj = types[j]
            r = reactions.index((tj.host, tj.dopant))
            if r in seen_reactions:
                continue
            seen_reactions.add(r)
            host, dop = reactions[r]
            # all host sites of this reaction may already be consumed by the
            # base defect (e.g. self-pair with a single host site): skip
            n_host = sum(1 for s in base_rec.structure if s.species_string == host)
            if n_host == 0:
                stats["jobs"].append(
                    {
                        "base": ti.label,
                        "reaction": f"{dop}@{host}",
                        "n_configs": 0,
                        "skipped_no_host": True,
                    }
                )
                continue
            configs = engine.enumerate_substitutions(
                base_rec.structure, host, dop, symprec=enum_symprec
            )
            stats["jobs"].append(
                {"base": ti.label, "reaction": f"{dop}@{host}", "n_configs": len(configs)}
            )
            for cfg in configs:
                t_c = type_by_reaction[r][cfg.site_index]
                if t_c is None:
                    stats["dropped"] += 1
                    continue
                c = label_index[t_c.label]
                if c < i:
                    stats["dropped"] += 1
                    continue
                key = f"{ti.label}+{t_c.label}"
                distance = engine.periodic_distance(
                    cfg.structure, base_rec.dopant_site, cfg.site_index
                )
                cand = PairCandidate(
                    config_index=cfg.config_index,
                    site_index=cfg.site_index,
                    structure=cfg.structure,
                    distance=distance,
                    shell=-1,  # assigned below
                    degeneracy=cfg.degeneracy,
                )
                complexes.setdefault(
                    key,
                    PairComplex(
                        combo=key,
                        type_a=ti,
                        type_b=t_c,
                        base_structure=base_rec.structure,
                        base_dopant_site=base_rec.dopant_site,
                    ),
                ).candidates.append(cand)

    # ---- step 4: shells and representatives ------------------------------
    keep_shells = {"nn": {0}, "nnn": {0, 1}}.get(shells)
    if keep_shells is None:
        raise ValueError(f"shells must be 'nn' or 'nnn', got {shells!r}")

    for key, comp in complexes.items():
        comp.candidates.sort(key=lambda c: (c.distance, c.config_index))
        groups = shell_groups([c.distance for c in comp.candidates], dist_tol)
        for cand, g in zip(comp.candidates, groups, strict=True):
            cand.shell = g
        reps: dict[int, PairCandidate] = {}
        for cand in comp.candidates:
            reps.setdefault(cand.shell, cand)
        comp.shells = [
            PairShell(
                shell=g,
                label=SHELL_LABELS[g],
                candidate=reps[g],
                distance=reps[g].distance,
            )
            for g in sorted(reps)
            if g in keep_shells
        ]
        if not comp.candidates:
            stats["empty_combos"].append(key)

    complexes = {k: v for k, v in complexes.items() if v.candidates}
    stats["n_combos"] = len(complexes)
    return PairsResult(
        bulk=bulk,
        reactions=reactions,
        types=types,
        bases=bases,
        complexes=complexes,
        dist_tol=dist_tol,
        site_type_maps=type_by_reaction,
        stats=stats,
    )


def pair_rows(result: PairsResult) -> list[dict]:
    """Flatten the pair representatives into manifest rows (dicts)."""
    rows: list[dict] = []
    for key in sorted(result.complexes):
        comp = result.complexes[key]
        for shell in comp.shells:
            rows.append(
                {
                    "combo": comp.combo,
                    "shell": shell.label,
                    "dist_AB": round(shell.distance, 4),
                    "config_index": shell.candidate.config_index,
                    "degeneracy": shell.candidate.degeneracy,
                }
            )
    return rows
