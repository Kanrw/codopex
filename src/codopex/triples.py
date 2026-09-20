"""Phase-2 pipeline: placing a third dopant on a defect pair.

Faithful, structure-level port of the ``gen_codope_3def`` notebook logic.

A triple complex is generated from its *canonical parent* pair: a triple
whose sorted defect types are (i, j, c) with i <= j <= c is built by placing
a type-c dopant around the representative structure of the pair (i, j), the
pair made of its two lexicographically first types.  Requesting
``pair="A+B"`` therefore returns exactly the triple complexes whose
canonical parent is that pair -- matching the notebook's per-parent jobs.

For each resulting triple complex, the notebook's three placement criteria
each select one representative structure:

* ``nearA_farB`` -- third dopant in the chosen distance shell of A,
  maximizing d_BC (near A, far from B);
* ``nearB_farA`` -- symmetric;
* ``both_nn`` -- third dopant in the chosen shell of *both* A and B,
  minimizing the mean distance; falls back to the overall mean-distance
  minimum when no configuration satisfies both constraints.

The shell ("nn" = first, "nnn" = second) is controlled by ``shells``.
All inequivalent configurations are kept in ``complexes[combo].candidates``
together with their d_AC / d_BC distances and degeneracies.

:func:`generate_all_triples` runs the same pipeline over every parent pair
(or a requested subset) in one call and reports pairs that lack the
requested representative instead of failing.
"""

from __future__ import annotations

from collections.abc import Sequence

from codopex import engine
from codopex.pipeline import SHELL_LABELS, shell_groups
from codopex.symmetry import equivalent_atoms
from codopex.types import (
    AllTriplesResult,
    PairsResult,
    TripleCandidate,
    TripleComplex,
    TripleSelection,
    TriplesResult,
)


def generate_triples(
    pairs_result: PairsResult,
    pair: str,
    shells: str = "nn",
    pair_shell: str = "nn",
    third: list | None = None,
    enum_symprec: float = 1e-3,
) -> TriplesResult:
    """Generate triple complexes on top of one defect-pair combo.

    Parameters
    ----------
    pairs_result : PairsResult
        Output of :func:`codopex.generate_pairs`.
    pair : str
        Combo label of the pair, e.g. ``"Cu_Pb_C3v+S_O_C3v"``; must be one
        of the combos generated in ``pairs_result``.
    shells : "nn" or "nnn"
        Distance shell used by the placement criteria (first or second
        shell of d_AC / d_BC).
    pair_shell : "nn" or "nnn"
        Which representative of the pair combo is used as the base structure
        ("nnn" requires the pair to have been generated with
        ``generate_pairs(shells="nnn")``).
    third : optional list of (host, dopant)
        Restrict the added dopants to these reactions.  Default: all
        reactions used in the pairs run.
    enum_symprec : float
        spglib symmetry tolerance for configuration enumeration.

    Returns
    -------
    TriplesResult
    """
    if shells not in SHELL_LABELS.values():
        raise ValueError(f"shells must be 'nn' or 'nnn', got {shells!r}")
    if pair_shell not in SHELL_LABELS.values():
        raise ValueError(f"pair_shell must be 'nn' or 'nnn', got {pair_shell!r}")

    types = pairs_result.types
    label_index = {t.label: t.index for t in types}
    reactions = pairs_result.reactions
    site_type_maps = pairs_result.site_type_maps
    dist_tol = pairs_result.dist_tol

    if pair not in pairs_result.complexes:
        raise ValueError(
            f"pair {pair!r} not among generated combos. Available: {sorted(pairs_result.complexes)}"
        )
    comp = pairs_result.complexes[pair]
    rep = next((s for s in comp.shells if s.label == pair_shell), None)
    if rep is None:
        raise ValueError(
            f"pair {pair!r} has no {pair_shell!r} representative "
            f"(generate_pairs(shells=...) must include it). "
            f"Available shells: {[s.label for s in comp.shells]}"
        )

    j = comp.type_b.index
    base_pair = rep.candidate.structure
    a_idx = comp.base_dopant_site  # type-A dopant site (index into base_pair)
    b_idx = rep.candidate.site_index  # type-B dopant site (index into base_pair)

    third_set = {tuple(t) for t in third} if third is not None else None

    # reactions of the owned third types, in canonical type order
    jobs: list[int] = []
    for c in range(j, len(types)):
        r = reactions.index((types[c].host, types[c].dopant))
        if third_set is not None and (types[c].host, types[c].dopant) not in third_set:
            continue
        if r not in jobs:
            jobs.append(r)

    stats: dict = {
        "pair": pair,
        "pair_shell": pair_shell,
        "shells": shells,
        "jobs": [],
        "dropped": 0,
        "empty": [],
    }

    # one spglib search on the pair base, reused by every reaction below
    equiv = equivalent_atoms(base_pair, symprec=enum_symprec)

    complexes: dict[str, TripleComplex] = {}
    for r in jobs:
        host, dop = reactions[r]
        # all host sites may already be consumed by the pair
        if not any(s.species_string == host for s in base_pair):
            stats["jobs"].append(
                {"reaction": f"{dop}@{host}", "n_configs": 0, "skipped_no_host": True}
            )
            continue
        configs = engine.enumerate_substitutions(
            base_pair, host, dop, symprec=enum_symprec, equiv=equiv
        )
        stats["jobs"].append({"reaction": f"{dop}@{host}", "n_configs": len(configs)})
        for cfg in configs:
            t_c = site_type_maps[r][cfg.site_index]
            if t_c is None:
                stats["dropped"] += 1
                continue
            c = label_index[t_c.label]
            if c < j:  # belongs to a triple whose canonical parent is earlier
                stats["dropped"] += 1
                continue
            key = "+".join([comp.type_a.label, comp.type_b.label, t_c.label])
            d_ac = engine.periodic_distance(cfg.structure, a_idx, cfg.site_index)
            d_bc = engine.periodic_distance(cfg.structure, b_idx, cfg.site_index)
            cand = TripleCandidate(
                config_index=cfg.config_index,
                site_index=cfg.site_index,
                structure=cfg.structure,
                d_ac=d_ac,
                d_bc=d_bc,
                degeneracy=cfg.degeneracy,
                type_c=t_c,
            )
            complexes.setdefault(
                key, TripleComplex(combo=key, pair=pair, type_c=t_c)
            ).candidates.append(cand)

    # ---- criteria selections ------------------------------------------------
    tier = {"nn": 0, "nnn": 1}[shells]
    for tc in complexes.values():
        _select_criteria(tc, tier, dist_tol)

    for key, tc in complexes.items():
        if not tc.candidates:
            stats["empty"].append(key)
    complexes = {k: v for k, v in complexes.items() if v.candidates}
    stats["n_triples"] = len(complexes)
    return TriplesResult(
        pairs=pairs_result,
        pair=pair,
        pair_shell=pair_shell,
        complexes=complexes,
        stats=stats,
    )


def generate_all_triples(
    pairs_result: PairsResult,
    pairs: Sequence[str] | None = None,
    shells: str = "nn",
    pair_shell: str = "nn",
    third: list | None = None,
    enum_symprec: float = 1e-3,
) -> AllTriplesResult:
    """Generate triple complexes for many parent pairs in one call.

    This is a convenience wrapper around :func:`generate_triples`: every
    requested parent pair is expanded with the same options.  Pairs that have
    no representative in ``pair_shell`` (for example ``nnn`` was not
    generated for them) are reported in ``AllTriplesResult.skipped`` instead
    of raising, so a batch over ``shells="nnn"`` pairs stays usable.

    Parameters
    ----------
    pairs_result : PairsResult
        Output of :func:`codopex.generate_pairs`.
    pairs : optional sequence of str
        Parent pair combos to expand.  Default: all combos in
        ``pairs_result``, in sorted order.
    shells, pair_shell, third, enum_symprec
        Forwarded to :func:`generate_triples` for every parent.

    Returns
    -------
    AllTriplesResult
        ``results`` maps parent pair label -> :class:`TriplesResult`;
        ``skipped`` lists the parents without the requested shell.
    """
    if shells not in SHELL_LABELS.values():
        raise ValueError(f"shells must be 'nn' or 'nnn', got {shells!r}")
    if pair_shell not in SHELL_LABELS.values():
        raise ValueError(f"pair_shell must be 'nn' or 'nnn', got {pair_shell!r}")

    if pairs is None:
        requested = sorted(pairs_result.complexes)
    else:
        requested = list(pairs)
        missing = [p for p in requested if p not in pairs_result.complexes]
        if missing:
            raise ValueError(
                f"pairs not among generated combos: {missing}. "
                f"Available: {sorted(pairs_result.complexes)}"
            )

    results: dict[str, TriplesResult] = {}
    skipped: list[str] = []
    for parent in requested:
        comp = pairs_result.complexes[parent]
        if not any(s.label == pair_shell for s in comp.shells):
            skipped.append(parent)
            continue
        results[parent] = generate_triples(
            pairs_result,
            parent,
            shells=shells,
            pair_shell=pair_shell,
            third=third,
            enum_symprec=enum_symprec,
        )
    return AllTriplesResult(results=results, skipped=skipped)


def _select_criteria(tc: TripleComplex, tier: int, dist_tol: float) -> None:
    """Apply the three placement criteria to a triple complex.

    ``tier`` selects the distance shell used by the criteria (0 = first
    shell "nn", 1 = second shell "nnn").
    """
    cands = tc.candidates
    # deterministic: candidates stay in enumeration order; shells are
    # computed over the distances sorted ascending (stable)
    order_ac = sorted(range(len(cands)), key=lambda k: (cands[k].d_ac, k))
    order_bc = sorted(range(len(cands)), key=lambda k: (cands[k].d_bc, k))
    ac_group = shell_groups([cands[k].d_ac for k in order_ac], dist_tol)
    bc_group = shell_groups([cands[k].d_bc for k in order_bc], dist_tol)
    ac_shell = {k: g for g, k in zip(ac_group, order_ac, strict=True)}
    bc_shell = {k: g for g, k in zip(bc_group, order_bc, strict=True)}

    avg = lambda c: (c.d_ac + c.d_bc) / 2  # noqa: E731

    def pick_max(seq, key):
        best, best_key = None, -1.0
        for c in seq:  # first maximum wins on ties (stable)
            k = key(c)
            if k > best_key:
                best, best_key = c, k
        return best

    near_a = [c for k, c in enumerate(cands) if ac_shell[k] == tier]
    near_b = [c for k, c in enumerate(cands) if bc_shell[k] == tier]
    both = [c for k, c in enumerate(cands) if ac_shell[k] == tier and bc_shell[k] == tier]

    tc.selections = [
        TripleSelection("nearA_farB", pick_max(near_a, lambda c: c.d_bc)),
        TripleSelection("nearB_farA", pick_max(near_b, lambda c: c.d_ac)),
        # notebook fallback: when nothing is NN to both, take the overall
        # minimum mean distance
        TripleSelection("both_nn", min(both, key=avg) if both else min(cands, key=avg)),
    ]


def triple_rows(result: TriplesResult) -> list[dict]:
    """Flatten triple representatives into manifest rows (dicts)."""
    rows: list[dict] = []
    for tc in (result.complexes[k] for k in sorted(result.complexes)):
        for sel in tc.selections:
            cand = sel.candidate
            if cand is None:
                continue
            rows.append(
                {
                    "combo": tc.combo,
                    "criterion": sel.criterion,
                    "d_AC": round(cand.d_ac, 4),
                    "d_BC": round(cand.d_bc, 4),
                    "avg_dist": round((cand.d_ac + cand.d_bc) / 2, 4),
                    "config_index": cand.config_index,
                    "degeneracy": cand.degeneracy,
                }
            )
    return rows
