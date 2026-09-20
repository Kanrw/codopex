"""Export and persistence helpers.

The generation pipeline itself is purely in-memory; these helpers turn a
``PairsResult`` / ``TriplesResult`` into the directory layout used for VASP
pre-processing (mirroring the original notebooks):

.. code-block:: text

    out_root/
      single_defect/{type}/POSCAR            # base structure of each type
      {combo}/{shell}/POSCAR{idx}            # representatives (idx = config)
      {combo}/{shell}/POSCAR                 # unindexed copy of the same
      work_{combo}/POSCAR{idx}               # (write_all=True) all candidates

where ``shell`` is "nn" or "nnn" for pairs and the criterion name
("nearA_farB", ...) for triples.

POSCAR files are written with species grouped by atomic number (VASP
convention of the original DefectMaker output).  Structures returned by the
pipeline itself keep the pristine site order.

Results can also be persisted with :func:`save` / :func:`load` (JSON, gzip
when the path ends in ``.gz``), so phase 1 and phase 2 can run in separate
sessions or machines::

    cp.io.save(pairs, "runs/pairs.json.gz")
    pairs = cp.io.load("runs/pairs.json.gz")   # later, elsewhere
"""

from __future__ import annotations

import csv
import gzip
import json
import os
from collections.abc import Sequence

from pymatgen.core import Element, Structure

from codopex._version import __version__
from codopex.types import (
    AllTriplesResult,
    DefectType,
    PairCandidate,
    PairComplex,
    PairShell,
    PairsResult,
    SingleDefectBase,
    TripleCandidate,
    TripleComplex,
    TripleSelection,
    TriplesResult,
)

Result = PairsResult | TriplesResult | AllTriplesResult

RESULT_FORMAT = "codopex-result"
RESULT_FORMAT_VERSION = 1


def write_structure(structure: Structure, path: str) -> str:
    """Write ``structure`` to ``path`` as POSCAR, species sorted by Z.

    Returns the path.
    """
    order = sorted(range(len(structure)), key=lambda i: (Element(structure[i].species_string).Z, i))
    ordered = Structure(
        structure.lattice,
        [structure[i].species_string for i in order],
        [structure[i].frac_coords for i in order],
        validate_proximity=False,
        to_unit_cell=False,
    )
    ordered.to(filename=path, fmt="poscar")
    return path


def export_pairs(
    result: PairsResult,
    out_root: str,
    write_all: bool = False,
) -> dict[str, list[str]]:
    """Write pair structures into ``out_root`` (see module docstring).

    Returns a mapping ``combo -> written representative paths``.
    """
    written: dict[str, list[str]] = {}
    os.makedirs(out_root, exist_ok=True)

    for label, base in result.bases.items():
        folder = os.path.join(out_root, "single_defect", label)
        os.makedirs(folder, exist_ok=True)
        write_structure(base.structure, os.path.join(folder, "POSCAR"))

    for combo, comp in result.complexes.items():
        paths: list[str] = []
        for shell in comp.shells:
            folder = os.path.join(out_root, combo, shell.label)
            os.makedirs(folder, exist_ok=True)
            cand = shell.candidate
            p1 = write_structure(cand.structure, os.path.join(folder, f"POSCAR{cand.config_index}"))
            p2 = write_structure(cand.structure, os.path.join(folder, "POSCAR"))
            paths.extend([p1, p2])
        written[combo] = paths
        if write_all:
            work = os.path.join(out_root, f"work_{combo}")
            os.makedirs(work, exist_ok=True)
            for cand in comp.candidates:
                write_structure(cand.structure, os.path.join(work, f"POSCAR{cand.config_index}"))
    return written


def export_triples(
    result: TriplesResult,
    out_root: str,
    write_all: bool = False,
) -> dict[str, list[str]]:
    """Write triple structures into ``out_root`` (see module docstring)."""
    written: dict[str, list[str]] = {}
    os.makedirs(out_root, exist_ok=True)

    for combo, tc in result.complexes.items():
        paths: list[str] = []
        for sel in tc.selections:
            cand = sel.candidate
            if cand is None:
                continue
            folder = os.path.join(out_root, combo, sel.criterion)
            os.makedirs(folder, exist_ok=True)
            p1 = write_structure(cand.structure, os.path.join(folder, f"POSCAR{cand.config_index}"))
            p2 = write_structure(cand.structure, os.path.join(folder, "POSCAR"))
            paths.extend([p1, p2])
        written[combo] = paths
        if write_all:
            work = os.path.join(out_root, f"work_{combo}")
            os.makedirs(work, exist_ok=True)
            for cand in tc.candidates:
                write_structure(cand.structure, os.path.join(work, f"POSCAR{cand.config_index}"))
    return written


def write_rows_csv(rows: Sequence[dict], path: str) -> str:
    """Write manifest rows (list of dicts) as UTF-8-BOM CSV (Excel friendly)."""
    if not rows:
        raise ValueError("no rows to write")
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# persistence (JSON, optional gzip)
# ---------------------------------------------------------------------------


def _open_text(path: str, mode: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def _index_of(items: Sequence, obj) -> int:
    for i, item in enumerate(items):
        if item is obj:
            return i
    raise ValueError("object is not part of the sequence (internal error)")


def _encode_structure(structure: Structure) -> dict:
    return structure.as_dict()


def _decode_structure(data: dict) -> Structure:
    return Structure.from_dict(data)


def _encode_type(t: DefectType) -> dict:
    return {
        "host": t.host,
        "dopant": t.dopant,
        "site_symmetry": t.site_symmetry,
        "site_tag": t.site_tag,
        "index": t.index,
        "label": t.label,
    }


def _decode_type(d: dict) -> DefectType:
    return DefectType(
        host=d["host"],
        dopant=d["dopant"],
        site_symmetry=d["site_symmetry"],
        site_tag=d["site_tag"],
        index=d["index"],
        label=d["label"],
    )


def _encode_pairs(result: PairsResult) -> dict:
    return {
        "bulk": _encode_structure(result.bulk),
        "reactions": [list(r) for r in result.reactions],
        "types": [_encode_type(t) for t in result.types],
        "bases": {
            label: {
                "structure": _encode_structure(base.structure),
                "dopant_site": base.dopant_site,
                "distance_to_center": base.distance_to_center,
                "degeneracy": base.degeneracy,
            }
            for label, base in result.bases.items()
        },
        "complexes": {
            key: {
                "combo": comp.combo,
                "type_a": comp.type_a.label,
                "type_b": comp.type_b.label,
                "base_structure": _encode_structure(comp.base_structure),
                "base_dopant_site": comp.base_dopant_site,
                "candidates": [
                    {
                        "config_index": cand.config_index,
                        "site_index": cand.site_index,
                        "structure": _encode_structure(cand.structure),
                        "distance": cand.distance,
                        "shell": cand.shell,
                        "degeneracy": cand.degeneracy,
                    }
                    for cand in comp.candidates
                ],
                "shells": [
                    {
                        "label": shell.label,
                        "candidate": _index_of(comp.candidates, shell.candidate),
                        "distance": shell.distance,
                    }
                    for shell in comp.shells
                ],
            }
            for key, comp in result.complexes.items()
        },
        "dist_tol": result.dist_tol,
        "site_type_maps": [
            {str(site): (t.label if t is not None else None) for site, t in mapping.items()}
            for mapping in result.site_type_maps
        ],
        "stats": result.stats,
    }


def _decode_pairs(data: dict) -> PairsResult:
    types = [_decode_type(t) for t in data["types"]]
    by_label = {t.label: t for t in types}

    bases = {
        label: SingleDefectBase(
            structure=_decode_structure(b["structure"]),
            dopant_site=b["dopant_site"],
            distance_to_center=b["distance_to_center"],
            degeneracy=b["degeneracy"],
        )
        for label, b in data["bases"].items()
    }

    complexes: dict[str, PairComplex] = {}
    for key, c in data["complexes"].items():
        candidates = [
            PairCandidate(
                config_index=cd["config_index"],
                site_index=cd["site_index"],
                structure=_decode_structure(cd["structure"]),
                distance=cd["distance"],
                shell=cd["shell"],
                degeneracy=cd["degeneracy"],
            )
            for cd in c["candidates"]
        ]
        shells = [
            PairShell(
                label=s["label"],
                candidate=candidates[s["candidate"]],
                distance=s["distance"],
            )
            for s in c["shells"]
        ]
        complexes[key] = PairComplex(
            combo=c["combo"],
            type_a=by_label[c["type_a"]],
            type_b=by_label[c["type_b"]],
            base_structure=_decode_structure(c["base_structure"]),
            base_dopant_site=c["base_dopant_site"],
            candidates=candidates,
            shells=shells,
        )

    site_type_maps = [
        {
            int(site): (by_label[label] if label is not None else None)
            for site, label in mapping.items()
        }
        for mapping in data["site_type_maps"]
    ]

    return PairsResult(
        bulk=_decode_structure(data["bulk"]),
        reactions=[tuple(r) for r in data["reactions"]],
        types=types,
        bases=bases,
        complexes=complexes,
        dist_tol=data["dist_tol"],
        site_type_maps=site_type_maps,
        stats=data["stats"],
    )


def _encode_triples(result: TriplesResult, include_pairs: bool = True) -> dict:
    data = {
        "pair": result.pair,
        "pair_shell": result.pair_shell,
        "complexes": {
            key: {
                "combo": tc.combo,
                "pair": tc.pair,
                "type_c": tc.type_c.label,
                "candidates": [
                    {
                        "config_index": cand.config_index,
                        "site_index": cand.site_index,
                        "structure": _encode_structure(cand.structure),
                        "d_ac": cand.d_ac,
                        "d_bc": cand.d_bc,
                        "degeneracy": cand.degeneracy,
                    }
                    for cand in tc.candidates
                ],
                "selections": [
                    {
                        "criterion": sel.criterion,
                        "candidate": (
                            _index_of(tc.candidates, sel.candidate)
                            if sel.candidate is not None
                            else None
                        ),
                    }
                    for sel in tc.selections
                ],
            }
            for key, tc in result.complexes.items()
        },
        "stats": result.stats,
    }
    if include_pairs:
        data["pairs"] = _encode_pairs(result.pairs)
    return data


def _decode_triples(data: dict, pairs: PairsResult | None = None) -> TriplesResult:
    if pairs is None:
        pairs = _decode_pairs(data["pairs"])
    by_label = {t.label: t for t in pairs.types}

    complexes: dict[str, TripleComplex] = {}
    for key, tc in data["complexes"].items():
        candidates = [
            TripleCandidate(
                config_index=cd["config_index"],
                site_index=cd["site_index"],
                structure=_decode_structure(cd["structure"]),
                d_ac=cd["d_ac"],
                d_bc=cd["d_bc"],
                degeneracy=cd["degeneracy"],
            )
            for cd in tc["candidates"]
        ]
        selections = [
            TripleSelection(
                criterion=sel["criterion"],
                candidate=(candidates[sel["candidate"]] if sel["candidate"] is not None else None),
            )
            for sel in tc["selections"]
        ]
        complexes[key] = TripleComplex(
            combo=tc["combo"],
            pair=tc["pair"],
            type_c=by_label[tc["type_c"]],
            candidates=candidates,
            selections=selections,
        )

    return TriplesResult(
        pairs=pairs,
        pair=data["pair"],
        pair_shell=data["pair_shell"],
        complexes=complexes,
        stats=data["stats"],
    )


def _encode_batch(result: AllTriplesResult) -> dict:
    parents = list(result.results)
    data: dict = {
        "results": {
            parent: _encode_triples(result.results[parent], include_pairs=False)
            for parent in parents
        },
        "skipped": list(result.skipped),
    }
    if parents:  # the shared phase-1 result is stored exactly once
        data["pairs"] = _encode_pairs(result.results[parents[0]].pairs)
    return data


def _decode_batch(data: dict) -> AllTriplesResult:
    pairs = _decode_pairs(data["pairs"]) if "pairs" in data else None
    results = {
        parent: _decode_triples(payload, pairs=pairs) for parent, payload in data["results"].items()
    }
    return AllTriplesResult(results=results, skipped=list(data["skipped"]))


def save(result: Result, path: str) -> str:
    """Save a result (:class:`PairsResult`, :class:`TriplesResult` or
    :class:`AllTriplesResult`) as JSON.

    The file is gzip-compressed when ``path`` ends in ``.gz``.  Structures are
    stored with pymatgen's own serialization, so a loaded result is fully
    equivalent to the in-memory one (including all candidates, degeneracies
    and stats).  Returns the path.
    """
    if isinstance(result, PairsResult):
        kind, payload = "pairs", _encode_pairs(result)
    elif isinstance(result, TriplesResult):
        kind, payload = "triples", _encode_triples(result)
    elif isinstance(result, AllTriplesResult):
        kind, payload = "triples-batch", _encode_batch(result)
    else:
        raise TypeError(
            "save() expects a PairsResult, TriplesResult or AllTriplesResult, "
            f"got {type(result).__name__}"
        )
    document = {
        "format": RESULT_FORMAT,
        "format_version": RESULT_FORMAT_VERSION,
        "codopex_version": __version__,
        "kind": kind,
        "data": payload,
    }
    with _open_text(path, "w") as f:
        json.dump(document, f)
    return path


def load(path: str) -> Result:
    """Load a result written by :func:`save`.

    Returns a :class:`PairsResult`, :class:`TriplesResult` or
    :class:`AllTriplesResult` depending on what was saved.  Raises
    ``ValueError`` for files that are not codopex result documents or use an
    unsupported format version.
    """
    with _open_text(path, "r") as f:
        document = json.load(f)
    if not isinstance(document, dict) or document.get("format") != RESULT_FORMAT:
        raise ValueError(f"{path} is not a codopex result file")
    version = document.get("format_version")
    if version != RESULT_FORMAT_VERSION:
        raise ValueError(
            f"unsupported codopex result format version {version!r} "
            f"(this version reads {RESULT_FORMAT_VERSION})"
        )
    kind = document.get("kind")
    if kind == "pairs":
        return _decode_pairs(document["data"])
    if kind == "triples":
        return _decode_triples(document["data"])
    if kind == "triples-batch":
        return _decode_batch(document["data"])
    raise ValueError(f"unknown codopex result kind {kind!r} in {path}")
