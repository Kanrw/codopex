"""Command-line interface for codopex.

The CLI mirrors the Python API:

.. code-block:: console

    # phase 1: all nearest pairs among Cu@Pb, S@O and S@P
    codopex pairs CONTCAR -r Cu@Pb -r S@O -r S@P --shells nnn \
        --out runs/pairs --save runs/pairs/pairs.json.gz

    # phase 2: expand every pair by an additional S@O dopant
    codopex triples --load-pairs runs/pairs/pairs.json.gz -r O@S \
        --out runs/triples --save runs/triples/triples.json.gz

    # label the defects of an externally produced structure
    codopex classify bulk.vasp defect.vasp

Reactions are written ``DOPANT@HOST`` (e.g. ``Cu@Pb`` = Cu on Pb sites).
"""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Sequence

from codopex import io
from codopex._version import __version__
from codopex.classify import classify_defects
from codopex.pairs import generate_pairs, pair_rows
from codopex.pipeline import load_bulk
from codopex.triples import generate_all_triples, triple_rows
from codopex.types import PairsResult, Reaction


def _reaction(text: str) -> Reaction:
    """Parse a ``DOPANT@HOST`` CLI token into a (host, dopant) reaction."""
    if "@" not in text:
        raise argparse.ArgumentTypeError(f"expected DOPANT@HOST (e.g. Cu@Pb), got {text!r}")
    dopant, host = (part.strip() for part in text.split("@", 1))
    if not dopant or not host:
        raise argparse.ArgumentTypeError(f"expected DOPANT@HOST (e.g. Cu@Pb), got {text!r}")
    return (host, dopant)


def _type_order(text: str | None) -> list[str] | None:
    if text is None:
        return None
    labels = [label.strip() for label in text.split(",") if label.strip()]
    if not labels:
        raise ValueError("--type-order must list at least one type label")
    return labels


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="codopex",
        description="Symmetry-inequivalent co-doping defect-complex generation.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    pairs = sub.add_parser("pairs", help="generate symmetry-inequivalent defect pairs")
    pairs.add_argument("structure", help="pristine supercell (CONTCAR/POSCAR/CIF/...)")
    pairs.add_argument(
        "-r",
        "--reaction",
        dest="reactions",
        action="append",
        required=True,
        type=_reaction,
        metavar="DOPANT@HOST",
        help="substitution reaction (repeatable), e.g. Cu@Pb",
    )
    pairs.add_argument("--shells", choices=("nn", "nnn"), default="nn")
    pairs.add_argument("--dist-tol", type=float, default=0.05)
    pairs.add_argument(
        "--symprec",
        type=float,
        default=1e-3,
        help="spglib tolerance (A) for configuration enumeration",
    )
    pairs.add_argument("--site-symprec", type=float, default=0.01)
    pairs.add_argument(
        "--type-order",
        metavar="L1,L2,...",
        help="comma-separated type labels overriding the canonical order",
    )
    pairs.add_argument("-o", "--out", help="write a POSCAR tree + manifest.csv here")
    pairs.add_argument(
        "--write-all",
        action="store_true",
        help="also write every candidate under work_{combo}/",
    )
    pairs.add_argument("--save", metavar="FILE", help="save the full result as JSON")
    pairs.add_argument("--csv", metavar="FILE", help="write the manifest to this CSV")
    pairs.set_defaults(handler=_cmd_pairs)

    triples = sub.add_parser("triples", help="expand parent pair(s) by a third dopant")
    triples.add_argument(
        "structure",
        nargs="?",
        help="pristine supercell (not needed with --load-pairs)",
    )
    triples.add_argument(
        "-r",
        "--reaction",
        dest="reactions",
        action="append",
        type=_reaction,
        metavar="DOPANT@HOST",
        help="substitution reaction of phase 1 (repeatable)",
    )
    triples.add_argument(
        "--load-pairs",
        metavar="FILE",
        help="reuse a phase-1 result saved with 'pairs --save'",
    )
    triples.add_argument(
        "--pair",
        dest="pairs",
        action="append",
        metavar="COMBO",
        help="parent pair combo (repeatable); default: every generated pair",
    )
    triples.add_argument("--shells", choices=("nn", "nnn"), default="nn")
    triples.add_argument("--pair-shell", choices=("nn", "nnn"), default="nn")
    triples.add_argument(
        "--third",
        dest="third",
        action="append",
        type=_reaction,
        metavar="DOPANT@HOST",
        help="restrict the third dopant to this reaction (repeatable)",
    )
    triples.add_argument("--dist-tol", type=float, default=0.05)
    triples.add_argument("--symprec", type=float, default=1e-3)
    triples.add_argument("--site-symprec", type=float, default=0.01)
    triples.add_argument(
        "--type-order",
        metavar="L1,L2,...",
        help="comma-separated type labels overriding the canonical order",
    )
    triples.add_argument("-o", "--out", help="write a POSCAR tree + manifest.csv here")
    triples.add_argument(
        "--write-all",
        action="store_true",
        help="also write every candidate under work_{combo}/",
    )
    triples.add_argument("--save", metavar="FILE", help="save the full result as JSON")
    triples.add_argument("--csv", metavar="FILE", help="write the manifest to this CSV")
    triples.set_defaults(handler=_cmd_triples)

    classify = sub.add_parser(
        "classify", help="label the defects of a structure against the pristine bulk"
    )
    classify.add_argument("bulk", help="pristine structure")
    classify.add_argument("defect", help="defective structure (same lattice)")
    classify.add_argument("--symprec", type=float, default=0.01)
    classify.add_argument("--stol", type=float, default=0.5)
    classify.set_defaults(handler=_cmd_classify)

    return parser


def _print_pairs_summary(pairs: PairsResult) -> None:
    print(f"defect types ({len(pairs.types)}):")
    for t in pairs.types:
        print(f"  {t.label:<26} {t.dopant}@{t.host:<4} {t.site_symmetry}")
    print(f"pair combos ({len(pairs.complexes)}):")
    for combo in sorted(pairs.complexes):
        shells = "  ".join(f"{s.label} {s.distance:.4f} A" for s in pairs.complexes[combo].shells)
        print(f"  {combo:<40} {shells}")


def _cmd_pairs(args: argparse.Namespace) -> int:
    pairs = generate_pairs(
        args.structure,
        args.reactions,
        shells=args.shells,
        dist_tol=args.dist_tol,
        enum_symprec=args.symprec,
        site_symprec=args.site_symprec,
        type_order=_type_order(args.type_order),
    )
    _print_pairs_summary(pairs)
    rows = pair_rows(pairs)
    print(f"representatives: {len(rows)}")

    if args.out:
        written = io.export_pairs(pairs, args.out, write_all=args.write_all)
        n_files = sum(len(paths) for paths in written.values())
        print(f"wrote {n_files} POSCAR files under {args.out}")
        if rows:
            manifest = os.path.join(args.out, "manifest.csv")
            io.write_rows_csv(rows, manifest)
            print(f"manifest: {manifest}")
    if args.csv:
        io.write_rows_csv(rows, args.csv)
        print(f"manifest: {args.csv}")
    if args.save:
        io.save(pairs, args.save)
        print(f"saved result: {args.save}")
    return 0


def _cmd_triples(args: argparse.Namespace) -> int:
    if args.load_pairs:
        if args.structure is not None or args.reactions:
            raise ValueError("pass either STRUCTURE with -r/--reaction, or --load-pairs, not both")
        pairs = io.load(args.load_pairs)
        if not isinstance(pairs, PairsResult):
            raise ValueError(f"{args.load_pairs} does not contain a PairsResult")
        print(f"loaded pair result: {args.load_pairs}")
    else:
        if args.structure is None or not args.reactions:
            raise ValueError(
                "STRUCTURE and at least one -r/--reaction are required unless --load-pairs is used"
            )
        pairs = generate_pairs(
            args.structure,
            args.reactions,
            shells=args.pair_shell,  # the pair base must contain this shell
            dist_tol=args.dist_tol,
            enum_symprec=args.symprec,
            site_symprec=args.site_symprec,
            type_order=_type_order(args.type_order),
        )

    batch = generate_all_triples(
        pairs,
        pairs=args.pairs,
        shells=args.shells,
        pair_shell=args.pair_shell,
        third=args.third,
        enum_symprec=args.symprec,
    )
    for parent in batch.skipped:
        print(f"skipped {parent}: no {args.pair_shell!r} representative", file=sys.stderr)
    if len(batch) == 0:
        print("codopex: error: no parent pair produced a triple complex", file=sys.stderr)
        return 1

    rows: list[dict] = []
    for parent in batch:
        print(f"  {parent:<40} {len(batch[parent].complexes)} triple combos")
        rows.extend(triple_rows(batch[parent]))
    print(f"parents: {len(batch)}; triple representatives: {len(rows)}")

    if args.out:
        written: dict[str, list[str]] = {}
        for parent in batch:
            written.update(io.export_triples(batch[parent], args.out, write_all=args.write_all))
        n_files = sum(len(paths) for paths in written.values())
        print(f"wrote {n_files} POSCAR files under {args.out}")
        if rows:
            manifest = os.path.join(args.out, "manifest.csv")
            io.write_rows_csv(rows, manifest)
            print(f"manifest: {manifest}")
    if args.csv:
        io.write_rows_csv(rows, args.csv)
        print(f"manifest: {args.csv}")
    if args.save:
        io.save(batch, args.save)
        print(f"saved result: {args.save}")
    return 0


def _cmd_classify(args: argparse.Namespace) -> int:
    bulk = load_bulk(args.bulk)
    defect = load_bulk(args.defect)
    records = classify_defects(bulk, defect, symprec=args.symprec, stol=args.stol)
    if not records:
        print("no defects found")
        return 0
    for record in records:
        print(f"{record.defect_type:<14} {record.defect_name}")
    print(f"{len(records)} defect(s)")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point of the ``codopex`` console script."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except (ValueError, FileNotFoundError, OSError) as exc:
        print(f"codopex: error: {exc}", file=sys.stderr)
        return 1
