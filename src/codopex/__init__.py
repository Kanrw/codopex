"""codopex: symmetry-inequivalent co-doping defect complexes.

Generate nearest-neighbour (and next-nearest) defect-pair complexes in a
host supercell, and expand a pair by a third dopant into triple complexes
placed near A, near B, or near both -- all configurations enumerated under
the symmetry of the host (pure pymatgen/spglib, no external enumerator).

Typical usage::

    import codopex as cp

    # Phase 1: all nearest pairs among Cu@Pb, S@O, S@P
    pairs = cp.generate_pairs("CONTCAR", [("Pb", "Cu"), ("O", "S"), ("P", "S")])

    # Phase 2: expand one pair by an additional S@O dopant
    triples = cp.generate_triples(pairs, "Cu_Pb_C3v+S_O_C3v", third=[("O", "S")])

Structures are returned in-memory; use :mod:`codopex.io` to export POSCAR
trees and CSV manifests for VASP runs.
"""

from codopex import io
from codopex.classify import DefectRecord, classify_defects
from codopex.pairs import generate_pairs, pair_rows
from codopex.triples import generate_triples, triple_rows
from codopex.types import (
    DefectType,
    PairCandidate,
    PairComplex,
    PairsResult,
    Reaction,
    SingleDefectBase,
    TripleCandidate,
    TripleComplex,
    TripleSelection,
    TriplesResult,
)

__version__ = "0.1.0"

__all__ = [
    "DefectRecord",
    "DefectType",
    "PairCandidate",
    "PairComplex",
    "PairsResult",
    "Reaction",
    "SingleDefectBase",
    "TripleCandidate",
    "TripleComplex",
    "TripleSelection",
    "TriplesResult",
    "classify_defects",
    "generate_pairs",
    "generate_triples",
    "pair_rows",
    "triple_rows",
]
