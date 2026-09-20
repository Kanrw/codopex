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

    # ... or expand every pair in one call
    batch = cp.generate_all_triples(pairs)

Structures are returned in-memory; use :mod:`codopex.io` to export POSCAR
trees and CSV manifests for VASP runs, or to save/load results as JSON
(``codopex.io.save`` / ``codopex.io.load``).  The same workflows are
available from the command line (``codopex --help``).
"""

from codopex import io as io
from codopex._version import __version__
from codopex.classify import DefectRecord, classify_defects
from codopex.pairs import generate_pairs, pair_rows
from codopex.triples import generate_all_triples, generate_triples, triple_rows
from codopex.types import (
    AllTriplesResult,
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

__all__ = [
    "AllTriplesResult",
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
    "__version__",
    "classify_defects",
    "generate_all_triples",
    "generate_pairs",
    "generate_triples",
    "pair_rows",
    "triple_rows",
]
