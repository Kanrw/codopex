"""Core data containers for codopex.

All structures handled by the pipeline are ``pymatgen.core.Structure``
objects whose site order is identical to that of the pristine input bulk:
substitutions never move atoms and never reorder sites.  This invariant lets
every atom keep a stable "pristine site index", which is how defects are
labelled and tracked through the pair / triple stages.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

from pymatgen.core import Structure

Reaction = tuple[str, str]  # (host_element, dopant_element)


@dataclass(frozen=True)
class DefectType:
    """A substitutional defect class: dopant on one symmetry class of host sites.

    The class is defined by the pristine (defect-free) symmetry of the host
    site, so that e.g. ``Cu`` on the two inequivalent ``Pb`` positions of
    Pb3(PO4)2 produces two types with different site symmetry.
    """

    host: str
    dopant: str
    site_symmetry: str  # Schoenflies label, e.g. "C3v"
    site_tag: str  # wyckoff-style tag used to disambiguate same-symmetry classes
    index: int  # canonical order index (see PairsResult.types)
    label: str = field(compare=False, hash=False)
    """Unique human-readable label.

    ``f"{dopant}_{host}_{site_symmetry}"``, extended by ``_{site_tag}`` only
    when the plain name would collide with another type (two inequivalent
    host classes sharing the same site-symmetry label).
    """


@dataclass
class SubstitutionConfig:
    """One symmetry-inequivalent single-substitution configuration.

    ``structure`` equals the enumeration base with the dopant placed on site
    ``site_index``; all other sites (positions and species) are untouched.
    """

    site_index: int  # index of the substituted host site in the base
    structure: Structure  # base with dopant on site_index
    degeneracy: int  # number of equivalent host sites in the same orbit
    config_index: int  # 0-based order within the enumeration task


@dataclass
class SingleDefectBase:
    """Base structure for one defect type (its own dopant, closest to center)."""

    structure: Structure  # bulk with one dopant of this type
    dopant_site: int  # index of the dopant atom in ``structure``
    distance_to_center: float  # periodic fractional distance to (0.5,0.5,0.5)
    degeneracy: int  # orbit size of the chosen host site in the bulk
    config_index: int


@dataclass
class PairCandidate:
    """A pair configuration: base type-i dopant fixed + second dopant at site_index."""

    config_index: int
    site_index: int  # host site (in the base structure) carrying the 2nd dopant
    structure: Structure
    distance: float  # minimum periodic distance between the two dopants
    shell: int  # 0-based distance-shell group within the combo
    degeneracy: int


@dataclass
class PairShell:
    """One selected representative (shell "nn", "nnn", ...)."""

    shell: int  # 0-based shell group
    label: str  # "nn", "nnn", ... (see shell_label())
    candidate: PairCandidate
    distance: float


@dataclass
class PairComplex:
    """All generated configurations of one defect-pair combo "A+B"."""

    combo: str  # canonical "A+B" label
    type_a: DefectType
    type_b: DefectType
    base_structure: Structure  # single-defect base of type A used for generation
    base_dopant_site: int  # index of the type-A dopant in base_structure
    candidates: list[PairCandidate] = field(default_factory=list)
    shells: list[PairShell] = field(default_factory=list)  # representatives, in order

    def shell_structure(self, label: str) -> Structure | None:
        for s in self.shells:
            if s.label == label:
                return s.candidate.structure
        return None


@dataclass
class TripleCandidate:
    """A triple configuration: fixed pair A-B + third dopant at site_index."""

    config_index: int
    site_index: int  # host site (in the pair structure) carrying the 3rd dopant
    structure: Structure
    d_ac: float
    d_bc: float
    degeneracy: int
    type_c: DefectType  # pristine type of the third dopant site


@dataclass
class TripleSelection:
    """One selected representative per criterion."""

    criterion: str  # "nearA_farB" | "nearB_farA" | "both_nn"
    candidate: TripleCandidate | None = None


@dataclass
class TripleComplex:
    """Candidates and selected representatives of one triple combo "A+B+C"."""

    combo: str
    pair: str  # the requested pair "A+B"
    type_c: DefectType
    candidates: list[TripleCandidate] = field(default_factory=list)
    selections: list[TripleSelection] = field(default_factory=list)

    def selection(self, criterion: str) -> TripleCandidate | None:
        for s in self.selections:
            if s.criterion == criterion:
                return s.candidate
        return None


# ---------------------------------------------------------------------------
# registry containers
# ---------------------------------------------------------------------------


@dataclass
class PairsResult:
    """Result of :func:`codopex.generate_pairs`.

    Holds the discovered defect types (canonical order), the single-defect
    base structure of every type, and the per-combo pair complexes.
    """

    bulk: Structure
    reactions: list[Reaction]
    types: list[DefectType]  # canonical order
    bases: dict[str, SingleDefectBase]  # type label -> base structure + metadata
    complexes: dict[str, PairComplex]  # combo label -> complex
    dist_tol: float = 0.05
    site_type_maps: list[dict[int, DefectType | None]] = field(default_factory=list)
    """Per reaction (index in ``reactions``): pristine site index -> type."""
    stats: dict = field(default_factory=dict)

    @property
    def type_labels(self) -> list[str]:
        return [t.label for t in self.types]

    def complex_for(self, combo: str) -> PairComplex | None:
        return self.complexes.get(combo)


@dataclass
class TriplesResult:
    """Result of :func:`codopex.generate_triples` for one requested pair."""

    pairs: PairsResult
    pair: str  # requested pair combo
    pair_shell: str  # tier ("nn"/"nnn") of the pair structure used as base
    complexes: dict[str, TripleComplex]  # triple combo label -> complex
    stats: dict = field(default_factory=dict)


@dataclass
class AllTriplesResult:
    """Result of :func:`codopex.generate_all_triples` (batch over parent pairs).

    Iterating the object yields the parent pair labels; ``batch[pair]`` gives
    the :class:`TriplesResult` of that parent.
    """

    results: dict[str, TriplesResult]  # parent pair combo -> result
    skipped: list[str] = field(default_factory=list)
    """Parent pairs without a representative in the requested pair shell."""

    @property
    def pair_labels(self) -> list[str]:
        """Parent pair combos that produced a result, in request order."""
        return list(self.results)

    def __getitem__(self, pair: str) -> TriplesResult:
        return self.results[pair]

    def __iter__(self) -> Iterator[str]:
        return iter(self.results)

    def __len__(self) -> int:
        return len(self.results)
