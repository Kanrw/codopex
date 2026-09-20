"""Small structures shared by the test modules."""

from pymatgen.core import Lattice, Structure


def cubic_2x2x2():
    """Simple cubic lattice (a=4 A), 8 host atoms, one symmetry class."""
    return Structure(
        Lattice.cubic(4.0),
        ["Na"] * 8,
        [[x / 2, y / 2, z / 2] for x in (0, 1) for y in (0, 1) for z in (0, 1)],
    )


def rocksalt_2x2x2():
    """Rocksalt supercell: 64 atoms, Na and Cl each one symmetry class."""
    base = Structure(Lattice.cubic(4.0), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    return base * [2, 2, 2]


def rocksalt_1x1x1():
    return Structure(Lattice.cubic(4.0), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
