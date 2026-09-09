"""Export helpers: write structures and manifests to disk.

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
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, Optional, Sequence

from pymatgen.core import Element, Structure

from codopex.types import PairsResult, TriplesResult


def write_structure(structure: Structure, path: str) -> str:
    """Write ``structure`` to ``path`` as POSCAR, species sorted by Z.

    Returns the path.
    """
    ordered = structure.copy()
    order = sorted(
        range(len(structure)), key=lambda i: (Element(structure[i].species_string).Z, i)
    )
    ordered = Structure(
        structure.lattice,
        [structure[i].species_string for i in order],
        [structure[i].frac_coords for i in order],
        validate_proximity=False,
        to_unit_cell=False,
    )
    ordered.to(filename=path, fmt="poscar")
    return path


def _shell_dir_name(label: str) -> str:
    return label


def export_pairs(
    result: PairsResult,
    out_root: str,
    write_all: bool = False,
) -> Dict[str, List[str]]:
    """Write pair structures into ``out_root`` (see module docstring).

    Returns a mapping ``combo -> written representative paths``.
    """
    written: Dict[str, List[str]] = {}
    os.makedirs(out_root, exist_ok=True)

    for label, base in result.bases.items():
        folder = os.path.join(out_root, "single_defect", label)
        os.makedirs(folder, exist_ok=True)
        write_structure(base.structure, os.path.join(folder, "POSCAR"))

    for combo, comp in result.complexes.items():
        paths: List[str] = []
        for shell in comp.shells:
            folder = os.path.join(out_root, combo, _shell_dir_name(shell.label))
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
) -> Dict[str, List[str]]:
    """Write triple structures into ``out_root`` (see module docstring)."""
    written: Dict[str, List[str]] = {}
    os.makedirs(out_root, exist_ok=True)

    for combo, tc in result.complexes.items():
        paths: List[str] = []
        for sel in tc.selections:
            cand = sel.candidate
            if cand is None:
                continue
            folder = os.path.join(out_root, combo, _shell_dir_name(sel.criterion))
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


def to_dataframe(rows: Sequence[dict]):
    """Return rows as a pandas DataFrame (pandas must be installed)."""
    try:
        import pandas as pd
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "pandas is required for to_dataframe(); "
            "install it or use write_rows_csv()"
        ) from exc
    return pd.DataFrame(rows)
