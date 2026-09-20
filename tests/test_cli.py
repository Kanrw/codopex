"""End-to-end tests of the ``codopex`` command line interface."""

import os

import pytest

import codopex as cp
from codopex.cli import main
from helpers import cubic_2x2x2


def _write(structure, path):
    structure.to(filename=str(path), fmt="poscar")
    return str(path)


def test_cli_version(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert cp.__version__ in capsys.readouterr().out


def test_cli_pairs_exports_and_saves(tmp_path, capsys):
    structure = _write(cubic_2x2x2(), tmp_path / "bulk.vasp")
    out = tmp_path / "out"
    save = tmp_path / "pairs.json.gz"
    rc = main(
        [
            "pairs",
            structure,
            "-r",
            "Li@Na",
            "--shells",
            "nnn",
            "--out",
            str(out),
            "--save",
            str(save),
        ]
    )
    assert rc == 0
    printed = capsys.readouterr().out
    assert "Li_Na_Oh" in printed
    assert os.path.exists(out / "Li_Na_Oh+Li_Na_Oh" / "nn" / "POSCAR")
    assert os.path.exists(out / "Li_Na_Oh+Li_Na_Oh" / "nnn" / "POSCAR")
    assert os.path.exists(out / "manifest.csv")

    loaded = cp.io.load(str(save))
    assert loaded.type_labels == ["Li_Na_Oh"]


def test_cli_triples_from_saved_pairs(tmp_path, capsys):
    pairs_path = tmp_path / "pairs.json"
    pairs = cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")], shells="nnn")
    cp.io.save(pairs, str(pairs_path))

    out = tmp_path / "triples"
    rc = main(
        [
            "triples",
            "--load-pairs",
            str(pairs_path),
            "--pair",
            "Li_Na_Oh+Li_Na_Oh",
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    printed = capsys.readouterr().out
    assert "Li_Na_Oh+Li_Na_Oh" in printed
    assert os.path.exists(out / "Li_Na_Oh+Li_Na_Oh+Li_Na_Oh" / "nearA_farB" / "POSCAR")
    assert os.path.exists(out / "manifest.csv")


def test_cli_triples_skips_pairs_without_shell(tmp_path, capsys):
    pairs_path = tmp_path / "pairs.json"
    cp.io.save(
        cp.generate_pairs(cubic_2x2x2(), [("Na", "Li")]),  # nn only
        str(pairs_path),
    )
    rc = main(["triples", "--load-pairs", str(pairs_path), "--pair-shell", "nnn"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "no 'nnn' representative" in err
    assert "no parent pair produced" in err


def test_cli_classify(tmp_path, capsys):
    bulk = cubic_2x2x2()
    defect = bulk.copy()
    defect.replace(0, "Li")
    rc = main(
        [
            "classify",
            _write(bulk, tmp_path / "bulk.vasp"),
            _write(defect, tmp_path / "defect.vasp"),
        ]
    )
    assert rc == 0
    printed = capsys.readouterr().out
    assert "substitution" in printed
    assert "Li_Na_Oh" in printed


def test_cli_rejects_bad_reaction(tmp_path):
    with pytest.raises(SystemExit) as exc:
        main(["pairs", str(tmp_path / "bulk.vasp"), "-r", "LiNa"])
    assert exc.value.code == 2


def test_cli_missing_structure_returns_error(tmp_path, capsys):
    rc = main(["pairs", str(tmp_path / "nope.vasp"), "-r", "Li@Na"])
    assert rc == 1
    assert "error" in capsys.readouterr().err
