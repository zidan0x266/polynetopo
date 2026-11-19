"""LammpsDataReader: the atom styles it accepts give the same network."""

from __future__ import annotations

import numpy as np
import pytest

from graph_analysis.lammps import LammpsDataReader, LammpsFormatError

ATOMS = [(3, 1, 1, 2.5, 0.5, 1.0), (1, 1, 1, 0.5, 0.5, 1.0), (2, 1, 1, 1.5, 0.5, 1.0), (4, 2, 1, 3.5, 2.5, 1.0)]
BONDS = [(1, 1, 1, 2), (2, 1, 2, 3), (3, 1, 3, 4)]


def write(path, header, style, images=True):
    rows = []
    for atom_id, mol, atom_type, x, y, z in ATOMS:
        charge = ["0.0"] if style == "full" else []
        flags = ["0", "0", "0"] if images else []
        rows.append(" ".join([str(atom_id), str(mol), str(atom_type), *charge, str(x), str(y), str(z), *flags]))
    text = "\n".join([
        "LAMMPS data file", "", f"{len(ATOMS)} atoms", "1 atom types", f"{len(BONDS)} bonds", "1 bond types", "",
        "0.0 4.0 xlo xhi", "0.0 4.0 ylo yhi", "0.0 4.0 zlo zhi", "", "Masses", "", "1 1.0", "", header, "", *rows, "",
        "Bonds", "", *(" ".join(map(str, bond)) for bond in BONDS), ""])
    path.write_text(text)
    return path


@pytest.mark.parametrize("header, style, images", [
    ("Atoms # molecular", "molecular", True), ("Atoms # molecular", "molecular", False),
    ("Atoms # bond", "bond", True), ("Atoms # angle", "angle", False), ("Atoms", "full", True)])
def test_atom_styles_give_the_same_network(tmp_path, header, style, images):
    reference = LammpsDataReader().read(write(tmp_path / "full.data", "Atoms # full", "full"))
    network = LammpsDataReader().read(write(tmp_path / f"{style}.data", header, style, images))
    assert np.array_equal(network.atom_ids, [1, 2, 3, 4])
    assert np.array_equal(network.coordinates, reference.coordinates)
    assert np.array_equal(network.coordinates[2], [2.5, 0.5, 1.0])
    assert np.array_equal(network.edges, reference.edges)
    assert np.array_equal(network.box_bounds, reference.box_bounds)


def test_unsupported_atom_style(tmp_path):
    with pytest.raises(LammpsFormatError, match="Unsupported atom style 'atomic'"):
        LammpsDataReader().read(write(tmp_path / "atomic.data", "Atoms # atomic", "molecular"))


def test_rows_of_another_style(tmp_path):
    with pytest.raises(LammpsFormatError, match="6 or 9 columns"):
        LammpsDataReader().read(write(tmp_path / "mixed.data", "Atoms # molecular", "full"))
