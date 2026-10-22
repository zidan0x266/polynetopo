"""The LAMMPS data reader shared by the spatial analyses."""

from __future__ import annotations

import numpy as np
import pytest

from graph_analysis.lammps import Direction, LammpsDataReader, LammpsFormatError

from conftest import write_lammps


def test_atom_ids_need_not_be_consecutive(tmp_path):
    coordinates = np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [3.0, 1.0, 1.0]])
    path = write_lammps(tmp_path / "a.data", coordinates, [(40, 7), (7, 12)], (10.0, 10.0, 10.0),
                        lo=(-5.0, 0.0, 0.0), atom_ids=[12, 7, 40])
    network = LammpsDataReader().read(path)
    assert network.atom_ids.tolist() == [7, 12, 40]  # sorted
    assert network.bond_atom_ids.tolist() == [[40, 7], [7, 12]]  # file order and endpoint order kept
    assert network.box_lengths.tolist() == [10.0, 10.0, 10.0]
    assert network.box_bounds[0].tolist() == [-5.0, 5.0]


def test_a_bond_to_an_unknown_atom_is_rejected(tmp_path):
    path = write_lammps(tmp_path / "a.data", np.zeros((2, 3)), [(1, 9)], (10.0, 10.0, 10.0))
    with pytest.raises(LammpsFormatError, match="unknown atom ID"):
        LammpsDataReader().read(path)


def test_triclinic_boxes_are_rejected(tmp_path):
    path = write_lammps(tmp_path / "a.data", np.zeros((2, 3)), [(1, 2)], (10.0, 10.0, 10.0))
    path.write_text(path.read_text().replace("zlo zhi\n", "zlo zhi\n0.0 0.0 0.0 xy xz yz\n"))
    with pytest.raises(LammpsFormatError, match="Triclinic"):
        LammpsDataReader().read(path)


def test_direction_parses_names_and_numbers():
    assert [Direction.parse(v) for v in ("x", "Y", "z", 0, 1, "2")] == [0, 1, 2, 0, 1, 2]
    with pytest.raises(ValueError):
        Direction.parse("w")
