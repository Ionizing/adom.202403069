#!/usr/bin/env python3

import numpy as np
from ase.io import read
from ase.build import make_supercell

atoms = read("POSCAR_prim.vasp")
M = np.array([[4, 2, 0],
              [0, 3, 0],
              [0, 0, 1]])

atoms_sup = make_supercell(atoms, M)
atoms_sup = atoms_sup[atoms_sup.numbers.argsort()[::-1]]
atoms_sup.write("POSCAR_sup.vasp", vasp5=True, direct=True)
