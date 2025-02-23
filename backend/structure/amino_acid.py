import numpy as np

import backend.structure.residue_constants as rc

from typing import Dict, Optional

from backend.structure.atom import Atom
from backend.structure.types import Position, Angle


class AminoAcid:
    def __init__(self, sequence_id: int, one_letter_code: str, angles: Angle, predecessor: Optional["AminoAcid"], is_last: bool) -> None:
        self._sequence_id = sequence_id
        self._three_letter_code = rc.restype_1to3[one_letter_code]
        self._angles = angles
        self._predecessor = predecessor
        self._is_last = is_last

        self._compute_structure()

    @property
    def three_letter_code(self) -> str:
        return self._three_letter_code

    @property
    def n(self):
        return self.atoms[0].get_position()

    @property
    def ca(self):
        return self.atoms[1].get_position()

    @property
    def c(self):
        return self.atoms[2].get_position()

    @property
    def phi(self) -> float:
        return self._angles[0]

    @property
    def psi(self) -> float:
        return self._angles[1]

    def _compute_structure(self):
        atom_positions = rc.rigid_group_atom_positions[self._three_letter_code].copy()

        if not self._is_last:  # only last amino acid gets termination oxygen
            del atom_positions['OXT']

        if self._predecessor:  # if it is not the first amino acid, the position and rotation has to be calculated
            atom_positions = self._transform_atom_positions(atom_positions)

        self.atoms = [Atom(atom_id, float(pos[0]), float(pos[1]), float(pos[2])) for atom_id, pos in atom_positions.items()]

    @staticmethod
    def _build_atom(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, bond_length: float, bond_angle: float, dihedral: float) -> np.ndarray:
        theta = np.radians(bond_angle)
        tau = np.radians(dihedral)

        v1 = p2 - p1
        v2 = p3 - p2

        e1 = v2 / np.linalg.norm(v2)

        n = np.cross(v1, v2)
        n /= np.linalg.norm(n)
        e2 = np.cross(n, e1)

        p4 = p3 + bond_length * (-np.cos(theta) * e1 + np.sin(theta) * (np.cos(tau) * e2 + np.sin(tau) * n))
        return p4

    def _transform_atom_positions(self, atom_positions: Dict[str, Position]) -> Dict[str, np.ndarray]:
        prev_N = np.array(self._predecessor.n)
        prev_CA = np.array(self._predecessor.ca)
        prev_C = np.array(self._predecessor.c)

        N_new = self._build_atom(prev_N, prev_CA, prev_C, bond_length=rc.n_c, bond_angle=116.2, dihedral=self._predecessor.psi)
        CA_new = self._build_atom(prev_CA, prev_C, N_new, bond_length=1.458, bond_angle=121.7, dihedral=0)
        C_new = self._build_atom(prev_C, N_new, CA_new, bond_length=1.525, bond_angle=110.4, dihedral=self.phi)

        transformed = {'N': N_new, 'CA': CA_new, 'C': C_new}

        template_N = np.array(atom_positions['N'])
        template_CA = np.array(atom_positions['CA'])
        template_C = np.array(atom_positions['C'])
        t_e1 = template_CA - template_N
        t_e1 /= np.linalg.norm(t_e1)
        t_e2 = template_C - template_CA
        t_e2 /= np.linalg.norm(t_e2)
        t_e3 = np.cross(t_e1, t_e2)

        g_e1 = CA_new - N_new
        g_e1 /= np.linalg.norm(g_e1)
        g_e2 = C_new - CA_new
        g_e2 /= np.linalg.norm(g_e2)
        g_e3 = np.cross(g_e1, g_e2)

        T = np.column_stack((g_e1, g_e2, g_e3)) @ np.column_stack((t_e1, t_e2, t_e3)).T
        translation = N_new - T @ template_N

        for atom_id, pos in atom_positions.items():
            if atom_id not in ['N', 'CA', 'C']:
                pos_arr = np.array(pos)
                new_pos = T @ pos_arr + translation
                transformed[atom_id] = new_pos

        return transformed