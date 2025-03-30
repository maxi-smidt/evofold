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

    @property
    def omega(self) -> float:
        return self._angles[2]

    def _compute_structure(self):
        atom_positions = rc.rigid_group_atom_positions[self._three_letter_code].copy()

        if self._predecessor:  # if it is not the first amino acid, the position and rotation has to be calculated
            if 'H' in atom_positions: del atom_positions['H']
            del atom_positions['H2']
            del atom_positions['H3']
            atom_positions = self._transform_atom_positions(atom_positions)

        self.atoms = [Atom(atom_id, float(pos[0]), float(pos[1]), float(pos[2])) for atom_id, pos in atom_positions.items()]

        if self._is_last:
            self._build_oxygens()

    @staticmethod
    def _build_atom(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, bond_length: float, bond_angle: float, dihedral_angle: float) -> np.ndarray:
        bond_angle_rad = np.radians(bond_angle)
        dihedral_angle_rad = np.radians(dihedral_angle)

        u_z = (p2 - p3) / np.linalg.norm(p2 - p3)

        BA = p1 - p2
        v_perp = BA - np.dot(BA, u_z) * u_z
        u_x = v_perp / np.linalg.norm(v_perp)
        u_y = np.cross(u_z, u_x)

        x = bond_length * np.sin(bond_angle_rad) * np.cos(dihedral_angle_rad)
        y = bond_length * np.sin(bond_angle_rad) * np.sin(dihedral_angle_rad)
        z = bond_length * np.cos(bond_angle_rad)

        return p3 + x * u_x + y * u_y + z * u_z

    def _transform_atom_positions(self, atom_positions: Dict[str, Position]) -> Dict[str, np.ndarray]:
        prev_n  = np.array(self._predecessor.n)
        prev_ca = np.array(self._predecessor.ca)
        prev_c  = np.array(self._predecessor.c)

        n_new   = self._build_atom(prev_n,  prev_ca, prev_c, bond_length=rc.n_c, bond_angle=rc.ca_c_n_angle, dihedral_angle=-self._predecessor.psi)
        ca_new  = self._build_atom(prev_ca, prev_c,  n_new,  bond_length=1.458,  bond_angle=rc.c_n_ca_angle, dihedral_angle=self.omega)
        c_new   = self._build_atom(prev_c,  n_new,   ca_new, bond_length=1.525,  bond_angle=rc.n_ca_c_angle, dihedral_angle=-self.phi)

        self._build_o(self._predecessor, prev_ca, prev_c, n_new)

        local_points = np.array([atom_positions['N'], atom_positions['CA'], atom_positions['C']])
        absolute_points = np.array([n_new, ca_new, c_new])

        rotation, translation = self._kabsch_algorithm(local_points, absolute_points)

        transformed = {'N': n_new, 'CA': ca_new, 'C': c_new}
        for atom_id, pos in atom_positions.items():
            if atom_id not in ['N', 'CA', 'C']:
                absolute_coord = np.dot(rotation, pos) + translation
                transformed[atom_id] = tuple(absolute_coord.round(3))
        if self.three_letter_code != 'PRO':
            transformed['H'] = self._build_specific_atom(prev_c, n_new, ca_new, rc.c_n_h_angle, rc.n_h)
        return transformed

    @staticmethod
    def _build_o(amino_acid: "AminoAcid", prev_ca: np.array, prev_c: np.array, n_new: np.array, angle: float=rc.ca_c_o_angle, label: str='O') -> None:
        o_position = AminoAcid._build_specific_atom(prev_ca, prev_c, n_new, angle, rc.c_o)
        amino_acid.atoms.append(Atom(label, *o_position))

    @staticmethod
    def _build_specific_atom(a1: np.array, a2: np.array, a3: np.array, angle: float, distance: float) -> np.ndarray:
        angle = np.radians(angle)

        vec_a1 = a1 - a2
        vec_a3 = a3 - a2

        u = vec_a1 / np.linalg.norm(vec_a1)
        proj_n_on_u = np.dot(vec_a3, u) * u
        w = vec_a3 - proj_n_on_u
        w = -w / np.linalg.norm(w)

        direction = distance * (np.cos(angle) * u + np.sin(angle) * w)
        return a2 + direction
    
    def _build_oxygens(self) -> None:
        n_new = self._build_atom(self.n,  self.ca, self.c, bond_length=rc.n_c, bond_angle=rc.ca_c_n_angle, dihedral_angle=180)
        self._build_o(self, self.ca, self.c, n_new)
        self._build_o(self, self.ca, self.c, n_new, angle=-rc.ca_c_o_angle, label='OXT')

    @staticmethod
    def _kabsch_algorithm(local_points, absolute_points):
        centroid_local = np.mean(local_points, axis=0)
        centroid_abs = np.mean(absolute_points, axis=0)

        local_centered = local_points - centroid_local
        abs_centered = absolute_points - centroid_abs

        H = np.dot(local_centered.T, abs_centered)

        U, S, Vt = np.linalg.svd(H)

        R = np.dot(Vt.T, U.T)

        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = np.dot(Vt.T, U.T)

        T = centroid_abs - np.dot(R, centroid_local)

        return R, T