import numpy as np
import numpy.typing as npt
import backend.structure.residue_constants as rc

from dataclasses import dataclass, field
from scipy.optimize import minimize
from typing import Union, List, Literal

from backend.structure.atom import Atom


@dataclass
class AminoAcid:
    sequence_id: int
    one_letter_code: str
    previous_aa: Union['AminoAcid', None]
    next_aa: Union['AminoAcid', None] = None
    three_letter_code: str = field(init=False)
    atoms: List[Atom] = field(init=False)

    def __post_init__(self):
        self.three_letter_code = rc.restype_1to3[self.one_letter_code]
        self.atoms = []
        atom_positions = rc.rigid_group_atom_positions[self.three_letter_code]
        ca = self._calculate_ca_pos(np.array(atom_positions[0][2]))
        for atom_position in atom_positions:
            if atom_position[0] == 'CA':
                self.atoms.append(ca)
            else:
                R_x_180 = np.array([
                    [ 1,  0,  0],
                    [ 0, -1,  0],
                    [ 0,  0, -1]
                ])

                x = atom_position[2][0] + ca.x
                y = atom_position[2][1] + ca.y
                z = atom_position[2][2] + ca.z
                self.atoms.append(Atom(atom_position[0], x, y, z))

    def _calculate_ca_pos(self, d_n) -> Atom:
        atom_id = 'CA'
        if self.previous_aa is None:
            return Atom(atom_id, 0.0, 0.0, 0.0)

        initial_guess = self.previous_aa.get_back_bone('CA').get_position() + np.array([rc.ca_ca, 0.0, 0.0])
        result = minimize(self._objective_start_position, initial_guess, d_n)
        return Atom(atom_id, result.x[0], result.x[1], result.x[2])

    def _objective_start_position(self, r_ca_2: npt.ArrayLike, d_n: npt.ArrayLike):
        r_prev_ca = self.previous_aa.get_back_bone('CA').get_position()
        r_prev_c = self.previous_aa.get_back_bone('C').get_position()
        r_n_2 = r_ca_2 + d_n
        constraint1 = np.linalg.norm(r_ca_2 - r_prev_ca) - rc.ca_ca
        constraint2 = np.linalg.norm(r_n_2 - r_prev_c) - rc.n_c

        v1 = r_prev_ca - r_prev_c
        v2 = r_n_2 - r_prev_c
        angle_cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        target_angle_cos = np.cos(np.radians(111))  # Target ~111°
        constraint3 = angle_cos - target_angle_cos

        return constraint1 ** 2 + constraint2 ** 2 + constraint3 ** 2

    def get_back_bone(self, element: Literal['CA', 'C', 'N']):
        if element == 'CA':
            return self.atoms[1]
        elif element == 'C':
            return self.atoms[2]
        elif element == 'N':
            return self.atoms[0]
        else:
            raise NotImplementedError