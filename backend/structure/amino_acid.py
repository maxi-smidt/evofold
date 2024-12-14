import numpy as np
import numpy.typing as npt
import backend.structure.residue_constants as rc

from dataclasses import dataclass, field
from typing import Union, List, Literal

from backend.structure.atom import Atom


@dataclass
class AminoAcid:
    sequence_id: int
    one_letter_code: str
    previous_aa: Union['AminoAcid', None]
    is_last: bool
    next_aa: Union['AminoAcid', None] = None
    three_letter_code: str = field(init=False)
    atoms: List[Atom] = field(init=False)

    def __post_init__(self):
        self.three_letter_code = rc.restype_1to3[self.one_letter_code]
        self.atoms = []
        atom_positions = rc.rigid_group_atom_positions[self.three_letter_code]
        atom_positions = atom_positions[:-1] if not self.is_last else atom_positions  # last on gets termination O

        n = self._calculate_n_pos()
        for atom_position in atom_positions:
            if atom_position[0] == 'N':
                self.atoms.append(n)
            else:
                x = atom_position[2][0] + n.x
                y = atom_position[2][1] + n.y
                z = atom_position[2][2] + n.z
                self.atoms.append(Atom(atom_position[0], x, y, z))

    def _calculate_n_pos(self) -> Atom:
        atom_id = 'N'
        if self.previous_aa is None:
            return Atom(atom_id, 0.0, 0.0, 0.0)

        prev_n  = self.previous_aa.get_back_bone_position('N')
        prev_ca = self.previous_aa.get_back_bone_position('CA')
        prev_c  = self.previous_aa.get_back_bone_position('C')

        nc = prev_c - prev_n
        u_nc = nc / np.linalg.norm(nc)
        u_cn = u_nc

        cca = prev_ca - prev_c
        u_cca = cca / np.linalg.norm(cca)

        u_cn = self.adjust_angle(u_cca, u_cn, 114)

        cn = prev_c + u_cn * rc.n_c

        return Atom(atom_id, cn[0], cn[1], cn[2])

    @staticmethod
    def adjust_angle(u_cca, u_cn, alpha):
        alpha = np.deg2rad(alpha)

        cos_theta = np.dot(u_cca, u_cn)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        theta = np.arccos(cos_theta)  # current angle

        delta_alpha = alpha - theta  # rotation angle

        r = np.cross(u_cca, u_cn)  # rotation axis

        r_unit = r / np.linalg.norm(r)

        u_cn_rotated = (
                u_cn * np.cos(delta_alpha) +
                np.cross(r_unit, u_cn) * np.sin(delta_alpha) +
                r_unit * np.dot(r_unit, u_cn) * (1 - np.cos(delta_alpha))
        )
        return u_cn_rotated


    def get_back_bone_position(self, element: Literal['CA', 'C', 'N']) -> npt.ArrayLike:
        if element == 'CA':
            return self.atoms[1].get_position()
        elif element == 'C':
            return self.atoms[2].get_position()
        elif element == 'N':
            return self.atoms[0].get_position()