import numpy as np
import numpy.typing as npt
import residue_constants as rc

from dataclasses import dataclass, field
from itertools import count
from scipy.optimize import minimize
from typing import List, Literal, Union, TextIO


@dataclass
class AtomStructure:
    atom_id: str
    x: float
    y: float
    z: float

    def get_position(self) -> npt.ArrayLike:
        return np.array([self.x, self.y, self.z])

@dataclass
class AminoAcidStructure:
    sequence_id: int
    one_letter_code: str
    previous_aa: Union['AminoAcidStructure', None]
    three_letter_code: str = field(init=False)
    atoms: List[AtomStructure] = field(init=False)

    def __post_init__(self):
        self.three_letter_code = rc.restype_1to3[self.one_letter_code]
        self.atoms = []
        atom_positions = rc.rigid_group_atom_positions[self.three_letter_code]
        ca = self._calculate_ca_pos(np.array(atom_positions[0][2]))
        for atom_position in atom_positions:
            if atom_position[0] == 'CA':
                self.atoms.append(ca)
            else:
                position = np.array(atom_position[2]) + ca.get_position()
                self.atoms.append(AtomStructure(atom_position[0], position[0], position[1], position[2]))

    def _calculate_ca_pos(self, d_n) -> AtomStructure:
        atom_id = 'CA'
        if self.previous_aa is None:
            return AtomStructure(atom_id, 0.0, 0.0, 0.0)

        initial_guess = self.previous_aa.get_back_bone('CA').get_position() + np.array([rc.ca_ca, 0.0, 0.0])
        result = minimize(self._objective_start_position, initial_guess, d_n)
        return AtomStructure(atom_id, result.x[0], result.x[1], result.x[2])

    def _objective_start_position(self, r_ca_2: npt.ArrayLike, d_n: npt.ArrayLike):
        r_prev_ca = self.previous_aa.get_back_bone('CA').get_position()
        r_prev_c = self.previous_aa.get_back_bone('C').get_position()
        r_n_2 = r_ca_2 + d_n
        constraint1 = np.linalg.norm(r_ca_2 - r_prev_ca) - rc.ca_ca
        constraint2 = np.linalg.norm(r_n_2 - r_prev_c) - rc.n_c
        return constraint1 ** 2 + constraint2 ** 2

    def get_back_bone(self, element: Literal['CA', 'C', 'N']):
        if element == 'CA':
            return self.atoms[1]
        elif element == 'C':
            return self.atoms[2]
        elif element == 'N':
            return self.atoms[0]
        else:
            raise NotImplementedError


@dataclass
class ProteinStructure:
    sequence: str
    amino_acids: List[AminoAcidStructure] = field(init=False)

    def __post_init__(self):
        previous_aa = None
        self.amino_acids = []
        for i, aa in enumerate(self.sequence):
            current_aa = AminoAcidStructure(i, aa, previous_aa)
            self.amino_acids.append(current_aa)
            previous_aa = current_aa

    def to_cif(self, file_name: str):
        with open(file_name, 'w') as f:
            self._print_cif_header(f)
            self._print_cif_positions(f)

    def _print_cif_positions(self, file: TextIO):
        atom_header = [
            'loop_', '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol', '_atom_site.label_atom_id',
            '_atom_site.label_alt_id', '_atom_site.label_comp_id', '_atom_site.label_asym_id',
            '_atom_site.label_entity_id', '_atom_site.label_seq_id', '_atom_site.pdbx_PDB_ins_code',
            '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
            '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge', '_atom_site.auth_seq_id',
            '_atom_site.auth_comp_id', '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
            '_atom_site.pdbx_PDB_model_num'
        ]

        file.write('\n'.join(atom_header))
        atom_nr = count(start=0)
        aa_nr = count(start=0)
        for aa in self.amino_acids:
            seq_id = next(aa_nr)
            for atom in aa.atoms:
                data = [
                    'ATOM',                 # group_PDB (here always ATOM)
                    f"{next(atom_nr)}",     # id
                    atom.atom_id[0],        # type_symbol (e.g. N, O, C)
                    atom.atom_id,           # label_atom_id (e.g. CA, O, H1)
                    '.',                    # label_alt_id (always .)
                    aa.three_letter_code,   # label_comp_id
                    'C',                    # label_asym_id (always C)
                    '3',                    # label_entity_id (always 3)
                    f"{seq_id}",            # label_seq_id
                    '?',                    # pdbx_PDB_ins_code (always ?)
                    f"{atom.x:.3f}",        # Cartn_x
                    f"{atom.y:.3f}",        # Cartn_y
                    f"{atom.z:.3f}",        # Cartn_z
                    '1.00',                 # occupancy (always 1.00)
                    '0.00',                 # B_iso_or_equiv (always 0.00)
                    '?',                    # pdbx_formal_charge (always ?)
                    f"{seq_id}",            # auth_seq_id (like label_seq_id)
                    aa.three_letter_code,   # auth_comp_id
                    'A',                    # auth_asym_id (always A)
                    f"{atom.atom_id}",      # auth_atom_id (like label_atom_id)
                    '1'                     # pdbx_PDB_model_num (1 because we always have only one chain)
                ]
                file.write('\n' + '\t'.join(data))

    def _print_cif_header(self, file: TextIO):
        file.write(f'data_EvoFold_{self.sequence}\n')

    @staticmethod
    def _print_cif_seperator(file: TextIO):
        file.write('#\n')


x = ProteinStructure("PEPTIDE")

x.to_cif('test_new.cif')




















