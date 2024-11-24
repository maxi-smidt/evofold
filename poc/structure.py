from dataclasses import dataclass, field
from typing import List

import residue_constants as rc


@dataclass
class AtomStructure:
    atom_id: str
    x: float
    y: float
    z: float


@dataclass()
class AminoAcidStructure:
    sequence_position: int
    one_letter_code: str
    following_aa: AtomStructure | None
    three_letter_code: str = field(init=False)
    atoms: List[AtomStructure] = field(init=False)

    def __post_init__(self):
        self.three_letter_code = rc.rigid_group_atom_positions[self.one_letter_code]
        self.atoms = []
        atom_positions = rc.rigid_group_atom_positions[self.three_letter_code]
        for atom_position in atom_positions:
            position = atom_position[2]
            self.atoms.append(AtomStructure(atom_position[0], position[0], position[1], position[2]))


@dataclass
class ProteinStructure:
    sequence: str
    amino_acids: List[AminoAcidStructure] = field(init=False)

    def __post_init__(self):
        following_aa = None
        self.amino_acids = []
        for i, aa in enumerate(self.sequence[::-1]):
            current_aa = AminoAcidStructure(i, aa, following_aa)
            self.amino_acids.append(current_aa)
            following_aa = current_aa

        self.amino_acids = self.amino_acids[::-1]

    def to_cif(self):
        ...