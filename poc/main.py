from dataclasses import dataclass
from typing import List

import residue_constants as rc

@dataclass
class Atom:
    atom_id: str
    comp_id: str
    seq_id: int
    x: float
    y: float
    z: float

    def __str__(self):
        return f"{self.atom_id} {self.comp_id} {self.seq_id} {self.x:.2f} {self.y:.2f} {self.z:.2f}"


def transform_sequence_to_line_structure(sequence: str) -> List[Atom]:
    atoms = []
    ca_pos = (0.0, 0.0, 0.0)
    seq_id = 0
    for aa in sequence:
        restype_3 = rc.restype_1to3[aa]
        atom_positions = rc.rigid_group_atom_positions[restype_3]
        for atom_position in atom_positions:
            position = atom_position[2]
            atoms.append(Atom(atom_position[0], restype_3, seq_id, position[0]+ca_pos[0], position[1]+ca_pos[1], position[2]+ca_pos[2]))
        ca_pos = (ca_pos[0] + rc.ca_ca, ca_pos[1], ca_pos[2])
        seq_id += 1
    return atoms


def atoms_to_cif(atoms: List[Atom], file_name) -> None:
    with open(file_name, 'w') as f:
        f.write(f"loop_\n")
        f.write(f"_atom_site.label_atom_id\n")
        f.write(f"_atom_site.label_comp_id\n")
        f.write(f"_atom_site.label_seq_id\n")
        f.write(f"_atom_site.Cartn_x\n")
        f.write(f"_atom_site.Cartn_y\n")
        f.write(f"_atom_site.Cartn_z\n")
        f.write('\n'.join(map(str, atoms)))


def main():
    structure = transform_sequence_to_line_structure("PEPTIDE")
    atoms_to_cif(structure, "test.cif")


if __name__ == '__main__':
    main()