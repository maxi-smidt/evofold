import io

from dataclasses import dataclass, field
from itertools import count
from typing import List
from openmm import VerletIntegrator
from openmm.app import Simulation, ForceField, NoCutoff
from openmm.unit import pico, Quantity

from amino_acid import AminoAcid
from fixer.pdbfixer import PDBFixer


@dataclass
class ProteinStructure:
    sequence: str
    amino_acids: List[AminoAcid] = field(init=False)

    def __post_init__(self):
        previous_aa = None
        self.amino_acids = []
        for i, aa in enumerate(self.sequence):
            current_aa = AminoAcid(i, aa, previous_aa)
            self.amino_acids.append(current_aa)
            previous_aa = current_aa
        previous_aa = None
        for aa in self.amino_acids[::-1]:
            if previous_aa is not None:
                aa.next_aa = previous_aa
            previous_aa = aa



    def to_cif(self, file_name: str | None = None) -> str:
        cif_string = self._cif_header_to_str()
        cif_string += self._cif_seperator_to_str()
        cif_string += self._cif_positions_to_str()

        if file_name is not None:
            with open(file_name, 'w') as f:
                f.write(cif_string)

        return cif_string

    def _cif_positions_to_str(self) -> str:
        atom_header = [
            'loop_', '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol', '_atom_site.label_atom_id',
            '_atom_site.label_alt_id', '_atom_site.label_comp_id', '_atom_site.label_asym_id',
            '_atom_site.label_entity_id', '_atom_site.label_seq_id', '_atom_site.pdbx_PDB_ins_code',
            '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
            '_atom_site.B_iso_or_equiv', '_atom_site.auth_seq_id', '_atom_site.auth_asym_id',
            '_atom_site.pdbx_PDB_model_num'
        ]
        positions_string = '\n'.join(atom_header)
        atom_nr = count(start=1)
        aa_nr = count(start=1)
        for aa in self.amino_acids:
            seq_id = next(aa_nr)
            for atom in aa.atoms:
                data = [
                    'ATOM',                 # group_PDB (here always ATOM)
                    f'{next(atom_nr)}',     # id
                    atom.atom_id[0],        # type_symbol (e.g. N, O, C)
                    atom.atom_id,           # label_atom_id (e.g. CA, O, H1)
                    '.',                    # label_alt_id (always .)
                    aa.three_letter_code,   # label_comp_id
                    'C',                    # label_asym_id (always C)
                    '3',                    # label_entity_id (always 3)
                    f'{seq_id}',            # label_seq_id
                    '?',                    # pdbx_PDB_ins_code (always ?)
                    f'{atom.x:.3f}',        # Cartn_x
                    f'{atom.y:.3f}',        # Cartn_y
                    f'{atom.z:.3f}',        # Cartn_z
                    '1.00',                 # occupancy (always 1.00)
                    '0.00',                 # B_iso_or_equiv (always 0.00)
                    f'{seq_id}',            # auth_seq_id (like label_seq_id)
                    'A',                    # auth_asym_id (always A)
                    '1'                     # pdbx_PDB_model_num (1 because we always have only one chain)
                ]
                positions_string += '\n' + '\t'.join(data)
        return positions_string

    def _cif_header_to_str(self):
        return f'data_EvoFold_{self.sequence}\n'

    @staticmethod
    def _cif_seperator_to_str():
        return '#\n'

    def fitness(self) -> Quantity:
        cif_file = io.StringIO(self.to_cif())
        cif_file.seek(0)
        fixer = PDBFixer(cif_file)

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens()

        forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

        system = forcefield.createSystem(fixer.topology, nonbondedMethod=NoCutoff)
        integrator = VerletIntegrator(0.001 * pico.factor)
        simulation = Simulation(fixer.topology, system, integrator)
        simulation.context.setPositions(fixer.positions)
        state = simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy()