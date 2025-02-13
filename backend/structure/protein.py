import io

from itertools import count
from random import randint
from typing import List, Optional, Tuple
from openmm import VerletIntegrator
from openmm.app import Simulation, ForceField, NoCutoff
from openmm.unit import pico, kilojoule_per_mole

from backend.structure.amino_acid import AminoAcid
from backend.structure.fixer.pdbfixer import PDBFixer

AngleList = List[Tuple[float, float, float]]

class Protein:
    ANGLE_MIN = -180
    ANGLE_MAX = 180
    FORCE_FIELD = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    def __init__(self, sequence: str, angles: Optional[AngleList] = None):
        self._sequence:     str             = sequence
        self._amino_acids:  List[AminoAcid] = []
        self._fitness:      Optional[float] = None
        self._cif_str:      Optional[str]   = None
        self._angles:       AngleList       = angles or self._get_random_angles(sequence)  # ϕ and ψ angles alternating

        self._compute_structure()

    def _compute_structure(self):
        previous_aa = None
        for i, aa in enumerate(self._sequence):
            current_aa = AminoAcid(i, aa, previous_aa, i == len(self._sequence) - 1)
            self._amino_acids.append(current_aa)
            previous_aa = current_aa
        previous_aa = None
        for aa in self._amino_acids[::-1]:
            if previous_aa is not None:
                aa.next_aa = previous_aa
            previous_aa = aa

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def angles(self) -> AngleList:
        return self._angles

    @property
    def fitness(self) -> float:
        if self._fitness is None:
            self._compute_fitness()
        return self._fitness

    def to_cif(self, file_name: Optional[str] = None) -> str:
        if self._cif_str is None:
            self._cif_str = self._cif_header_to_str()
            self._cif_str += self._cif_seperator_to_str()
            self._cif_str += self._cif_positions_to_str()

        if file_name is not None:
            with open(file_name, 'w') as f:
                f.write(self._cif_str)

        return self._cif_str

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
        for aa in self._amino_acids:
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
        return f'data_EvoFold_{self._sequence}\n'

    @staticmethod
    def _cif_seperator_to_str():
        return '#\n'

    @staticmethod
    def _get_random_angles(sequence: str) -> AngleList:
        rand = lambda: randint(Protein.ANGLE_MIN, Protein.ANGLE_MAX)
        return [(rand(), rand(), 180) for _ in range(len(sequence))]

    def _compute_fitness(self) -> None:
        cif_file = io.StringIO(self.to_cif())
        cif_file.seek(0)
        fixer = PDBFixer(cif_file)
        fixer.addMissingHydrogens()

        system = self.FORCE_FIELD.createSystem(fixer.topology, nonbondedMethod=NoCutoff)
        integrator = VerletIntegrator(0.001 * pico.factor)
        simulation = Simulation(fixer.topology, system, integrator)
        simulation.context.setPositions(fixer.positions)
        state = simulation.context.getState(getEnergy=True)
        self._fitness = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)