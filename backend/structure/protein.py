import time

import numpy as np

from io import StringIO
from itertools import count
from random import randint
from openmm import VerletIntegrator
from openmm.app import Simulation, NoCutoff, PDBxFile, ForceField
from openmm.unit import pico, kilojoule_per_mole
from typing import List, Optional, Dict

from backend.structure.amino_acid import AminoAcid
from backend.structure.types import AngleList


class Protein:
    TIME_STRUCTURE = []
    TIME_CIF = []
    TIME_FITNESS = []

    ANGLE_MIN = -180
    ANGLE_MAX = 180
    FORCE_FIELDS = {
        'amber': ForceField('amber14-all.xml'),
        'charmm': ForceField('charmm36.xml'),
    }

    SIMULATION_CACHE: Dict[str, Simulation] = {}

    def __init__(self, sequence: str, force_field: str='amber', *, angles: Optional[AngleList]=None, sigma: Optional[np.array]=None, flat_angles: Optional[np.array]=None):
        self._sequence:       str                = sequence
        self._force_field:    str                = force_field
        self._amino_acids:    List[AminoAcid]    = []
        self._atom_positions: List               = []
        self._fitness:        Optional[float]    = None
        self._cif_str:        Optional[str]      = None
        self._sigma:          Optional[np.array] = sigma # for self-adaptive evolution strategy

        assert not angles or len(angles) == len(sequence)

        if angles is not None:
            self._angles = angles
        if flat_angles is not None:
            self._angles = [(flat_angles[i], flat_angles[i+1], 180) for i in range(0, len(flat_angles), 2)]
        if angles is None and flat_angles is None:
            self._angles = self._get_random_angles(sequence)

        start = time.process_time_ns()
        self._compute_structure()
        end_structure = time.process_time_ns()
        self._compute_cif()
        end_cif = time.process_time_ns()
        self._compute_fitness()
        end_fitness = time.process_time_ns()

        self.TIME_STRUCTURE.append(end_structure - start)
        self.TIME_CIF.append(end_cif - end_structure)
        self.TIME_FITNESS.append(end_fitness - end_cif)

    def _compute_structure(self):
        predecessor = None
        for i, (aa, angles) in enumerate(zip(self._sequence, self._angles)):
            current_aa = AminoAcid(aa, angles, predecessor, i == len(self._sequence) - 1)
            self._amino_acids.append(current_aa)
            predecessor = current_aa

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def angles(self) -> AngleList:
        return self._angles

    @property
    def angles_flat(self) -> np.array:
        return np.array([angle for angles in self.angles for angle in angles[:2]])

    @property
    def fitness(self) -> float:
        return self._fitness

    @property
    def cif_str(self) -> str:
        return self._cif_str

    @property
    def atom_positions(self):
        return self._atom_positions

    @property
    def sigma(self) -> np.array:
        return self._sigma

    @sigma.setter
    def sigma(self, value):
        self._sigma = value

    def _compute_atom_positions(self):
        aa_nr = count(start=1)
        atom_nr = count(start=1)
        for aa in self._amino_acids:
            self._atom_positions.append(
                [
                    next(aa_nr),
                    aa.three_letter_code,
                    [[next(atom_nr), atom.atom_id, round(atom.x, 3), round(atom.y, 3), round(atom.z, 3)] for atom in aa.atoms]
                ]
            )

    def _compute_cif(self) -> None:
        self._compute_atom_positions()
        self._cif_str = self._cif_header_to_str()
        self._cif_str += self._cif_seperator_to_str()
        self._cif_str += self._cif_positions_to_str()

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
        for aa in self._atom_positions:
            for atom in aa[2]:
                data = [
                    'ATOM',       # group_PDB (here always ATOM)
                    f'{atom[0]}', # id
                    atom[1][0],   # type_symbol (e.g. N, O, C)
                    atom[1],      # label_atom_id (e.g. CA, O, H1)
                    '.',          # label_alt_id (always .)
                    aa[1],        # label_comp_id
                    'C',          # label_asym_id (always C)
                    '3',          # label_entity_id (always 3)
                    f'{aa[0]}',   # label_seq_id
                    '?',          # pdbx_PDB_ins_code (always ?)
                    f'{atom[2]}', # Cartn_x
                    f'{atom[3]}', # Cartn_y
                    f'{atom[4]}', # Cartn_z
                    '1.00',       # occupancy (always 1.00)
                    '0.00',       # B_iso_or_equiv (always 0.00)
                    f'{aa[0]}',   # auth_seq_id (like label_seq_id)
                    'A',          # auth_asym_id (always A)
                    '1'           # pdbx_PDB_model_num (1 because we always have only one chain)
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
        pdbx = PDBxFile(StringIO(self._cif_str))
        simulation = self._get_simulation(pdbx)
        state = simulation.context.getState(getEnergy=True)
        self._fitness = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)

    def _get_simulation(self, pdbx: PDBxFile) -> Simulation:
        cache_key = f"{self._sequence}_{self._force_field}"
        if cache_key in Protein.SIMULATION_CACHE:
            simulation = Protein.SIMULATION_CACHE[cache_key]
            simulation.context.setPositions(pdbx.positions)
        else:
            force_field = self.FORCE_FIELDS[self._force_field]
            system = force_field.createSystem(pdbx.topology, nonbondedMethod=NoCutoff)
            integrator = VerletIntegrator(0.001 * pico.factor)
            simulation = Simulation(pdbx.topology, system, integrator)
            simulation.context.setPositions(pdbx.positions)
            Protein.SIMULATION_CACHE[cache_key] = simulation

            if len(Protein.SIMULATION_CACHE) > 100:
                keys_to_remove = list(Protein.SIMULATION_CACHE.keys())[:20]
                for key in keys_to_remove:
                    del Protein.SIMULATION_CACHE[key]
        return simulation