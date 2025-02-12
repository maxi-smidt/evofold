# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Constants used in AlphaFold."""

import collections

ca_ca = 3.80209737096
n_c = 1.33

# Format: The list for each AA type contains chi1, chi2, chi3, chi4 in
# this order (or a relevant subset from chi1 onwards). ALA and GLY don't have
# chi angles so their chi angle lists are empty.
chi_angles_atoms = {
    'ALA': [],
    # Chi5 in arginine is always 0 +- 5 degrees, so ignore it.
    'ARG': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'],
            ['CB', 'CG', 'CD', 'NE'], ['CG', 'CD', 'NE', 'CZ']],
    'ASN': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
    'ASP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
    'CYS': [['N', 'CA', 'CB', 'SG']],
    'GLN': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'],
            ['CB', 'CG', 'CD', 'OE1']],
    'GLU': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'],
            ['CB', 'CG', 'CD', 'OE1']],
    'GLY': [],
    'HIS': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'ND1']],
    'ILE': [['N', 'CA', 'CB', 'CG1'], ['CA', 'CB', 'CG1', 'CD1']],
    'LEU': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
    'LYS': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'],
            ['CB', 'CG', 'CD', 'CE'], ['CG', 'CD', 'CE', 'NZ']],
    'MET': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'SD'],
            ['CB', 'CG', 'SD', 'CE']],
    'PHE': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
    'PRO': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD']],
    'SER': [['N', 'CA', 'CB', 'OG']],
    'THR': [['N', 'CA', 'CB', 'OG1']],
    'TRP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
    'TYR': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
    'VAL': [['N', 'CA', 'CB', 'CG1']],
}

# Atoms positions relative to the 8 rigid groups, defined by the pre-omega, phi,
# psi and chi angles:
# 0: 'backbone group',
# 1: 'pre-omega-group', (empty)
# 2: 'phi-group', (currently empty, because it defines only hydrogens)
# 3: 'psi-group',
# 4,5,6,7: 'chi1,2,3,4-group'
# The atom positions are relative to the axis-end-atom of the corresponding
# rotation axis. The x-axis is in direction of the rotation axis, and the y-axis
# is defined such that the dihedral-angle-defining atom (the last entry in
# chi_angles_atoms above) is in the xy-plane (with a positive y-coordinate).
# format: [atomname, group_idx, rel_position (relative to ca)]
rigid_group_atom_positions = {
    'ALA': [
        ['N',   0, (-0.525,  1.363,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.526, -0.000, -0.000)],
        ['CB',  0, (-0.529, -0.774, -1.205)],
        # ['O',   3, ( 0.627,  1.062,  0.000)],
    ],
    'ARG': [
        ['N',   0, (-0.524,  1.362, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.525, -0.000, -0.000)],
        ['CB',  0, (-0.524, -0.778, -1.209)],
        # ['O',   3, ( 0.626,  1.062,  0.000)],
        # ['CG',  4, ( 0.616,  1.390, -0.000)],
        # ['CD',  5, ( 0.564,  1.414,  0.000)],
        # ['NE',  6, ( 0.539,  1.357, -0.000)],
        # ['NH1', 7, ( 0.206,  2.301,  0.000)],
        # ['NH2', 7, ( 2.078,  0.978, -0.000)],
        # ['CZ',  7, ( 0.758,  1.093, -0.000)],
    ],
    'ASN': [
        ['N',   0, (-0.536,  1.357,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.526, -0.000, -0.000)],
        ['CB',  0, (-0.531, -0.787, -1.200)],
        # ['O',   3, ( 0.625,  1.062,  0.000)],
        # ['CG',  4, ( 0.584,  1.399,  0.000)],
        # ['ND2', 5, ( 0.593, -1.188,  0.001)],
        # ['OD1', 5, ( 0.633,  1.059,  0.000)],
    ],
    'ASP': [
        ['N',   0, (-0.525,  1.362, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.527,  0.000, -0.000)],
        ['CB',  0, (-0.526, -0.778, -1.208)],
        # ['O',   3, ( 0.626,  1.062, -0.000)],
        # ['CG',  4, ( 0.593,  1.398, -0.000)],
        # ['OD1', 5, ( 0.610,  1.091,  0.000)],
        # ['OD2', 5, ( 0.592, -1.101, -0.003)],
    ],
    'CYS': [
        ['N',   0, (-0.522,  1.362, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.524,  0.000,  0.000)],
        ['CB',  0, (-0.519, -0.773, -1.212)],
        # ['O',   3, ( 0.625,  1.062, -0.000)],
        # ['SG',  4, ( 0.728,  1.653,  0.000)],
    ],
    'GLN': [
        ['N',   0, (-0.526,  1.361, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.526,  0.000,  0.000)],
        ['CB',  0, (-0.525, -0.779, -1.207)],
        # ['O',   3, ( 0.626,  1.062, -0.000)],
        # ['CG',  4, ( 0.615,  1.393,  0.000)],
        # ['CD',  5, ( 0.587,  1.399, -0.000)],
        # ['NE2', 6, ( 0.593, -1.189, -0.001)],
        # ['OE1', 6, ( 0.634,  1.060,  0.000)],
    ],
    'GLU': [
        ['N',   0, (-0.528,  1.361,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.526, -0.000, -0.000)],
        ['CB',  0, (-0.526, -0.781, -1.207)],
        # ['O',   3, ( 0.626,  1.062,  0.000)],
        # ['CG',  4, ( 0.615,  1.392,  0.000)],
        # ['CD',  5, ( 0.600,  1.397,  0.000)],
        # ['OE1', 6, ( 0.607,  1.095, -0.000)],
        # ['OE2', 6, ( 0.589, -1.104, -0.001)],
    ],
    'GLY': [
        ['N',   0, (-0.572,  1.337,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.517, -0.000, -0.000)],
        # ['O',   3, ( 0.626,  1.062, -0.000)],
    ],
    'HIS': [
        ['N',   0, (-0.527,  1.360,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.525,  0.000,  0.000)],
        ['CB',  0, (-0.525, -0.778, -1.208)],
        # ['O',   3, ( 0.625,  1.063,  0.000)],
        # ['CG',  4, ( 0.600,  1.370, -0.000)],
        # ['CD2', 5, ( 0.889, -1.021,  0.003)],
        # ['ND1', 5, ( 0.744,  1.160, -0.000)],
        # ['CE1', 5, ( 2.030,  0.851,  0.002)],
        # ['NE2', 5, ( 2.145, -0.466,  0.004)],
    ],
    'ILE': [
        ['N',   0, (-0.493,  1.373, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.527, -0.000, -0.000)],
        ['CB',  0, (-0.536, -0.793, -1.213)],
        # ['O',   3, ( 0.627,  1.062, -0.000)],
        # ['CG1', 4, ( 0.534,  1.437, -0.000)],
        # ['CG2', 4, ( 0.540, -0.785, -1.199)],
        # ['CD1', 5, ( 0.619,  1.391,  0.000)],
    ],
    'LEU': [
        ['N',   0, (-0.520,  1.363,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.525, -0.000, -0.000)],
        ['CB',  0, (-0.522, -0.773, -1.214)],
        # ['O',   3, ( 0.625,  1.063, -0.000)],
        # ['CG',  4, ( 0.678,  1.371,  0.000)],
        # ['CD1', 5, ( 0.530,  1.430, -0.000)],
        # ['CD2', 5, ( 0.535, -0.774,  1.200)],
    ],
    'LYS': [
        ['N',   0, (-0.526,  1.362, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.526,  0.000,  0.000)],
        ['CB',  0, (-0.524, -0.778, -1.208)],
        # ['O',   3, ( 0.626,  1.062, -0.000)],
        # ['CG',  4, ( 0.619,  1.390,  0.000)],
        # ['CD',  5, ( 0.559,  1.417,  0.000)],
        # ['CE',  6, ( 0.560,  1.416,  0.000)],
        # ['NZ',  7, ( 0.554,  1.387,  0.000)],
    ],
    'MET': [
        ['N',   0, (-0.521,  1.364, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.525,  0.000,  0.000)],
        ['CB',  0, (-0.523, -0.776, -1.210)],
        # ['O',   3, ( 0.625,  1.062, -0.000)],
        # ['CG',  4, ( 0.613,  1.391, -0.000)],
        # ['SD',  5, ( 0.703,  1.695,  0.000)],
        # ['CE',  6, ( 0.320,  1.786, -0.000)],
    ],
    'PHE': [
        ['N',   0, (-0.518,  1.363,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.524,  0.000, -0.000)],
        ['CB',  0, (-0.525, -0.776, -1.212)],
        # ['O',   3, ( 0.626,  1.062, -0.000)],
        # ['CG',  4, ( 0.607,  1.377,  0.000)],
        # ['CD1', 5, ( 0.709,  1.195, -0.000)],
        # ['CD2', 5, ( 0.706, -1.196,  0.000)],
        # ['CE1', 5, ( 2.102,  1.198, -0.000)],
        # ['CE2', 5, ( 2.098, -1.201, -0.000)],
        # ['CZ',  5, ( 2.794, -0.003, -0.001)],
    ],
    'PRO': [
        ['N',   0, (-0.566,  1.351, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.527, -0.000,  0.000)],
        ['CB',  0, (-0.546, -0.611, -1.293)],
        # ['O',   3, ( 0.621,  1.066,  0.000)],
        # ['CG',  4, ( 0.382,  1.445,  0.000)],
        # ['CD',  5, ( 0.427,  1.440,  0.000)],
        # ['CD',  5, ( 0.477,  1.424,  0.000)],  # manually made angle 2 degrees larger
    ],
    'SER': [
        ['N',   0, (-0.529,  1.360, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.525, -0.000, -0.000)],
        ['CB',  0, (-0.518, -0.777, -1.211)],
        # ['O',   3, ( 0.626,  1.062, -0.000)],
        # ['OG',  4, ( 0.503,  1.325,  0.000)],
    ],
    'THR': [
        ['N',   0, (-0.517,  1.364,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.526,  0.000, -0.000)],
        ['CB',  0, (-0.516, -0.793, -1.215)],
        # ['O',   3, ( 0.626,  1.062,  0.000)],
        # ['CG2', 4, ( 0.550, -0.718, -1.228)],
        # ['OG1', 4, ( 0.472,  1.353,  0.000)],
    ],
    'TRP': [
        ['N',   0, (-0.521,  1.363,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.525, -0.000,  0.000)],
        ['CB',  0, (-0.523, -0.776, -1.212)],
        # ['O',   3, ( 0.627,  1.062,  0.000)],
        # ['CG',  4, ( 0.609,  1.370, -0.000)],
        # ['CD1', 5, ( 0.824,  1.091,  0.000)],
        # ['CD2', 5, ( 0.854, -1.148, -0.005)],
        # ['CE2', 5, ( 2.186, -0.678, -0.007)],
        # ['CE3', 5, ( 0.622, -2.530, -0.007)],
        # ['NE1', 5, ( 2.140,  0.690, -0.004)],
        # ['CH2', 5, ( 3.028, -2.890, -0.013)],
        # ['CZ2', 5, ( 3.283, -1.543, -0.011)],
        # ['CZ3', 5, ( 1.715, -3.389, -0.011)],
    ],
    'TYR': [
        ['N',   0, (-0.522,  1.362,  0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.524, -0.000, -0.000)],
        ['CB',  0, (-0.522, -0.776, -1.213)],
        # ['O',   3, ( 0.627,  1.062, -0.000)],
        # ['CG',  4, ( 0.607,  1.382, -0.000)],
        # ['CD1', 5, ( 0.716,  1.195, -0.000)],
        # ['CD2', 5, ( 0.713, -1.194, -0.001)],
        # ['CE1', 5, ( 2.107,  1.200, -0.002)],
        # ['CE2', 5, ( 2.104, -1.201, -0.003)],
        # ['OH',  5, ( 4.168, -0.002, -0.005)],
        # ['CZ',  5, ( 2.791, -0.001, -0.003)],
    ],
    'VAL': [
        ['N',   0, (-0.494,  1.373, -0.000)],
        ['CA',  0, ( 0.000,  0.000,  0.000)],
        ['C',   0, ( 1.527, -0.000, -0.000)],
        ['CB',  0, (-0.533, -0.795, -1.213)],
        # ['O',   3, ( 0.627,  1.062, -0.000)],
        # ['CG1', 4, ( 0.540,  1.429, -0.000)],
        # ['CG2', 4, ( 0.533, -0.776,  1.203)],
    ],
}

# A list of atoms (excluding hydrogen) for each AA type. PDB naming convention.
residue_atoms = {
    'ALA': ['C', 'CA', 'CB', 'N', 'O'],
    'ARG': ['C', 'CA', 'CB', 'CG', 'CD', 'CZ', 'N', 'NE', 'O', 'NH1', 'NH2'],
    'ASP': ['C', 'CA', 'CB', 'CG', 'N', 'O', 'OD1', 'OD2'],
    'ASN': ['C', 'CA', 'CB', 'CG', 'N', 'ND2', 'O', 'OD1'],
    'CYS': ['C', 'CA', 'CB', 'N', 'O', 'SG'],
    'GLU': ['C', 'CA', 'CB', 'CG', 'CD', 'N', 'O', 'OE1', 'OE2'],
    'GLN': ['C', 'CA', 'CB', 'CG', 'CD', 'N', 'NE2', 'O', 'OE1'],
    'GLY': ['C', 'CA', 'N', 'O'],
    'HIS': ['C', 'CA', 'CB', 'CG', 'CD2', 'CE1', 'N', 'ND1', 'NE2', 'O'],
    'ILE': ['C', 'CA', 'CB', 'CG1', 'CG2', 'CD1', 'N', 'O'],
    'LEU': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'N', 'O'],
    'LYS': ['C', 'CA', 'CB', 'CG', 'CD', 'CE', 'N', 'NZ', 'O'],
    'MET': ['C', 'CA', 'CB', 'CG', 'CE', 'N', 'O', 'SD'],
    'PHE': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'N', 'O'],
    'PRO': ['C', 'CA', 'CB', 'CG', 'CD', 'N', 'O'],
    'SER': ['C', 'CA', 'CB', 'N', 'O', 'OG'],
    'THR': ['C', 'CA', 'CB', 'CG2', 'N', 'O', 'OG1'],
    'TRP': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'N', 'NE1', 'O'],
    'TYR': ['C', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'N', 'O', 'OH'],
    'VAL': ['C', 'CA', 'CB', 'CG1', 'CG2', 'N', 'O']
}

# Naming swaps for ambiguous atom names.
# Due to symmetries in the amino acids the naming of atoms is ambiguous in
# 4 of the 20 amino acids.
# (The LDDT paper lists 7 amino acids as ambiguous, but the naming ambiguities
# in LEU, VAL and ARG can be resolved by using the 3d constellations of
# the 'ambiguous' atoms and their neighbours)
residue_atom_renaming_swaps = {
    'ASP': {'OD1': 'OD2'},
    'GLU': {'OE1': 'OE2'},
    'PHE': {'CD1': 'CD2', 'CE1': 'CE2'},
    'TYR': {'CD1': 'CD2', 'CE1': 'CE2'},
}

# Van der Waals radii [Angstroem] of the atoms (from Wikipedia)
van_der_waals_radius = {
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'S': 1.8,
}

Bond = collections.namedtuple(
    'Bond', ['atom1_name', 'atom2_name', 'length', 'stddev'])
BondAngle = collections.namedtuple(
    'BondAngle',
    ['atom1_name', 'atom2_name', 'atom3name', 'angle_rad', 'stddev'])


# Between-residue bond lengths for general bonds (first element) and for Proline
# (second element).
between_res_bond_length_c_n = [1.329, 1.341]
between_res_bond_length_stddev_c_n = [0.014, 0.016]

# Between-residue cos_angles.
between_res_cos_angles_c_n_ca = [-0.5203, 0.0353]  # degrees: 121.352 +- 2.315
between_res_cos_angles_ca_c_n = [-0.4473, 0.0311]  # degrees: 116.568 +- 1.995


# A compact atom encoding with 14 columns
restype_name_to_atom14_names = {
    'ALA': ['N', 'CA', 'C', 'O', 'CB', '',    '',    '',    '',    '',    '',    '',    '',    ''],
    'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'NE',  'CZ',  'NH1', 'NH2', '',    '',    ''],
    'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'OD1', 'ND2', '',    '',    '',    '',    '',    ''],
    'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'OD1', 'OD2', '',    '',    '',    '',    '',    ''],
    'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG',  '',    '',    '',    '',    '',    '',    '',    ''],
    'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'OE1', 'NE2', '',    '',    '',    '',    ''],
    'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'OE1', 'OE2', '',    '',    '',    '',    ''],
    'GLY': ['N', 'CA', 'C', 'O', '',   '',    '',    '',    '',    '',    '',    '',    '',    ''],
    'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'ND1', 'CD2', 'CE1', 'NE2', '',    '',    '',    ''],
    'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', '',    '',    '',    '',    '',    ''],
    'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', '',    '',    '',    '',    '',    ''],
    'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  'CE',  'NZ',  '',    '',    '',    '',    ''],
    'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'SD',  'CE',  '',    '',    '',    '',    '',    ''],
    'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', 'CE1', 'CE2', 'CZ',  '',    '',    ''],
    'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD',  '',    '',    '',    '',    '',    '',    ''],
    'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG',  '',    '',    '',    '',    '',    '',    '',    ''],
    'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', '',    '',    '',    '',    '',    '',    ''],
    'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG',  'CD1', 'CD2', 'CE1', 'CE2', 'CZ',  'OH',  '',    ''],
    'VAL': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', '',    '',    '',    '',    '',    '',    '']
}

restype_1to3 = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'Q': 'GLN',
    'E': 'GLU',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
}

def atom_id_to_type(atom_id: str) -> str:
  if atom_id.startswith('C'):
    return 'C'
  elif atom_id.startswith('N'):
    return 'N'
  elif atom_id.startswith('O'):
    return 'O'
  elif atom_id.startswith('H'):
    return 'H'
  elif atom_id.startswith('S'):
    return 'S'
  raise ValueError('Atom ID not recognized.')


restype_3to1 = {v: k for k, v in restype_1to3.items()}