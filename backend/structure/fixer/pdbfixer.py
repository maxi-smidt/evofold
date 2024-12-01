"""
pdbfixer.py: Fixes problems in PDB files

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2013-2024 Stanford University and the Authors.
Authors: Peter Eastman
Modified by: Maximilian Smidt
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import absolute_import
from dataclasses import dataclass
from io import StringIO
from typing import Literal, Optional
from urllib.request import urlopen

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
from openmm.app.element import hydrogen, oxygen
from openmm.app.forcefield import NonbondedGenerator

from openmm.app.internal.compiled import matchResidueToTemplate as matchResidue

import numpy as np
import numpy.linalg as lin
import sys
import os
import os.path
import math
from collections import defaultdict

substitutions = {
    '2AS':'ASP', '3AH':'HIS', '5HP':'GLU', '5OW':'LYS', 'ACL':'ARG', 'AGM':'ARG', 'AIB':'ALA', 'ALM':'ALA', 'ALO':'THR', 'ALY':'LYS', 'ARM':'ARG',
    'ASA':'ASP', 'ASB':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASQ':'ASP', 'AYA':'ALA', 'BCS':'CYS', 'BHD':'ASP', 'BMT':'THR', 'BNN':'ALA',
    'BUC':'CYS', 'BUG':'LEU', 'C5C':'CYS', 'C6C':'CYS', 'CAS':'CYS', 'CCS':'CYS', 'CEA':'CYS', 'CGU':'GLU', 'CHG':'ALA', 'CLE':'LEU', 'CME':'CYS',
    'CSD':'ALA', 'CSO':'CYS', 'CSP':'CYS', 'CSS':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'CXM':'MET', 'CY1':'CYS', 'CY3':'CYS', 'CYG':'CYS',
    'CYM':'CYS', 'CYQ':'CYS', 'DAH':'PHE', 'DAL':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DCY':'CYS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA',
    'DHI':'HIS', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLY':'LYS', 'DNP':'ALA', 'DPN':'PHE', 'DPR':'PRO', 'DSN':'SER', 'DSP':'ASP',
    'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'EFC':'CYS', 'FLA':'ALA', 'FME':'MET', 'GGL':'GLU', 'GL3':'GLY', 'GLZ':'GLY',
    'GMA':'GLU', 'GSC':'GLY', 'HAC':'ALA', 'HAR':'ARG', 'HIC':'HIS', 'HIP':'HIS', 'HMR':'ARG', 'HPQ':'PHE', 'HTR':'TRP', 'HYP':'PRO',
    'IAS':'ASP', 'IIL':'ILE', 'IYR':'TYR', 'KCX':'LYS', 'LLP':'LYS', 'LLY':'LYS', 'LTR':'TRP', 'LYM':'LYS', 'LYZ':'LYS', 'MAA':'ALA', 'MEN':'ASN',
    'MHS':'HIS', 'MIS':'SER', 'MK8':'LEU', 'MLE':'LEU', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MVA':'VAL', 'NEM':'HIS', 'NEP':'HIS', 'NLE':'LEU',
    'NLN':'LEU', 'NLP':'LEU', 'NMC':'GLY', 'OAS':'SER', 'OCS':'CYS', 'OMT':'MET', 'PAQ':'TYR', 'PCA':'GLU', 'PEC':'CYS', 'PHI':'PHE',
    'PHL':'PHE', 'PR3':'CYS', 'PRR':'ALA', 'PTR':'TYR', 'PYX':'CYS', 'SAC':'SER', 'SAR':'GLY', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS',
    'SEL':'SER', 'SEP':'SER', 'SET':'SER', 'SHC':'CYS', 'SHR':'LYS', 'SMC':'CYS', 'SOC':'CYS', 'STY':'TYR', 'SVA':'SER', 'TIH':'ALA',
    'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR'
}
proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
rnaResidues = ['A', 'G', 'C', 'U', 'I']
dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']

class Sequence(object):
    """Sequence holds the sequence of a chain, as specified by SEQRES records."""
    def __init__(self, chainId, residues):
        self.chainId = chainId
        self.residues = residues

class ModifiedResidue(object):
    """ModifiedResidue holds information about a modified residue, as specified by a MODRES record."""
    def __init__(self, chainId, number, residueName, standardName):
        self.chainId = chainId
        self.number = number
        self.residueName = residueName
        self.standardName = standardName

class Template:
    """Template represents a standard residue, or a nonstandard one registered with registerTemplate()."""
    def __init__(self, topology, positions, terminal=None):
        self.topology = topology
        self.positions = positions
        if terminal is None:
            terminal = [False]*topology.getNumAtoms()
        self.terminal = terminal

@dataclass
class CCDAtomDefinition:
    """
    Description of an atom in a residue from the Chemical Component Dictionary (CCD).
    """
    atomName: str
    symbol: str
    leaving: bool
    coords: mm.Vec3
    charge: int
    aromatic: bool

@dataclass
class CCDBondDefinition:
    """
    Description of a bond in a residue from the Chemical Component Dictionary (CCD).
    """
    atom1: str
    atom2: str
    order: Literal['SING', 'DOUB', 'TRIP', 'QUAD', 'AROM', 'DELO', 'PI', 'POLY']
    aromatic: bool

@dataclass
class CCDResidueDefinition:
    """
    Description of a residue from the Chemical Component Dictionary (CCD).
    """
    residueName: str
    atoms: list[CCDAtomDefinition]
    bonds: list[CCDBondDefinition]

    @classmethod
    def fromReader(cls, reader: PdbxReader) -> 'CCDResidueDefinition':
        """
        Create a CCDResidueDefinition by parsing a CCD CIF file.
        """
        data = []
        reader.read(data)
        block = data[0]

        residueName = block.getObj('chem_comp').getValue("id")

        descriptorsData = block.getObj("pdbx_chem_comp_descriptor")
        typeCol = descriptorsData.getAttributeIndex("type")

        atomData = block.getObj('chem_comp_atom')
        atomNameCol = atomData.getAttributeIndex('atom_id')
        symbolCol = atomData.getAttributeIndex('type_symbol')
        leavingCol = atomData.getAttributeIndex('pdbx_leaving_atom_flag')
        xCol = atomData.getAttributeIndex('pdbx_model_Cartn_x_ideal')
        yCol = atomData.getAttributeIndex('pdbx_model_Cartn_y_ideal')
        zCol = atomData.getAttributeIndex('pdbx_model_Cartn_z_ideal')
        chargeCol = atomData.getAttributeIndex('charge')
        aromaticCol = atomData.getAttributeIndex('pdbx_aromatic_flag')

        atoms = [
            CCDAtomDefinition(
                atomName=row[atomNameCol],
                symbol=row[symbolCol],
                leaving=row[leavingCol] == 'Y',
                coords=mm.Vec3(float(row[xCol]), float(row[yCol]), float(row[zCol]))*0.1,
                charge=row[chargeCol],
                aromatic=row[aromaticCol] == 'Y'
            ) for row in atomData.getRowList()
        ]

        bondData = block.getObj('chem_comp_bond')
        if bondData is not None:
            atom1Col = bondData.getAttributeIndex('atom_id_1')
            atom2Col = bondData.getAttributeIndex('atom_id_2')
            orderCol = bondData.getAttributeIndex('value_order')
            aromaticCol = bondData.getAttributeIndex('pdbx_aromatic_flag')
            bonds = [
                CCDBondDefinition(
                    atom1=row[atom1Col],
                    atom2=row[atom2Col],
                    order=row[orderCol],
                    aromatic=row[aromaticCol] == 'Y',
                ) for row in bondData.getRowList()
            ]
        else:
            bonds = []

        return cls(residueName=residueName, atoms=atoms, bonds=bonds)

def _overlayPoints(points1, points2):
    if len(points1) == 0:
        return mm.Vec3(0, 0, 0), np.identity(3), mm.Vec3(0, 0, 0)
    if len(points1) == 1:
        return points1[0], np.identity(3), -1 * points2[0]

    # Compute centroids.

    center1 = unit.sum(points1)/float(len(points1))
    center2 = unit.sum(points2)/float(len(points2))

    # Compute R matrix.

    R = np.zeros((3, 3))
    for p1, p2 in zip(points1, points2):
        x = p1-center1
        y = p2-center2
        for i in range(3):
            for j in range(3):
                R[i][j] += y[i]*x[j]

    # Use an SVD to compute the rotation matrix.

    (u, s, v) = lin.svd(R)
    return -1 * center2, np.dot(u, v).transpose(), center1

def _dihedralRotation(points, angle):
    """Given four points that form a dihedral, compute the matrix that rotates the last point around the axis to
    produce the desired dihedral angle."""
    points = [p.value_in_unit(unit.nano) for p in points]
    v0 = points[0]-points[1]
    v1 = points[2]-points[1]
    v2 = points[2]-points[3]
    cp0 = np.cross(v0, v1)
    cp1 = np.cross(v1, v2)
    axis = v1/unit.norm(v1)
    currentAngle = np.arctan2(np.dot(np.cross(cp0, cp1), axis), np.dot(cp0, cp1))
    ct = np.cos(angle-currentAngle)
    st = np.sin(angle-currentAngle)
    return np.array([
        [axis[0]*axis[0]*(1-ct) + ct,            axis[0]*axis[1]*(1-ct) - axis[2]*st,    axis[0]*axis[2]*(1-ct) + axis[1]*st],
        [axis[1]*axis[0]*(1-ct) + axis[2]*st,    axis[1]*axis[1]*(1-ct) + ct,            axis[1]*axis[2]*(1-ct) - axis[0]*st],
        [axis[2]*axis[0]*(1-ct) - axis[1]*st,    axis[2]*axis[1]*(1-ct) + axis[0]*st,    axis[2]*axis[2]*(1-ct) + ct        ]
    ])

def _findUnoccupiedDirection(point, positions):
    """Given a point in space and a list of atom positions, find the direction in which the local density of atoms is lowest."""

    point = point.value_in_unit(unit.nano)
    direction = mm.Vec3(0, 0, 0)
    for pos in positions.value_in_unit(unit.nano):
        delta = pos-point
        distance = unit.norm(delta)
        if distance > 0.1:
            distance2 = distance*distance
            direction -= delta/(distance2*distance2)
    direction /= unit.norm(direction)
    return direction

class PDBFixer(object):
    def __init__(self, file_data: StringIO):
        self._initializeFromPDBx(file_data)

        atoms = list(self.topology.atoms())
        if len(atoms) == 0:
            raise Exception("Structure contains no atoms.")

        self._ccdCache = {}

        self.templates = {}
        templatesPath = os.path.join(os.path.dirname(__file__), 'templates')
        for file in os.listdir(templatesPath):
            templatePdb = app.PDBFile(os.path.join(templatesPath, file))
            name = next(templatePdb.topology.residues()).name
            self.templates[name] = Template(templatePdb.topology, templatePdb.positions)

    def _initializeFromPDBx(self, file):
        pdbx = app.PDBxFile(file)
        self.topology = pdbx.topology
        self.positions = pdbx.positions

        # PDBxFile doesn't record the information about sequence or modified residues, so we need to read them separately.

        file.seek(0)
        reader = PdbxReader(file)
        data = []
        reader.read(data)
        block = data[0]

        # Load the sequence data.

        sequenceData = block.getObj('entity_poly_seq')
        sequences = {}
        if sequenceData is not None:
            entityIdCol = sequenceData.getAttributeIndex('entity_id')
            residueCol = sequenceData.getAttributeIndex('mon_id')
            for row in sequenceData.getRowList():
                entityId = row[entityIdCol]
                residue = row[residueCol]
                if entityId not in sequences:
                    sequences[entityId] = []
                sequences[entityId].append(residue)

        # Sequences are stored by "entity".  There could be multiple chains that are all the same entity, so we need to
        # convert from entities to chains.

        asymData = block.getObj('struct_asym')
        self.sequences = []
        if asymData is not None:
            asymIdCol = asymData.getAttributeIndex('id')
            entityIdCol = asymData.getAttributeIndex('entity_id')
            for row in asymData.getRowList():
                asymId = row[asymIdCol]
                entityId = row[entityIdCol]
                if entityId in sequences:
                    self.sequences.append(Sequence(asymId, sequences[entityId]))

        # Load the modified residues.

        modData = block.getObj('pdbx_struct_mod_residue')
        self.modifiedResidues = []
        if modData is not None:
            asymIdCol = modData.getAttributeIndex('label_asym_id')
            resNameCol = modData.getAttributeIndex('label_comp_id')
            resNumCol = modData.getAttributeIndex('auth_seq_id')
            standardResCol = modData.getAttributeIndex('parent_comp_id')
            if -1 not in (asymIdCol, resNameCol, resNumCol, standardResCol):
                for row in modData.getRowList():
                    self.modifiedResidues.append(ModifiedResidue(row[asymIdCol], int(row[resNumCol]), row[resNameCol], row[standardResCol]))


    def _downloadCCDDefinition(self, residueName: str) -> Optional[CCDResidueDefinition]:
        residueName = residueName.upper()

        if residueName in self._ccdCache:
            return self._ccdCache[residueName]

        try:
            file = urlopen(f'https://files.rcsb.org/ligands/download/{residueName}.cif')
            contents = file.read().decode('utf-8')
            file.close()
        except:
            # None represents that the residue has been looked up and could not
            # be found. This is distinct from an entry simply not being present
            # in the cache.
            self._ccdCache[residueName] = None
            return None

        reader = PdbxReader(StringIO(contents))
        ccdDefinition = CCDResidueDefinition.fromReader(reader)
        self._ccdCache[residueName] = ccdDefinition
        return ccdDefinition

    def _getTemplate(self, name):
        """Return the template with a name.  If none has been registered, this will return None."""
        if name in self.templates:
            return self.templates[name]
        return None

    def _addAtomsToTopology(self, heavyAtomsOnly, omitUnknownMolecules):
        """Create a new Topology in which missing atoms have been added.

        Parameters
        ----------
        heavyAtomsOnly : bool
            If True, only heavy atoms will be added to the topology.
        omitUnknownMolecules : bool
            If True, unknown molecules will be omitted from the topology.

        Returns
        -------
        newTopology : openmm.app.Topology
            A new Topology object containing atoms from the old.
        newPositions : list of openmm.unit.Quantity with units compatible with nanometers
            Atom positions for the new Topology object.
        newAtoms : openmm.app.Topology.Atom
            New atom objects.
        existingAtomMap : dict
            Mapping from old atoms to new atoms.

        """

        newTopology = app.Topology()
        newPositions = []*unit.nanometer
        newAtoms = []
        existingAtomMap = {}
        addedAtomMap = {}
        addedHeterogenBonds = []
        residueCenters = [self._computeResidueCenter(res).value_in_unit(unit.nanometers) for res in self.topology.residues()]*unit.nanometers
        for chain in self.topology.chains():
            if omitUnknownMolecules and all(self._getTemplate(residue.name) is None for residue in chain.residues()):
                continue
            chainResidues = list(chain.residues())
            newChain = newTopology.addChain(chain.id)
            for indexInChain, residue in enumerate(chain.residues()):

                # Insert missing residues here.

                if (chain.index, indexInChain) in self.missingResidues:
                    insertHere = self.missingResidues[(chain.index, indexInChain)]
                    endPosition = self._computeResidueCenter(residue)
                    if indexInChain > 0:
                        startPosition = self._computeResidueCenter(chainResidues[indexInChain-1])
                        loopDirection = _findUnoccupiedDirection((startPosition+endPosition)/2, residueCenters)
                    else:
                        outward = _findUnoccupiedDirection(endPosition, residueCenters)*unit.nanometers
                        norm = unit.norm(outward)
                        if norm > 0*unit.nanometer:
                            outward *= len(insertHere)*0.5*unit.nanometer/norm
                        startPosition = endPosition+outward
                        loopDirection = None
                    firstIndex = int(residue.id)-len(insertHere)
                    self._addMissingResiduesToChain(newChain, insertHere, startPosition, endPosition, loopDirection, residue, newAtoms, newPositions, firstIndex)

                # Create the new residue and add existing heavy atoms.

                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                for atom in residue.atoms():
                    if not heavyAtomsOnly or (atom.element is not None and atom.element != hydrogen):
                        if atom.name == 'OXT' and (chain.index, indexInChain+1) in self.missingResidues:
                            continue # Remove terminal oxygen, since we'll add more residues after this one
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        existingAtomMap[atom] = newAtom
                        newPositions.append(self.positions[atom.index])
                if residue in self.missingAtoms:

                    # Find corresponding atoms in the residue and the template.

                    template = self._getTemplate(residue.name)
                    atomPositions = dict((atom.name, self.positions[atom.index]) for atom in residue.atoms())
                    points1 = []
                    points2 = []
                    for atom in template.topology.atoms():
                        if atom.name in atomPositions:
                            points1.append(atomPositions[atom.name].value_in_unit(unit.nanometer))
                            points2.append(template.positions[atom.index].value_in_unit(unit.nanometer))

                    # Compute the optimal transform to overlay them.

                    (translate2, rotate, translate1) = _overlayPoints(points1, points2)

                    # Add the missing atoms.

                    addedAtomMap[residue] = {}
                    for atom in self.missingAtoms[residue]:
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms.append(newAtom)
                        addedAtomMap[residue][atom] = newAtom
                        templatePosition = template.positions[atom.index].value_in_unit(unit.nanometer)
                        newPositions.append((mm.Vec3(*np.dot(rotate, templatePosition+translate2))+translate1)*unit.nanometer)

                    if residue.name not in app.Topology._standardBonds:

                        # This is a heterogen.  Make sure bonds will get added for any new atoms.

                        addedAtomNames = set(atom.name for atom in addedAtomMap[residue])
                        newResidueAtoms = {atom.name: atom for atom in newResidue.atoms()}
                        for atom1, atom2 in template.topology.bonds():
                            if atom1.name in addedAtomNames or atom2.name in addedAtomNames:
                                if atom1.name in newResidueAtoms and atom2.name in newResidueAtoms:
                                    addedHeterogenBonds.append((newResidueAtoms[atom1.name], newResidueAtoms[atom2.name]))

                if residue in self.missingTerminals:
                    terminalsToAdd = self.missingTerminals[residue]
                else:
                    terminalsToAdd = None

                # If this is the end of the chain, add any missing residues that come after it.

                if residue == chainResidues[-1] and (chain.index, indexInChain+1) in self.missingResidues:
                    insertHere = self.missingResidues[(chain.index, indexInChain+1)]
                    if len(insertHere) > 0:
                        startPosition = self._computeResidueCenter(residue)
                        outward = _findUnoccupiedDirection(startPosition, residueCenters)*unit.nanometers
                        norm = unit.norm(outward)
                        if norm > 0*unit.nanometer:
                            outward *= len(insertHere)*0.5*unit.nanometer/norm
                        endPosition = startPosition+outward
                        firstIndex = int(residue.id)+1
                        self._addMissingResiduesToChain(newChain, insertHere, startPosition, endPosition, None, residue, newAtoms, newPositions, firstIndex)
                        newResidue = list(newChain.residues())[-1]
                        if newResidue.name in proteinResidues:
                            terminalsToAdd = ['OXT']
                        else:
                            terminalsToAdd = None

                # If a terminal OXT is missing, add it.

                if terminalsToAdd is not None:
                    atomPositions = dict((atom.name, newPositions[atom.index].value_in_unit(unit.nanometer)) for atom in newResidue.atoms())
                    if 'OXT' in terminalsToAdd:
                        newAtom = newTopology.addAtom('OXT', oxygen, newResidue)
                        newAtoms.append(newAtom)
                        d_ca_o = atomPositions['O']-atomPositions['CA']
                        d_ca_c = atomPositions['C']-atomPositions['CA']
                        d_ca_c /= unit.sqrt(unit.dot(d_ca_c, d_ca_c))
                        v = d_ca_o - d_ca_c*unit.dot(d_ca_c, d_ca_o)
                        newPositions.append((atomPositions['O']+2*v)*unit.nanometer)
        newTopology.setUnitCellDimensions(self.topology.getUnitCellDimensions())
        newTopology.createStandardBonds()
        newTopology.createDisulfideBonds(newPositions)

        # Add the existing bonds between atoms in heterogens.

        for a1,a2 in self.topology.bonds():
            if a1 in existingAtomMap and a2 in existingAtomMap and (a1.residue.name not in app.Topology._standardBonds or a2.residue.name not in app.Topology._standardBonds):
                newTopology.addBond(existingAtomMap[a1], existingAtomMap[a2])

        # Add any new bonds within heterogens.

        bonds = set((atom1.index, atom2.index) for atom1, atom2 in newTopology.bonds())
        for atom1, atom2 in addedHeterogenBonds:
            if (atom1.index, atom2.index) not in bonds and (atom2.index, atom1.index) not in bonds:
                newTopology.addBond(atom1, atom2)
                bonds.add((atom1.index, atom2.index))

        # Return the results.

        return newTopology, newPositions, newAtoms, existingAtomMap

    def _computeResidueCenter(self, residue):
        """Compute the centroid of a residue."""
        return unit.sum([self.positions[atom.index] for atom in residue.atoms()])/len(list(residue.atoms()))

    def _addMissingResiduesToChain(self, chain, residueNames, startPosition, endPosition, loopDirection, orientTo, newAtoms, newPositions, firstIndex):
        """Add a series of residues to a chain."""
        orientToPositions = dict((atom.name, self.positions[atom.index]) for atom in orientTo.atoms())
        if loopDirection is None:
            loopDirection = mm.Vec3(0, 0, 0)

        # We'll add the residues in an arc connecting the endpoints.  Figure out the height of that arc.

        length = unit.norm(endPosition-startPosition)
        numResidues = len(residueNames)
        if length > numResidues*0.3*unit.nano:
            loopHeight = 0*unit.nano
        else:
            loopHeight = (numResidues*0.3*unit.nano-length)/2

        # Add the residues.

        try:
            prevResidue = list(chain.residues())[-1]
        except:
            prevResidue = None
        for i, residueName in enumerate(residueNames):
            template = self._getTemplate(residueName)

            # Find a translation that best matches the adjacent residue.

            points1 = []
            points2 = []
            for atom in template.topology.atoms():
                if atom.name in orientToPositions:
                    points1.append(orientToPositions[atom.name].value_in_unit(unit.nano))
                    points2.append(template.positions[atom.index].value_in_unit(unit.nano))
            (translate2, rotate, translate1) = _overlayPoints(points1, points2)

            # Create the new residue.

            newResidue = chain.topology.addResidue(residueName, chain, "%d" % ((firstIndex+i)%10000))
            fraction = (i+1.0)/(numResidues+1.0)
            translate = startPosition + (endPosition-startPosition)*fraction + loopHeight*math.sin(fraction*math.pi)*loopDirection
            templateAtoms = list(template.topology.atoms())
            if newResidue == next(chain.residues()):
                templateAtoms = [atom for atom in templateAtoms if atom.name not in ('P', 'OP1', 'OP2')]
            for atom in templateAtoms:
                newAtom = chain.topology.addAtom(atom.name, atom.element, newResidue)
                newAtoms.append(newAtom)
                templatePosition = template.positions[atom.index].value_in_unit(unit.nano)
                newPositions.append(mm.Vec3(*np.dot(rotate, templatePosition))*unit.nano+translate)
            if prevResidue is not None:
                atoms1 = {atom.name: atom for atom in prevResidue.atoms()}
                atoms2 = {atom.name: atom for atom in newResidue.atoms()}
                if 'CA' in atoms1 and 'C' in atoms1 and 'N' in atoms2 and 'CA' in atoms2:

                    # We're adding a peptide bond between this residue and the previous one.  Rotate it to try to
                    # put the peptide bond into the trans conformation.

                    atoms = (atoms1['CA'], atoms1['C'], atoms2['N'], atoms2['CA'])
                    points = [newPositions[a.index] for a in atoms]
                    rotation = _dihedralRotation(points, np.pi)
                    for atom in newResidue.atoms():
                        d = (newPositions[atom.index]-points[2]).value_in_unit(unit.nano)
                        newPositions[atom.index] = mm.Vec3(*np.dot(rotation, d))*unit.nano + points[2]

            prevResidue = newResidue

    def findMissingResidues(self):
        chains = [c for c in self.topology.chains() if len(list(c.residues())) > 0]
        chainWithGaps = {}

        # Find the sequence of each chain, with gaps for missing residues.

        for chain in chains:
            residues = list(chain.residues())
            ids = [int(r.id) for r in residues]
            for i, res in enumerate(residues):
                if res.insertionCode not in ('', ' '):
                    for j in range(i, len(residues)):
                        ids[j] += 1
            minResidue = min(ids)
            maxResidue = max(ids)
            chainWithGaps[chain] = [None]*(maxResidue-minResidue+1)
            for r, id in zip(residues, ids):
                chainWithGaps[chain][id-minResidue] = r.name

        # Try to find the chain that matches each sequence.

        chainSequence = {}
        chainOffset = {}
        for sequence in self.sequences:
            for chain in chains:
                if chain.id != sequence.chainId:
                    continue
                if chain in chainSequence:
                    continue
                for offset in range(len(sequence.residues)-len(chainWithGaps[chain])+1):
                    if all(a == b or b == None for a,b in zip(sequence.residues[offset:], chainWithGaps[chain])):
                        chainSequence[chain] = sequence
                        chainOffset[chain] = offset
                        break
                if chain in chainSequence:
                    break

        # Now build the list of residues to add.

        self.missingResidues = {}
        for chain in self.topology.chains():
            if chain in chainSequence:
                offset = chainOffset[chain]
                sequence = chainSequence[chain].residues
                gappedSequence = chainWithGaps[chain]
                index = 0
                for i in range(len(sequence)):
                    if i < offset or i >= len(gappedSequence)+offset or gappedSequence[i-offset] is None:
                        key = (chain.index, index)
                        if key not in self.missingResidues:
                            self.missingResidues[key] = []
                        residueName = sequence[i]
                        if residueName in substitutions:
                            residueName = substitutions[sequence[i]]
                        self.missingResidues[key].append(residueName)
                    else:
                        index += 1

    def findMissingAtoms(self):
        missingAtoms = {}
        missingTerminals = {}

        # Determine which atoms have an external bond to another residue.

        hasExternal = defaultdict(bool)
        for atom1, atom2 in self.topology.bonds():
            if atom1.residue != atom2.residue:
                hasExternal[(atom1.residue, atom1.name)] = True
                hasExternal[(atom2.residue, atom2.name)] = True
        for chain in self.topology.chains():
            chainResidues = list(chain.residues())
            for residue in chain.residues():
                atomNames = [atom.name for atom in residue.atoms()]
                if all(name in atomNames for name in ['C', 'O', 'CA']):
                    # We'll be adding peptide bonds.
                    if residue != chainResidues[0]:
                        hasExternal[(residue, 'N')] = True
                    if residue != chainResidues[-1]:
                        hasExternal[(residue, 'C')] = True

        # Loop over residues.

        for chain in self.topology.chains():
            nucleic = any(res.name in dnaResidues or res.name in rnaResidues for res in chain.residues())
            chainResidues = list(chain.residues())
            for residue in chain.residues():
                template = self._getTemplate(residue.name)
                if template is not None:
                    # If an atom is marked as terminal only, and if it is bonded to any atom that has an external bond
                    # to another residue, we need to omit that atom and any other terminal-only atom bonded to it.

                    bondedTo = defaultdict(set)
                    for atom1, atom2 in template.topology.bonds():
                        bondedTo[atom1].add(atom2)
                        bondedTo[atom2].add(atom1)
                    skip = set()
                    for atom, terminal in zip(template.topology.atoms(), template.terminal):
                        if terminal:
                            for atom2 in bondedTo[atom]:
                                if hasExternal[(residue, atom2.name)]:
                                    skip.add(atom)
                    for atom, terminal in zip(template.topology.atoms(), template.terminal):
                        if terminal:
                            for atom2 in bondedTo[atom]:
                                if atom2 in skip:
                                    skip.add(atom)
                    atomNames = set(atom.name for atom in residue.atoms())
                    templateAtoms = [atom for atom in template.topology.atoms() if atom not in skip]
                    if nucleic and residue == chainResidues[0] and (chain.index, 0) not in self.missingResidues:
                        templateAtoms = [atom for atom in templateAtoms if atom.name not in ('P', 'OP1', 'OP2')]

                    # Add atoms from the template that are missing.

                    missing = []
                    for atom in templateAtoms:
                        if atom.name not in atomNames and atom.element != app.element.hydrogen:
                            missing.append(atom)
                    if len(missing) > 0:
                        missingAtoms[residue] = missing

                    # Add missing terminal atoms.

                    terminals = []
                    if residue == chainResidues[-1] and (chain.index, len(chainResidues)) not in self.missingResidues:
                        templateNames = set(atom.name for atom in template.topology.atoms())
                        if 'OXT' not in atomNames and all(name in templateNames for name in ['C', 'O', 'CA']):
                            terminals.append('OXT')
                        if len(terminals) > 0:
                            missingTerminals[residue] = terminals
        self.missingAtoms = missingAtoms
        self.missingTerminals = missingTerminals

    def addMissingAtoms(self, seed=None):
        # Create a Topology that 1) adds missing atoms, 2) removes all hydrogens, and 3) removes unknown molecules.

        (newTopology, newPositions, newAtoms, existingAtomMap) = self._addAtomsToTopology(True, True)
        if len(newAtoms) == 0:

            # No atoms were added, but new bonds might have been created.

            newBonds = set(newTopology.bonds())
            for atom1, atom2 in self.topology.bonds():
                if atom1 in existingAtomMap and atom2 in existingAtomMap:
                    a1 = existingAtomMap[atom1]
                    a2 = existingAtomMap[atom2]
                    if (a1, a2) in newBonds:
                        newBonds.remove((a1, a2))
                    elif (a2, a1) in newBonds:
                        newBonds.remove((a2, a1))

            # Add the new bonds to the original Topology.

            inverseAtomMap = dict((y,x) for (x,y) in existingAtomMap.items())
            for atom1, atom2 in newBonds:
                self.topology.addBond(inverseAtomMap[atom1], inverseAtomMap[atom2])
        else:

            # Create a System for energy minimizing it.

            forcefield = self._createForceField(newTopology, False)
            system = forcefield.createSystem(newTopology)

            # Set any previously existing atoms to be massless, they so won't move.

            for atom in existingAtomMap.values():
                system.setParticleMass(atom.index, 0.0)

            # If any heavy atoms were omitted, add them back to avoid steric clashes.

            nonbonded = [f for f in system.getForces() if isinstance(f, mm.CustomNonbondedForce)][0]
            for atom in self.topology.atoms():
                if atom.element not in (None, hydrogen) and atom not in existingAtomMap:
                    system.addParticle(0.0)
                    nonbonded.addParticle([])
                    newPositions.append(self.positions[atom.index])

            # For efficiency, only compute interactions that involve a new atom.

            nonbonded.addInteractionGroup([atom.index for atom in newAtoms], range(system.getNumParticles()))

            # Do an energy minimization.

            integrator = mm.LangevinIntegrator(300*unit.kelvin, 10/unit.picosecond, 5*unit.femtosecond)
            if seed is not None:
                integrator.setRandomNumberSeed(seed)
            context = mm.Context(system, integrator)
            context.setPositions(newPositions)
            mm.LocalEnergyMinimizer.minimize(context)
            state = context.getState(getPositions=True)
            if newTopology.getNumResidues() > 1:
                # When looking for pairs of atoms that are too close to each other, exclude pairs that
                # are in the same residue or are directly bonded to each other.

                exclusions = dict((atom, {a.index for a in atom.residue.atoms()}) for atom in newAtoms)
                for a1, a2 in newTopology.bonds():
                    if a1 in exclusions:
                        exclusions[a1].add(a2.index)
                    if a2 in exclusions:
                        exclusions[a2].add(a1.index)
                cutoff = 0.13
                nearest = self._findNearestDistance(context, newAtoms, cutoff, exclusions)
                if nearest < cutoff:

                    # Some atoms are very close together.  Run some dynamics while slowly increasing the strength of the
                    # repulsive interaction to try to improve the result.

                    for i in range(10):
                        context.setParameter('C', 0.15*(i+1))
                        integrator.step(200)
                        d = self._findNearestDistance(context, newAtoms, cutoff, exclusions)
                        if d > nearest:
                            nearest = d
                            state = context.getState(getPositions=True)
                            if nearest >= cutoff:
                                break
                    context.setState(state)
                    context.setParameter('C', 1.0)
                    mm.LocalEnergyMinimizer.minimize(context)
                    state = context.getState(getPositions=True)

            # Now create a new Topology, including all atoms from the original one and adding the missing atoms.

            (newTopology2, newPositions2, newAtoms2, existingAtomMap2) = self._addAtomsToTopology(False, False)

            # Copy over the minimized positions for the new atoms.

            for a1, a2 in zip(newAtoms, newAtoms2):
                newPositions2[a2.index] = state.getPositions()[a1.index]
            self.topology = newTopology2
            self.positions = newPositions2

    def addMissingHydrogens(self, pH=7.0, forcefield=None):
        extraDefinitions = self._downloadNonstandardDefinitions()
        variants = [self._describeVariant(res, extraDefinitions) for res in self.topology.residues()]
        modeller = app.Modeller(self.topology, self.positions)
        modeller.addHydrogens(pH=pH, forcefield=forcefield, variants=variants)
        self.topology = modeller.topology
        self.positions = modeller.positions

    def _downloadNonstandardDefinitions(self):
        """If the file contains any nonstandard residues, download their definitions and build
        the information needed to add hydrogens to them.
        """
        app.Modeller._loadStandardHydrogenDefinitions()
        resnames = set(residue.name for residue in self.topology.residues())
        definitions = {}
        for name in resnames:
            if name not in app.Modeller._residueHydrogens:
                # Try to download the definition.
                ccdDefinition = self._downloadCCDDefinition(name)
                if ccdDefinition is None:
                    continue

                # Record the atoms and bonds.
                atoms = [(atom.atomName, atom.symbol.upper(), atom.leaving) for atom in ccdDefinition.atoms]
                bonds = [(bond.atom1, bond.atom2) for bond in ccdDefinition.bonds]
                definitions[name] = (atoms, bonds)
        return definitions

    def _describeVariant(self, residue, definitions):
        if residue.name not in app.PDBFile._standardResidues and self._getTemplate(residue.name) is not None:
            # The user has registered a template for this residue.  Use the hydrogens from it.
            template = self._getTemplate(residue.name)
            atoms = [(atom.name, atom.element.symbol.upper(), terminal) for atom, terminal in zip(template.topology.atoms(), template.terminal)]
            resAtoms = dict((atom.name, atom) for atom in residue.atoms())
            bonds = []
            for atom1, atom2 in template.topology.bonds():
                if atom1.element == app.element.hydrogen and atom2.name in resAtoms:
                    bonds.append((atom1.name, atom2.name))
                elif atom2.element == app.element.hydrogen and atom1.name in resAtoms:
                    bonds.append((atom2.name, atom1.name))
        elif residue.name in definitions:
            # We downloaded a definition.
            atoms, bonds = definitions[residue.name]
        else:
            return None

        # See if the heavy atoms are identical.

        topologyHeavy = dict((atom.name, atom) for atom in residue.atoms() if atom.element is not None and atom.element != app.element.hydrogen)
        definitionHeavy = dict((atom[0], atom) for atom in atoms if atom[1] != '' and atom[1] != 'H')
        for name in topologyHeavy:
            if name not in definitionHeavy or definitionHeavy[name][1] != topologyHeavy[name].element.symbol.upper():
                # This atom isn't present in the definition
                return None
        for name in definitionHeavy:
            if name not in topologyHeavy and not definitionHeavy[name][2]:
                # This isn't a leaving atom, and it isn't present in the topology.
                return None

        # Build the list of hydrogens.

        variant = []
        definitionAtoms = dict((atom[0], atom) for atom in atoms)
        topologyBonds = list(residue.bonds())
        for name1, name2 in bonds:
            if definitionAtoms[name1][1] == 'H':
                h, parent = name1, name2
            elif definitionAtoms[name2][1] == 'H':
                h, parent = name2, name1
            else:
                continue
            if definitionAtoms[h][2]:
                # The hydrogen is marked as a leaving atom.  Omit it if the parent is not present,
                # or if the parent has an external bond.
                if parent not in topologyHeavy:
                    continue
                parentAtom = topologyHeavy[parent]
                exclude = False
                for atom1, atom2 in topologyBonds:
                    if atom1 == parentAtom and atom2.residue != residue:
                        exclude = True
                        break
                    if atom2 == parentAtom and atom1.residue != residue:
                        exclude = True
                        break
                if exclude:
                    continue
            variant.append((h, parent))
        return variant

    def _downloadFormalCharges(self, resName: str, includeLeavingAtoms: bool = True) -> dict[str, int]:
        # Try to download the definition.
        ccdDefinition = self._downloadCCDDefinition(resName.upper())
        if ccdDefinition is None:
            return {}

        # Record the formal charges.
        return {
            atom.atomName: atom.charge
            for atom in ccdDefinition.atoms
            if includeLeavingAtoms or not atom.leaving
        }

    def _createForceField(self, newTopology, water):
        if water:
            forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
            nonbonded = [f for f in forcefield._forces if isinstance(f, NonbondedGenerator)][0]
            radii = {'H':0.198, 'Li':0.203, 'C':0.340, 'N':0.325, 'O':0.299, 'F':0.312, 'Na':0.333, 'Mg':0.141,
                     'P':0.374, 'S':0.356, 'Cl':0.347, 'K':0.474, 'Br':0.396, 'Rb':0.527, 'I':0.419, 'Cs':0.605}
        else:
            forcefield = app.ForceField(os.path.join(os.path.dirname(__file__), 'soft.xml'))

        # The Topology may contain residues for which the ForceField does not have a template.
        # If so, we need to create new templates for them.

        atomTypes = {}
        bondedToAtom = []
        for atom in newTopology.atoms():
            bondedToAtom.append(set())
        for atom1, atom2 in newTopology.bonds():
            bondedToAtom[atom1.index].add(atom2.index)
            bondedToAtom[atom2.index].add(atom1.index)
        for residue in newTopology.residues():

            # Make sure the ForceField has a template for this residue.

            signature = app.forcefield._createResidueSignature([atom.element for atom in residue.atoms()])
            if signature in forcefield._templateSignatures:
                if any(matchResidue(residue, t, bondedToAtom) is not None for t in forcefield._templateSignatures[signature]):
                    continue

            # Create a new template.

            resName = "extra_"+residue.name
            template = app.ForceField._TemplateData(resName)
            forcefield._templates[resName] = template
            indexInResidue = {}
            # If we can't find formal charges in the CCD, make everything uncharged
            formalCharges = defaultdict(int)
            # See if we can get formal charges from the CCD
            if water:
                # The formal charges in the CCD can only be relied on if the
                # residue has all and only the same atoms, with the caveat that
                # leaving atoms are optional
                downloadedFormalCharges = self._downloadFormalCharges(residue.name)
                essentialAtoms = set(
                    self._downloadFormalCharges(residue.name, includeLeavingAtoms=False)
                )
                atomsInResidue = {atom.name for atom in residue.atoms()}
                if (
                    atomsInResidue.issuperset(essentialAtoms)
                    and atomsInResidue.issubset(downloadedFormalCharges)
                ):
                    # We got formal charges and the atom names match, so we can use them
                    formalCharges = downloadedFormalCharges
            for atom in residue.atoms():
                element = atom.element
                formalCharge = formalCharges.get(atom.name, 0)
                typeName = 'extra_'+element.symbol+'_'+str(formalCharge)
                if (element, formalCharge) not in atomTypes:
                    atomTypes[(element, formalCharge)] = app.ForceField._AtomType(typeName, '', 0.0, element)
                    forcefield._atomTypes[typeName] = atomTypes[(element, formalCharge)]
                    if water:
                        # Select a reasonable vdW radius for this atom type.

                        if element.symbol in radii:
                            sigma = radii[element.symbol]
                        else:
                            sigma = 0.5
                        nonbonded.registerAtom({
                            'type':typeName,
                            'charge':str(formalCharge),
                            'sigma':str(sigma),
                            'epsilon':'0'
                        })
                indexInResidue[atom.index] = len(template.atoms)
                template.atoms.append(app.ForceField._TemplateAtomData(atom.name, typeName, element))
            for atom in residue.atoms():
                for bondedTo in bondedToAtom[atom.index]:
                    if bondedTo in indexInResidue:
                        b = (indexInResidue[atom.index], indexInResidue[bondedTo])
                        if b[0] < b[1]:
                            template.bonds.append(b)
                            template.atoms[b[0]].bondedTo.append(b[1])
                            template.atoms[b[1]].bondedTo.append(b[0])
                    else:
                        b = indexInResidue[atom.index]
                        template.externalBonds.append(b)
                        template.atoms[b].externalBonds += 1
            if signature in forcefield._templateSignatures:
                forcefield._templateSignatures[signature].append(template)
            else:
                forcefield._templateSignatures[signature] = [template]
        return forcefield

    def _findNearestDistance(self, context, newAtoms, cutoff, exclusions):
        """Given a set of newly added atoms, find the closest distance between one of those atoms and another atom."""

        positions = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        boxSize = np.max(positions, axis=0)-np.min(positions, axis=0)
        boxVectors = [(boxSize[0], 0, 0), (0, boxSize[1], 0), (0, 0, boxSize[2])]
        cells = app.modeller._CellList(positions, cutoff, boxVectors, False)
        nearest_squared = sys.float_info.max
        for atom in newAtoms:
            excluded = exclusions[atom]
            for i in cells.neighbors(positions[atom.index]):
                if i not in excluded:
                    p = positions[atom.index]-positions[i]
                    dist_squared = np.dot(p, p)
                    if dist_squared < nearest_squared:
                        nearest_squared = dist_squared
        return np.sqrt(nearest_squared)
