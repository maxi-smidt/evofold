import io

from openmm import VerletIntegrator
from openmm.app import ForceField, NoCutoff, Simulation
from openmm.unit import pico

from poc.fixer.pdbfixer import PDBFixer

with open("test_new.cif", "r") as f:
    data = f.read()

cif_file = io.StringIO(data)
cif_file.seek(0)
fixer = PDBFixer(cif_file)

fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens()

forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

system = forcefield.createSystem(fixer.topology, nonbondedMethod=NoCutoff)

integrator = VerletIntegrator(0.001*pico.factor)
simulation = Simulation(fixer.topology, system, integrator)
simulation.context.setPositions(fixer.positions)
state = simulation.context.getState(getEnergy=True)

print(f"Potential Energy:{state.getPotentialEnergy()}")