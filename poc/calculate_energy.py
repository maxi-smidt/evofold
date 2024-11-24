from openmm import app, VerletIntegrator
from openmm import unit
from openmm.app import PDBFile
from openmm import Platform

# Load the protein structure from a PDB file
pdb = PDBFile('test.cif')

# Apply a force field (use appropriate force field files)
forcefield = app.ForceField('amber99sb.xml', 'amber99_obc.xml')

# Create the system using the PDB topology and the force field
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)

# Set up the simulation context
platform = Platform.getPlatformByName('CPU')  # Or 'CUDA' if you want to use a GPU
integrator = VerletIntegrator(0.001)  # A simple integrator, not needed for just energy calculation
simulation = app.Simulation(pdb.topology, system, integrator, platform)

# Set the atom positions (if you've modified them)
# If you haven't modified the structure, this is not necessary
positions = pdb.getPositions()
simulation.context.setPositions(positions)

# Calculate the potential energy of the system
state = simulation.context.getState(getEnergy=True)
energy = state.getPotentialEnergy()

# Print the energy in kilocalories per mole (kcal/mol)
print(f"Potential Energy: {energy.value_in_unit(unit.kilocalories_per_mole)} kcal/mol")
