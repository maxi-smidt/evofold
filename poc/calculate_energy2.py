from simtk.openmm import app
from simtk.openmm import openmm
from simtk import unit

# Example: Define a molecule with atoms and positions
app.PDBxFile
pdb = app.PDBFile('example.pdb')  # Load your structure (you can generate a PDB file)
forcefield = app.ForceField('amber99sb.xml')  # Specify force field
system = forcefield.createSystem(pdb.topology)

# Define integrator and simulation context
integrator = openmm.LangevinIntegrator(
    300*unit.kelvin,       # Temperature
    1/unit.picoseconds,    # Friction coefficient
    0.002*unit.picoseconds # Time step
)
simulation = app.Simulation(pdb.topology, system, integrator)

# Set initial positions
simulation.context.setPositions(pdb.positions)

# Calculate potential energy
state = simulation.context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy()
print("Potential Energy:", potential_energy)
