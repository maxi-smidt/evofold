from openmm import app
from openmm import openmm
from simtk import unit

# Example: Define a molecule with atoms and positions
pdb = app.PDBFile('1lcd.pdb')  # Load your structure (you can generate a PDB file)
force_field = app.ForceField('amber99sb.xml')  # Specify force field
system = force_field.createSystem(pdb.topology)

# Define integrator and simulation context
integrator = openmm.LangevinIntegrator(
    300*unit.kelvin,       # Temperature
    1/unit.picoseconds,    # Friction coefficient
    0.002*unit.picoseconds # Time step
)
simulation = app.Simulation(pdb.topology, system, integrator)

simulation.context.setPositions(pdb.positions)

state = simulation.context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy()
print("Potential Energy:", potential_energy)
