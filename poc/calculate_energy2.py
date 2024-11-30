from openmm import LangevinIntegrator
from openmm.app import ForceField, Simulation, PME, HBonds
from openmm.unit import kelvin, pico, nano
from pdbfixer import PDBFixer

def fix_pdb(pdb_id):
    path = os.getcwd()
    if len(pdb_id) != 4:
        print("Creating PDBFixer...")
        fixer = PDBFixer(pdb_id)
        print("Finding missing residues...")
        fixer.findMissingResidues()

        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                print("ok")
                del fixer.missingResidues[key]

        print("Finding nonstandard residues...")
        fixer.findNonstandardResidues()
        print("Replacing nonstandard residues...")
        fixer.replaceNonstandardResidues()
        print("Removing heterogens...")
        fixer.removeHeterogens(keepWater=True)

        print("Finding missing atoms...")
        fixer.findMissingAtoms()
        print("Adding missing atoms...")
        fixer.addMissingAtoms()
        print("Adding missing hydrogens...")
        fixer.addMissingHydrogens(7)
        print("Writing PDB file...")

        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            open(os.path.join(path, "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)),
                 "w"),
            keepIds=True)
        return "%s_fixed_pH_%s.pdb" % (pdb_id.split('.')[0], 7)
fix_pdb('1stp.pdb')

# Prompt manual checking of chain ends for terminal groups
for chain in fixer.topology.chains():
    print(f"First residue in chain {chain.id}: {next(chain.residues()).name}")
    print(f"Last residue in chain {chain.id}: {list(chain.residues())[-1].name}")

# Load the force field
forcefield = ForceField('forcefields/amber99sb.xml')


system = forcefield.createSystem(fixer.topology, nonbondedMethod=PME, nonbondedCutoff=nano.factor,
                                     constraints=HBonds)

# Set up the simulation
integrator = LangevinIntegrator(300 * kelvin, 1 / pico.factor, 0.002 * pico.factor)
simulation = Simulation(fixer.topology, system, integrator)
simulation.context.setPositions(fixer.positions)

# Compute the potential energy
state = simulation.context.getState(getEnergy=True)
potential_energy = state.getPotentialEnergy()
print(f"Potential Energy: {potential_energy}")
