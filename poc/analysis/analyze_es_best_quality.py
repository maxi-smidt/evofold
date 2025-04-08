import csv

import matplotlib.pyplot as plt

from dataclasses import replace
from functools import partial
from typing import List, Tuple

from openmm import VerletIntegrator
from openmm.app import ForceField, PDBxFile, Simulation
from openmm.app.forcefield import NoCutoff
from openmm.unit import kilojoule_per_mole, pico

from backend.algorithms.es import ES
from backend.algorithms.params.es_params import ESParams
from backend.algorithms.self_adaptive_es import SelfAdaptiveES
from backend.structure.protein import Protein


def handle_callback(generation: int, p: Protein, sigma: float, _: bool, best_proteins: List[Tuple[int, Protein, float]]):
    print('generation: ', generation)
    best_proteins.append((generation, p, sigma))

def measure_es(es: ES, sequence: str) -> List[Tuple[int, Protein, float]]:
    best_qualities: List[Tuple[int, Protein, float]] = []
    f = partial(handle_callback, best_proteins=best_qualities)
    es.run(sequence, f)
    return best_qualities

def measure_and_save(params: List, sequence: str):
    for p, fn in params:
        best_qualities = measure_es(es=SelfAdaptiveES(p), sequence=sequence)
        with open(f'best_qual/{fn}.csv', 'w') as file:
            file.write(f"seqeunce: {sequence}\n")
            file.write(f"generation,fitness,sigma,cif\n")
            for g, pr, sigma in best_qualities:
                file.write(f'{g},{pr.fitness},{sigma},"{pr.cif_str}"\n')


def get_fitness_and_rmsd(fn: str) -> Tuple[List[float], List[float]]:
    fitness = []
    rmsd = []
    with open(f'best_qual/{fn}.csv', 'r') as file:
        reader = csv.reader(file)
        next(reader, None) # skip sequence
        next(reader, None) # skip header
        print(fn, end='')
        for i, row in enumerate(reader):
            fitness.append(float(row[1]))
            rmsd.append(float(row[4]))
    assert len(fitness) == len(rmsd)
    return fitness, rmsd

def get_af_fitness(fn: str, ff: str):
    ff = ForceField('amber14-all.xml') if 'amber' in ff else ForceField('charmm36.xml')
    pdbx = PDBxFile(fn)
    system = ff.createSystem(pdbx.topology, nonbondedMethod=NoCutoff)
    integrator = VerletIntegrator(0.001 * pico.factor)
    simulation = Simulation(pdbx.topology, system, integrator)
    simulation.context.setPositions(pdbx.positions)
    state = simulation.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)


def plot_plus_comma_over_rmsd(file: str, title: str, show_rmsd=False, show_af_fitness=False):
    fitness, rmsd = get_fitness_and_rmsd(file)

    generations = range(len(fitness))

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(generations, fitness, 'b-')
    ax1.set_xlabel('Generation')
    ax1.set_ylabel('Fitness [kJ/mol]')
    ax1.tick_params(axis='y')

    lines_1, labels_1 = ax1.get_legend_handles_labels()

    if show_rmsd:
        ax2 = ax1.twinx()
        ax2.plot(generations, rmsd, 'r--')
        ax2.set_ylabel('RMSD to AlphaFold')
        ax2.tick_params(axis='y')
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')
    else:
        ax1.legend(lines_1, labels_1, loc='upper right')

    if show_af_fitness:
        f = get_af_fitness('best_qual/fold_peptide_v2_model_0_h.cif', file)
        plt.text(0.95, 0.95, f'AF fitness: {f:.2f} [kJ/mol]', ha='right', va='top', transform=plt.gca().transAxes,
                 fontsize=12, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))


    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def main():
    sequence = 'PEPTIDE'

    base_params = ESParams(generations=100,
                           population_size=100,
                           children_size=600,
                           plus_selection=True,
                           force_field='amber',
                           sigma=36,
                           premature_termination=None)

    params = [
        # (replace(base_params, plus_selection=False, force_field='amber'), 'comma-amber'),
        # (replace(base_params, plus_selection=True, force_field='amber'), 'plus-amber'),
        # (replace(base_params, plus_selection=False, force_field='charmm'), 'comma-charmm'),
        # (replace(base_params, plus_selection=True, force_field='charmm'), 'plus-charmm'),
    ]

    measure_and_save(params, sequence)

    plot_plus_comma_over_rmsd('plus-amber', 'AMBER - (µ+λ)', True, True)
    plot_plus_comma_over_rmsd('comma-amber', 'AMBER - (µ,λ)', True, True)
    plot_plus_comma_over_rmsd('comma-charmm', 'CHARMM - (µ,λ)', True, True)
    # plot_plus_comma_over_rmsd('plus-charmm', 'comma-charmm', False)


if __name__ == "__main__":
    main()