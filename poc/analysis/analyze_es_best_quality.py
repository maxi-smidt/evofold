import csv

import matplotlib.pyplot as plt

from dataclasses import replace
from functools import partial
from typing import List, Tuple

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

def plot_plus_comma_over_rmsd(plus_file: str, comma_file: str, show_rmsd=False):
    c_fitness, c_rmsd = get_fitness_and_rmsd(comma_file)
    p_fitness, p_rmsd = get_fitness_and_rmsd(plus_file)
    assert len(c_fitness) == len(p_fitness)

    generations = range(len(c_fitness))

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Remove symbols on the data points, just plot lines
    ax1.plot(generations, c_fitness, 'b-', label='Fitness (µ,λ)')  # 'b-' for blue solid line
    ax1.plot(generations, p_fitness, 'c-', label='Fitness (µ+λ)')  # 'c-' for cyan solid line
    ax1.set_xlabel('Generation')
    ax1.set_ylabel('Fitness [kJ/mol]')
    ax1.tick_params(axis='y')

    lines_1, labels_1 = ax1.get_legend_handles_labels()

    if show_rmsd:
        ax2 = ax1.twinx()
        # RMSD lines with dashed style (remove symbols)
        ax2.plot(generations, c_rmsd, 'r--', label='RMSD (µ,λ)')  # 'r--' for red dashed line
        ax2.plot(generations, p_rmsd, 'm--', label='RMSD (µ+λ)')  # 'm--' for magenta dashed line
        ax2.set_ylabel('RMSD to AlphaFold')
        ax2.tick_params(axis='y')
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')
    else:
        ax1.legend(lines_1, labels_1, loc='upper right')


    plt.title('Evolution of Best Fitness and RMSD Across Generations')
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
        (replace(base_params, plus_selection=True, force_field='charmm'), 'plus-charmm'),
    ]

    measure_and_save(params, sequence)

    # plot_plus_comma_over_rmsd('plus-amber', 'comma-amber', True)
    # plot_plus_comma_over_rmsd('plus-amber', 'comma-amber', False)
    # plot_plus_comma_over_rmsd('plus-charmm', 'comma-charmm', True)
    # plot_plus_comma_over_rmsd('plus-charmm', 'comma-charmm', False)


if __name__ == "__main__":
    main()