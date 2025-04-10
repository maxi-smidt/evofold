import csv
import matplotlib.pyplot as plt

from functools import partial
from typing import List, Tuple

from backend.algorithms.adaptive_es import AdaptiveES
from backend.algorithms.es import ES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams
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
    for es1, es2, fn in params:
        best_qualities1 = measure_es(es=es1, sequence=sequence)
        best_qualities2 = measure_es(es=es2, sequence=sequence)
        with open(f'best_qual/{fn}.csv', 'w') as file:
            file.write(f"seqeunce: {sequence}\n")
            file.write(f"generation,fitness,sigma,cif\n")
            for (g1, pr1, sigma1), (g2, pr2, sigma2) in zip(best_qualities1, best_qualities2):
                file.write(f'{g1},{pr1.fitness},{sigma1},"{pr1.cif_str}"\n')
                file.write(f'{g2},{pr2.fitness},{sigma2},"{pr2.cif_str}"\n')


def get_fitness(fn: str) -> Tuple[List[float], List[float]]:
    fitness1 = []
    fitness2 = []
    with open(f'best_qual/{fn}.csv', 'r') as file:
        reader = csv.reader(file)
        next(reader, None) # skip sequence
        next(reader, None) # skip header
        print(fn, end='')
        for i, row in enumerate(reader):
            if i % 2 == 0:
                fitness1.append(float(row[1]))
            else:
                fitness2.append(float(row[1]))
    return fitness1, fitness2


def plot(file: str, title: str, l1: str, l2: str):
    fitness1, fitness2 = get_fitness(file)
    generations = range(len(fitness1))

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(generations, fitness1, 'b-', label=l1)
    ax1.plot(generations, fitness2, 'g-', label=l2)
    ax1.set_xlabel('Generation')
    ax1.set_ylabel('Fitness [kJ/mol]')
    ax1.tick_params(axis='y')

    lines_1, labels_1 = ax1.get_legend_handles_labels()

    ax1.legend(lines_1, labels_1, loc='upper right')

    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def main():
    sequence = 'PEPTIDE'

    base_params = { 'generations':           100,
                    'population_size':       100,
                    'children_size':         600,
                    'plus_selection':        True,
                    'force_field':           'amber',
                    'sigma':                 36,
                    'premature_termination': None }

    es1 = AdaptiveES(AdaptiveESParams(**base_params))
    es2 = SelfAdaptiveES(ESParams(**base_params))

    base_params['plus_selection'] = False

    es3 = AdaptiveES(AdaptiveESParams(**base_params))
    es4 = SelfAdaptiveES(ESParams(**base_params))

    params = [
        (es1, es2, 'compare-es-plus'),
        (es3, es4, 'compare-es-comma'),
    ]

    # measure_and_save(params, sequence)

    plot('compare-es-plus', 'PEPTIDE - (µ+λ)', 'Adaptive ES', 'Self-adaptive ES')
    plot('compare-es-comma', 'PEPTIDE - (µ,λ)', 'Adaptive ES', 'Self-adaptive ES')


if __name__ == "__main__":
    main()