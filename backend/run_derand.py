import math

from typing import Literal

from backend.algorithms.derandomized_es import DerandomizedES
from backend.algorithms.params.derandomized_es_params import DerandomizedESParams


def get_es(force_field: Literal['charmm', 'amber'], sequence: str) -> DerandomizedES:
    derandomized_params = get_params(force_field)
    derandomized_params.tau = math.sqrt(len(sequence) * 2)
    derandomized_params.alpha = 1 / math.sqrt(len(sequence) * 2)
    return DerandomizedES(derandomized_params)


def get_params(ff: str):
    return DerandomizedESParams(
        force_field=ff,
        plus_selection= {'amber': False,               'charmm': True}[ff],
        population_size={'amber': 1,                   'charmm': 14}[ff],
        children_size=  {'amber': 645,                 'charmm': 978}[ff],
        sigma=          {'amber': 178.546855,          'charmm': 178.929899}[ff],
        crossover=      {'amber': 'global-arithmetic', 'charmm': 'global-uniform'}[ff],
        tau=            {'amber': 6.066264,            'charmm': 13.342222}[ff],
        alpha=          {'amber': 0.901668,            'charmm': 0.287913}[ff],
    )


def main():
    sequence: str = 'PEPTIDE'

    start = 0
    runs = 20
    for ff in ['amber', 'charmm']:
        print(f'Running {ff}...')
        for i in range(start, runs):
            derandomized = get_es(ff, sequence)
            p_derand = derandomized.run(sequence)
            print(f'{p_derand.fitness}')

        print(f'----')


if __name__ == "__main__":
    main()