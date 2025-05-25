from typing import Literal

from backend.algorithms.adaptive_es import AdaptiveES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams


def get_es(force_field: Literal['charmm', 'amber']) -> AdaptiveES:
    return AdaptiveES(get_params(force_field))


def get_params(ff: str):
    return AdaptiveESParams(
        force_field=ff,
        plus_selection= {'amber': True,      'charmm': True}[ff],
        population_size={'amber': 170,       'charmm': 129}[ff],
        children_size=  {'amber': 848,       'charmm': 955}[ff],
        sigma=          {'amber': 53.698234, 'charmm': 36.882018}[ff],
        theta=          {'amber': 0.239167,  'charmm': 0.391328}[ff],
        alpha=          {'amber': 1.381874,  'charmm': 1.286563}[ff],
        mod_frequency=  {'amber': 4,         'charmm': 7}[ff],
    )


def main():
    sequence: str = 'PEPTIDE'

    start = 0
    runs = 20
    for ff in ['amber', 'charmm']:
        print(f'Running {ff}...')
        for i in range(start, runs):
            adaptive = get_es(ff)
            p_ada = adaptive.run(sequence)
            print(f'{p_ada.fitness}')

        print(f'----')


if __name__ == "__main__":
    main()