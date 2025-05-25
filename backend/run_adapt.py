from typing import Literal

from backend.algorithms.adaptive_es import AdaptiveES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams


def get_es(force_field: Literal['charmm', 'amber']) -> AdaptiveES:
    return AdaptiveES(get_params(force_field))


def get_params(ff: str):
    return AdaptiveESParams(
        force_field=ff,
        plus_selection= {'amber': True,      'charmm': ...}[ff],
        population_size={'amber': 170,       'charmm': ...}[ff],
        children_size=  {'amber': 848,       'charmm': ...}[ff],
        sigma=          {'amber': 53.698234, 'charmm': ...}[ff],
        theta=          {'amber': 0.239167,  'charmm': ...}[ff],
        alpha=          {'amber': 1.381874,  'charmm': ...}[ff],
        mod_frequency=  {'amber': 4,         'charmm': ...}[ff],
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