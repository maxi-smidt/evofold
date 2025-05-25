import math

from typing import Tuple, Literal

from backend.algorithms.adaptive_es import AdaptiveES
from backend.algorithms.derandomized_es import DerandomizedES
from backend.algorithms.es import ES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams
from backend.algorithms.params.derandomized_es_params import DerandomizedESParams
from backend.algorithms.params.self_adaptive_es_params import SelfAdaptiveESParams
from backend.algorithms.self_adaptive_es import SelfAdaptiveES


def get_es(force_field: Literal['charmm', 'amber']) ->SelfAdaptiveES:
    return SelfAdaptiveES(get_params(force_field))


def get_params(ff: str):
    return SelfAdaptiveESParams(
        force_field=ff,
        plus_selection= {'amber': True,        'charmm': True}[ff],
        population_size={'amber': 134,         'charmm': 195}[ff],
        children_size=  {'amber': 955,         'charmm': 829}[ff],
        sigma=          {'amber': 30.212228,   'charmm': 70.164387}[ff],
        strategy_param= {'amber': 'gene-wise', 'charmm': 'genome-wise'}[ff],
    )


def main():
    sequence: str = 'PEPTIDE'

    start = 0
    runs = 20
    for ff in ['amber', 'charmm']:
        print(f'Running {ff}...')
        for i in range(start, runs):
            self_adaptive = get_es(ff)
            p_self_ada = self_adaptive.run(sequence)
            print(f'{p_self_ada.fitness}')

        print(f'----')


if __name__ == "__main__":
    main()