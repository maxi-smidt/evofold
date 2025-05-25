import math

from typing import Tuple, Literal

from backend.algorithms.adaptive_es import AdaptiveES
from backend.algorithms.derandomized_es import DerandomizedES
from backend.algorithms.es import ES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams
from backend.algorithms.params.derandomized_es_params import DerandomizedESParams
from backend.algorithms.params.self_adaptive_es_params import SelfAdaptiveESParams
from backend.algorithms.self_adaptive_es import SelfAdaptiveES


def get_es(force_field: Literal['charmm', 'amber'], sequence: str) -> Tuple[ES, ES, ES]:
    adaptive_params, self_adaptive_params, derandomized_params = get_params(force_field)
    derandomized_params.tau = math.sqrt(len(sequence) * 2)
    derandomized_params.alpha = 1 / math.sqrt(len(sequence) * 2)
    return AdaptiveES(adaptive_params), SelfAdaptiveES(self_adaptive_params), DerandomizedES(derandomized_params)


def get_params(ff: str):
    p1 = AdaptiveESParams(
        force_field=ff,
        plus_selection= {'amber': True,      'charmm': ...}[ff],
        population_size={'amber': 170,       'charmm': ...}[ff],
        children_size=  {'amber': 848,       'charmm': ...}[ff],
        sigma=          {'amber': 53.698234, 'charmm': ...}[ff],
        theta=          {'amber': 0.239167,  'charmm': ...}[ff],
        alpha=          {'amber': 1.381874,  'charmm': ...}[ff],
        mod_frequency=  {'amber': 4,         'charmm': ...}[ff],
    ),
    p2 = SelfAdaptiveESParams(
        force_field=ff,
        plus_selection= {'amber': True,        'charmm': True}[ff],
        population_size={'amber': 134,         'charmm': 195}[ff],
        children_size=  {'amber': 955,         'charmm': 829}[ff],
        sigma=          {'amber': 30.212228,   'charmm': 70.164387}[ff],
        strategy_param= {'amber': 'gene-wise', 'charmm': 'genome-wise'}[ff],
    ),
    p3 = DerandomizedESParams(
        force_field=ff,
        plus_selection= {'amber': False,               'charmm': True}[ff],
        population_size={'amber': 1,                   'charmm': 14}[ff],
        children_size=  {'amber': 645,                 'charmm': 978}[ff],
        sigma=          {'amber': 178.546855,          'charmm': 178.929899}[ff],
        crossover=      {'amber': 'global-arithmetic', 'charmm': 'global-uniform'}[ff],
        tau=            {'amber': 6.066264,            'charmm': 13.342222}[ff],
        alpha=          {'amber': 0.901668,            'charmm': 0.287913}[ff],
    )
    return p1, p2, p3


def main():
    sequence: str = 'PEPTIDE'

    start = 0
    runs = 20
    for ff in ['amber', 'charmm']:
        print(f'Running {ff}...')
        for i in range(start, runs):
            adaptive, self_adaptive, derandomized = get_es(ff, sequence)
            #p_ada = adaptive.run(sequence)
            #p_self_ada = self_adaptive.run(sequence)
            p_derand = derandomized.run(sequence)
            # print(f'{p_ada.fitness},{p_self_ada.fitness},{p_derand.fitness}')
            print(f'{p_derand.fitness}')

        print(f'----')


if __name__ == "__main__":
    main()