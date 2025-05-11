import math

from typing import Tuple, Literal

from backend.algorithms.adaptive_es import AdaptiveES
from backend.algorithms.derandomized_es import DerandomizedES
from backend.algorithms.es import ES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams
from backend.algorithms.params.derandomized_es_params import DerandomizedESParams
from backend.algorithms.params.self_adaptive_es_params import SelfAdaptiveESParams
from backend.algorithms.self_adaptive_es import SelfAdaptiveES


def get_es(force_field: Literal['charmm', 'amber'], plus_selection: bool, sequence: str) -> Tuple[ES, ES, ES]:
    adaptive_params = AdaptiveESParams(force_field=force_field, plus_selection=plus_selection)
    self_adaptive_params = SelfAdaptiveESParams(force_field=force_field, plus_selection=plus_selection)
    derandomized_params = DerandomizedESParams(force_field=force_field, plus_selection=plus_selection)
    derandomized_params.tau = math.sqrt(len(sequence) * 2)
    derandomized_params.alpha = 1 / math.sqrt(len(sequence) * 2)
    return AdaptiveES(adaptive_params), SelfAdaptiveES(self_adaptive_params), DerandomizedES(derandomized_params)

def map_selection(selection: bool):
    return 'plus' if selection else 'comma'

def main():
    sequence: str = 'PEPTIDE'

    configs = [
        ['amber', True],
        ['charmm', True],
        ['amber', False],
        ['charmm', False],
    ]
    start = 0
    runs = 20
    for ff, selection in configs:
        print(f'Running {ff} {map_selection(selection)}...')
        for i in range(start, runs):
            adaptive, self_adaptive, derandomized = get_es(ff, selection, sequence)
            p_ada = adaptive.run(sequence)
            p_self_ada = self_adaptive.run(sequence)
            p_derand = derandomized.run(sequence)
            print(f'{p_ada.fitness},{p_self_ada},{p_derand}')

        print(f'----')


if __name__ == "__main__":
    main()