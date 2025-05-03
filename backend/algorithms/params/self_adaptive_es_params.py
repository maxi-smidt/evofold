from dataclasses import dataclass
from typing import Literal

from backend.algorithms.params.es_params import ESParams


@dataclass
class SelfAdaptiveESParams(ESParams):
    strategy_param: Literal['gene-wise', 'genome-wise'] = 'gene-wise' # the way the strategy parameter sigma is tracked