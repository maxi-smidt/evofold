from dataclasses import dataclass

from backend.algorithms.params.es_params import ESParams


@dataclass
class AdaptiveESParams(ESParams):
    sigma:         int   = 36    # the initial sigma for gaussian mutation
    theta:         float = 0.2   # threshold for adaptive adaption (1 / 5 success rule)
    alpha:         float = 1.224 # modification factor
    mod_frequency: int   = 2     # the frequency of how often sigma should be adapted