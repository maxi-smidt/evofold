from dataclasses import dataclass

from backend.algorithms.params.es_params import ESParams


@dataclass
class DerandomizedESParams(ESParams):
    tau:           float = 0.2   # threshold for adaptive adaption (1 / 5 success rule)
    alpha:         float = 1.224 # modification factor