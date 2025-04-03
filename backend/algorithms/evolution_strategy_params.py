from dataclasses import dataclass
from typing import Optional


@dataclass
class EvolutionStrategyParams:
    generations:           int           = 500     # the amount of generations that will be produced
    population_size:       int           = 100     # the size of the population
    children_size:         int           = 600     # the amount of children produced by one population
    plus_selection:        bool          = True    # if true (µ+λ) else (µ,λ)
    force_field:           str           = 'amber' # the force field that is used for energy evaluation
    sigma:                 int           = 36      # the initial sigma for gaussian mutation
    theta:                 float         = 0.2     # threshold for adaptive adaption (1 / 5 success rule)
    alpha:                 float         = 1.224   # modification factor
    premature_termination: Optional[int] = 10      # if not none it terminates the algorithm if the best protein has not changed for n generations

    def __post_init__(self):
        assert self.population_size <= self.children_size