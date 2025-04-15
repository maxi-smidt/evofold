from dataclasses import dataclass
from typing import Optional, Literal

StrategyType = Literal['none', 'terminate', 'restart']
ForceFieldType = Literal['amber', 'charmm']

@dataclass
class ESParams:
    generations:           int             = 500     # the amount of generations that will be produced
    population_size:       int             = 100     # the size of the population
    children_size:         int             = 600     # the amount of children produced by one population
    plus_selection:        bool            = True    # if true (µ+λ) else (µ,λ)
    force_field:           ForceFieldType  = 'amber' # the force field that is used for energy evaluation
    sigma:                 int             = 36      # the initial sigma for gaussian mutation
    premature_strategy:    StrategyType    = 'none'  # strategy that is used if any of the premature cases occur
    premature_stagnation:  Optional[int]   = 15      # if not none it initiates the chosen strategy if the best fitness has not changed for n generations
    premature_sigma:       Optional[float] = 0.01    # if not none it initiates the chosen strategy if sigma falls below n
    premature_fitness:     Optional[float] = 2000    # if not none it initiates the chosen strategy if the best fitness falls below n

    def __post_init__(self):
        assert self.population_size <= self.children_size