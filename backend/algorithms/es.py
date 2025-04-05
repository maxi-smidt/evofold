from abc import ABC, abstractmethod
from heapq import nsmallest
from typing import Callable, List

from backend.algorithms.params.es_params import ESParams
from backend.structure.protein import Protein


class ES(ABC):
    def __init__(self, params: ESParams):
        self._params = params
        self._previous_best = None
        self._previous_best_count = 0

    @abstractmethod
    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None] = None, callback_frequency: int = 1) -> Protein:
        ...

    def _make_selection(self, children: List[Protein]) -> List[Protein]:
        fitness: Callable[[Protein], float] = lambda p: p.fitness
        return nsmallest(self._params.population_size, children, fitness)

    def _create_initial_population(self, sequence) -> List[Protein]:
        return [Protein(sequence, self._params.force_field) for _ in range(self._params.population_size)]

    def _is_premature_termination(self, best_offspring: Protein) -> bool:
        return self._params.premature_termination is not None and self._reached_premature_termination(best_offspring)

    def _reached_premature_termination(self, best_offspring: Protein) -> bool:
        if self._previous_best is not None:
            if best_offspring.fitness == self._previous_best.fitness:
                self._previous_best_count += 1
            else:
                self._previous_best_count = 1
            if self._previous_best_count == self._params.premature_termination:
                return True
        self._previous_best = best_offspring
        return False