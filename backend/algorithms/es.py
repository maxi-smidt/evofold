from abc import ABC, abstractmethod
from heapq import nsmallest
from typing import Callable, List

from backend.algorithms.params.es_params import ESParams
from backend.structure.protein import Protein


class ES(ABC):
    Callback = Callable[[int, Protein, float, bool], None]

    def __init__(self, params: ESParams):
        self._params = params
        self._previous_best = None
        self._previous_best_count = 0

    @abstractmethod
    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None] = None, callback_frequency: int = 1) -> Protein:
        ...

    def _initialize_run(self, sequence: str, callback_frequency) -> None:
        self._generation:            int           = 0
        self._population:            List[Protein] = self._create_initial_population(sequence)
        self._local_best_offspring:  Protein       = min(self._population, key=lambda p: p.fitness)
        self._global_best_offspring: Protein       = self._local_best_offspring
        self._callback_frequency:    int           = callback_frequency

    def _finalize_generation(self, children: List[Protein], callback: Callback, sigma: float | None) -> bool:
        self._population = self._make_selection(children)
        self._local_best_offspring = min(self._population, key=lambda p: p.fitness)
        self._global_best_offspring = min(self._local_best_offspring, self._global_best_offspring, key=lambda p: p.fitness)

        sigma = sigma or sum(p.sigma for p in self._population) / len(self._population)

        start_strategy = self._has_reached_stagnation() or self._has_reached_sigma(sigma) or self._has_reached_fitness()
        will_terminate = start_strategy and self._params.premature_strategy == 'terminate'

        if callback is not None and self._generation % self._callback_frequency == 0:
            callback(self._generation, self._local_best_offspring, sigma,
                     will_terminate or (self._generation >= self._params.generations))

        return start_strategy and self._params.premature_strategy != 'none'

    def _make_selection(self, children: List[Protein]) -> List[Protein]:
        fitness: Callable[[Protein], float] = lambda p: p.fitness
        return nsmallest(self._params.population_size, children, fitness)

    def _create_initial_population(self, sequence) -> List[Protein]:
        return [Protein(sequence, self._params.force_field) for _ in range(self._params.population_size)]

    def _has_reached_stagnation(self) -> bool:
        return self._params.premature_stagnation is not None and self._reached_stagnation_inner()

    def _reached_stagnation_inner(self) -> bool:
        if self._previous_best is not None:
            if abs(self._local_best_offspring.fitness - self._previous_best.fitness) < 0.0001:
                self._previous_best_count += 1
            else:
                self._previous_best_count = 1
            if self._previous_best_count == self._params.premature_stagnation:
                return True
        self._previous_best = self._local_best_offspring
        return False

    def _has_reached_sigma(self, sigma: float) -> bool:
        return self._params.premature_sigma is not None and sigma <= self._params.premature_sigma

    def _has_reached_fitness(self) -> bool:
        return self._params.premature_fitness is not None and self._local_best_offspring.fitness <= self._params.premature_sigma