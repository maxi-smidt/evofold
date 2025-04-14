import numpy as np

from multiprocessing import get_context
from overrides import overrides
from typing import Callable, List

from backend.algorithms.es import ES
from backend.algorithms.params.derandomized_es_params import DerandomizedESParams
from backend.structure.protein import Protein


class DerandomizedES(ES):
    def __init__(self, params: DerandomizedESParams):
        super().__init__(params)

    @staticmethod
    def _gaussian_mutation(angles: np.array, sigma: float) -> np.array:
        mutations = np.random.normal(0, sigma, size=len(angles))
        return [np.clip(angle + mutations, Protein.ANGLE_MIN, Protein.ANGLE_MAX)
                for angle, mutations in zip(angles, mutations)]

    def _process_child(self, sequence: str, mean: np.array, sigma: float):
        angles = DerandomizedES._gaussian_mutation(mean, sigma)
        return Protein(sequence, self._params.force_field, flat_angles=angles)

    @overrides
    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None]=None, callback_frequency: int=1) -> Protein:
        generation:     int           = 0
        sigma:          float         = self._params.sigma
        population:     List[Protein] = self._create_initial_population(sequence)
        best_offspring: Protein       = min(population, key=lambda p: p.fitness)
        l:              int           = len(sequence) * 2
        s:              np.array      = np.zeros(l)

        with get_context("spawn").Pool() as pool:
            while generation < self._params.generations:
                generation += 1

                mean = sum(p.angles_flat for p in population) / len(population)
                children = population if self._params.plus_selection else []

                results = pool.starmap(self._process_child, [(sequence, mean, sigma) for _ in range(self._params.children_size)])
                children.extend(results)

                population = self._make_selection(children)

                z_mean = sum((p.angles_flat - mean) / sigma  for p in population) / len(population)

                s = (1 - self._params.alpha) * s + np.sqrt(self._params.alpha * (2 - self._params.alpha) * self._params.population_size) * z_mean

                sigma = sigma * np.exp((np.linalg.norm(s) ** 2 - l) / (2 * self._params.tau * l))

                best_offspring = min(population, key=lambda p: p.fitness)

                is_premature_termination = self._is_premature_termination(best_offspring)

                if callback is not None and generation % callback_frequency == 0:
                    callback(generation, best_offspring, sigma, is_premature_termination or (generation >= self._params.generations))
                if is_premature_termination:
                    return best_offspring

        return best_offspring