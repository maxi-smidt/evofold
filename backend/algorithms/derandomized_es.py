import numpy as np

from multiprocessing import get_context
from overrides import overrides
from typing import List, Callable

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

    @staticmethod
    def _global_arithmetic_crossover(population: List[Protein]) -> np.array:
        return sum(p.angles_flat for p in population) / len(population)

    @staticmethod
    def _global_uniform_crossover(population: List[Protein]) -> np.array:
        l = len(population[0].angles_flat)
        all_genomes = np.array([p.angles_flat for p in population])
        selection = np.random.randint(len(population), size=l)
        return all_genomes[selection, np.arange(l)]

    @overrides
    def run(self, sequence: str, callback: ES.Callback=None, callback_frequency: int=1) -> Protein:
        self._initialize_run(sequence, callback_frequency)
        sigma:         float    = self._params.sigma
        genome_length: int      = len(sequence) * 2
        s:             np.array = np.zeros(genome_length)
        crossover:     Callable = self._global_arithmetic_crossover if self._params.crossover == 'arithmetic' else self._global_uniform_crossover

        with get_context("spawn").Pool() as pool:
            while self._generation < self._params.generations:
                self._generation += 1
                children = self._population if self._params.plus_selection else []

                mutated = crossover(self._population)

                results = pool.starmap(self._process_child, [(sequence, mutated, sigma) for _ in range(self._params.children_size)])
                children.extend(results)

                if self._finalize_generation(children, callback, sigma):
                    if self._params.premature_strategy == 'terminate': return self._global_best_offspring
                    if self._params.premature_strategy == 'restart':
                        sigma = self._params.sigma
                        s = np.zeros(genome_length)
                        self._population = self._create_initial_population(sequence)
                else:
                    z_mean = sum((p.angles_flat - mutated) / sigma  for p in self._population) / len(self._population)
                    s = (1 - self._params.alpha) * s + np.sqrt(self._params.alpha * (2 - self._params.alpha) * self._params.population_size) * z_mean
                    sigma = sigma * np.exp((np.linalg.norm(s) ** 2 - genome_length) / (2 * self._params.tau * genome_length))

        return self._global_best_offspring