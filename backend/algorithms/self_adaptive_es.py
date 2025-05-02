import math

import numpy as np

from typing import List, Tuple
from multiprocessing import get_context
from overrides import overrides

from backend.algorithms.es import ES
from backend.algorithms.params.es_params import ESParams
from backend.structure.protein import Protein, AngleList


class SelfAdaptiveES(ES):
    def __init__(self, params: ESParams):
        super().__init__(params)

    def _mutate_protein(self, protein: Protein) -> Protein:
        angles, sigma = self._self_adaptive_gaussian_mutation(protein)
        return Protein(protein.sequence, self._params.force_field, flat_angles=angles, sigma=sigma)

    @staticmethod
    def _self_adaptive_gaussian_mutation(protein: Protein) -> Tuple[AngleList, float]:
        sigma = protein.sigma * math.exp((1 / math.sqrt(len(protein.angles_flat))) * np.random.normal(0, 1))
        mutated_angles = np.clip(
            protein.angles_flat + np.random.normal(0, sigma, len(protein.angles_flat)),
            Protein.ANGLE_MIN,
            Protein.ANGLE_MAX
        )
        return mutated_angles, sigma

    @overrides
    def run(self, sequence: str, callback: ES.Callback, callback_frequency: int=1) -> Protein:
        self._initialize_run(sequence, callback_frequency)

        for protein in self._population:
            protein.sigma = self._params.sigma

        with get_context("spawn").Pool() as pool:
            while self._generation < self._params.generations:
                self._generation += 1
                children: List[Protein] = self._population if self._params.plus_selection else []
                parents_to_mutate = [self._population[np.random.randint(self._params.population_size)] for _ in range(self._params.children_size)]

                results = pool.map(self._mutate_protein, parents_to_mutate)
                children.extend(results)

                if self._finalize_generation(children, callback, None):
                    if self._params.premature_strategy == 'terminate': return self._global_best_offspring
                    if self._params.premature_strategy == 'restart':
                        self._population = self._create_initial_population(sequence)
                        for protein in self._population:
                            protein.sigma = self._params.sigma

        return self._global_best_offspring