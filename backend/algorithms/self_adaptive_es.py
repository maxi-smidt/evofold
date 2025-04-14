import math

import numpy as np

from typing import List, Callable, Tuple
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
        return Protein(protein.sequence, self._params.force_field, angles=angles, sigma=sigma)

    @staticmethod
    def _self_adaptive_gaussian_mutation(protein: Protein) -> Tuple[AngleList, float]:
        sigma = protein.sigma * math.exp((1 / math.sqrt(len(protein.angles))) * np.random.normal(0, 1))
        mutations = np.random.normal(0, sigma, size=(len(protein.angles), 2))
        mutated_angles = [
            (
                np.clip(phi + d_phi, Protein.ANGLE_MIN, Protein.ANGLE_MAX),
                np.clip(psi + d_psi, Protein.ANGLE_MIN, Protein.ANGLE_MAX),
                omega
            )
            for (phi, psi, omega), (d_phi, d_psi) in zip(protein.angles, mutations)
        ]
        return mutated_angles, sigma

    @overrides
    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None]=None, callback_frequency: int=1) -> Protein:
        generation:     int           = 0
        population:     List[Protein] = self._create_initial_population(sequence)
        best_offspring: Protein       = min(population, key=lambda p: p.fitness)

        for protein in population:
            protein.sigma = self._params.sigma

        with get_context("spawn").Pool() as pool:
            while generation < self._params.generations:
                generation += 1
                children: List[Protein] = population if self._params.plus_selection else []
                parents_to_mutate = [population[np.random.randint(self._params.population_size)] for _ in range(self._params.children_size)]

                results = pool.map(self._mutate_protein, parents_to_mutate)
                children.extend(results)

                population = self._make_selection(children)

                best_offspring = min(population, key=lambda p: p.fitness)

                is_premature_termination = self._is_premature_termination(best_offspring)

                if callback is not None and generation % callback_frequency == 0:
                    callback(generation, best_offspring, best_offspring.sigma, is_premature_termination or (generation >= self._params.generations))
                if is_premature_termination:
                    return best_offspring

        return best_offspring