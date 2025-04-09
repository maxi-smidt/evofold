import math

import numpy as np

from typing import List, Callable
from multiprocessing import get_context
from overrides import overrides

from backend.algorithms.es import ES
from backend.algorithms.params.derandomized_es_params import DerandomizedESParams
from backend.structure.protein import Protein, AngleList


class DerandomizedES(ES):
    def __init__(self, params: DerandomizedESParams):
        super().__init__(params)

    def _mutate_protein(self, protein: Protein, sigma: float) -> Protein:
        return Protein(protein.sequence, self._params.force_field, angles=self._gaussian_mutation(protein.angles, sigma))

    @staticmethod
    def _gaussian_mutation(angles: AngleList, sigma: float) -> AngleList:
        mutations = np.random.normal(0, sigma, size=(len(angles), 2))
        mutated_angles = [
            (
                np.clip(phi + d_phi, Protein.ANGLE_MIN, Protein.ANGLE_MAX),
                np.clip(psi + d_psi, Protein.ANGLE_MIN, Protein.ANGLE_MAX),
                omega
            )
            for (phi, psi, omega), (d_phi, d_psi) in zip(angles, mutations)
        ]
        return mutated_angles

    @staticmethod
    def _average_angles(angle_lists: List[List[float]]) -> List[float]:
        summed_angles = [0 for _ in range(len(angle_lists[0]))]
        for angel_list in angle_lists:
            for i, angle in enumerate(angel_list):
                summed_angles[i] += angle

        return [angle / len(summed_angles) for angle in summed_angles]

    def _global_arithmetic_crossover(self, population: List[Protein]) -> Protein:
        return Protein(population[0].sequence, self._params.force_field, flat_angles=self._average_angles([a.angles_flat for a in population]))

    @staticmethod
    def _apply_to_each_angle(angle_list: List[float], operation: Callable[[float], float]) -> List[float]:
        return [operation(angle) for angle in angle_list]

    @staticmethod
    def _add_angle_lists(angle_list1: List[float], angle_list2: List[float]) -> List[float]:
        return [a1+a2 for a1, a2 in zip(angle_list1, angle_list2)]

    @staticmethod
    def _distance_of_angle_list(angle_list: List[float]) -> float:
        return math.sqrt(sum(angle ** 2 for angle in angle_list)) # leave omega

    @overrides
    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None]=None, callback_frequency: int=1) -> Protein:
        generation:     int           = 0
        sigma:          float         = self._params.sigma
        population:     List[Protein] = self._create_initial_population(sequence)
        best_offspring: Protein       = min(population, key=lambda p: p.fitness)
        genome_length:  int = len(sequence) * 3
        direction_vec:  List[float]   = [0 for _ in range(genome_length)]

        with get_context("spawn").Pool() as pool:
            while generation < self._params.generations:
                generation += 1
                children: List[Protein] = population if self._params.plus_selection else []
                global_crossover = self._global_arithmetic_crossover(population)

                results = pool.starmap(self._mutate_protein, [(global_crossover, sigma) for _ in range(self._params.children_size)])
                children.extend(results)

                population = self._make_selection(children)

                mutation_directions: List[List[float]] = [
                    [
                        pa - gca for pa, gca in
                        zip(protein.angles_flat, global_crossover.angles_flat)
                    ]
                    for protein in population
                ]

                avg_mutation_direction = self._average_angles(mutation_directions)

                direction_factor = 1 - self._params.alpha
                avg_mutation_direction_factor = math.sqrt(self._params.alpha * (2 - self._params.alpha) * self._params.population_size)
                direction_vec: List[float] = self._add_angle_lists(
                    self._apply_to_each_angle(direction_vec, lambda a: a * direction_factor),
                    self._apply_to_each_angle(avg_mutation_direction, lambda a: a * avg_mutation_direction_factor)
                )

                sigma = sigma * math.exp((self._distance_of_angle_list(direction_vec) ** 2 - genome_length) / (2 * self._params.tau * genome_length))
                best_offspring = min(population, key=lambda p: p.fitness)

                is_premature_termination = self._is_premature_termination(best_offspring)

                if callback is not None and generation % callback_frequency == 0:
                    callback(generation, best_offspring, sigma, is_premature_termination or (generation >= self._params.generations))
                if is_premature_termination:
                    return best_offspring

        return best_offspring