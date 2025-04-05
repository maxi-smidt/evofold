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
        return Protein(protein.sequence, self._params.force_field, self._gaussian_mutation(protein.angles, sigma))

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
    def _average_angles(angle_lists: List[AngleList]) -> AngleList:
        summed_angles = [[a[0], a[1], a[2]] for a in angle_lists[0]]
        for angle_list in angle_lists[1:]:
            for new_angles, angles in zip(angle_list, summed_angles):
                angles[0] += new_angles[0]
                angles[1] += new_angles[1]
                angles[2] += new_angles[2]

        return [
            (angle[0] / len(angle_lists), angle[1] / len(angle_lists), angle[2] / len(angle_lists))
            for angle in summed_angles
        ]

    def _global_arithmetic_crossover(self, population: List[Protein]) -> Protein:
        return Protein(population[0].sequence, self._params.force_field, self._average_angles([a.angles for a in population]))

    @staticmethod
    def _apply_to_each_angle(angle_list: AngleList, operation: Callable[[float], float]) -> AngleList:
        return [(operation(angles[0]), operation(angles[1]), operation(angles[2])) for angles in angle_list]

    @staticmethod
    def _add_angle_lists(angle_list1: AngleList, angle_list2: AngleList) -> AngleList:
        return [(a1[0]+a2[0], a1[1]+a2[1], a1[2]+a2[2]) for a1, a2 in zip(angle_list1, angle_list2)]

    @staticmethod
    def _distance_of_angle_list(angle_list: AngleList) -> float:
        return math.sqrt(sum(angles[0] ** 2 + angles[1] ** 2 for angles in angle_list)) # leave omega

    @overrides
    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None]=None, callback_frequency: int=1) -> Protein:
        generation:     int           = 0
        sigma:          float         = self._params.sigma
        population:     List[Protein] = self._create_initial_population(sequence)
        best_offspring: Protein       = min(population, key=lambda p: p.fitness)
        direction_vec:  List          = [(0, 0, 0) for _ in range(len(sequence))] # s

        with get_context("spawn").Pool() as pool:
            while generation < self._params.generations:
                generation += 1
                children: List[Protein] = population if self._params.plus_selection else []
                global_crossover = self._global_arithmetic_crossover(population)

                results = pool.starmap(self._mutate_protein, [(global_crossover, sigma) for _ in range(self._params.children_size)])
                children.extend(results)

                population = self._make_selection(children)

                mutation_directions: List[AngleList] = [
                    [
                        (pa[0] - gca[0], pa[1] - gca[1], pa[2] - gca[2]) for pa, gca in
                        zip(protein.angles, global_crossover.angles)
                    ]
                    for protein in population
                ]

                avg_mutation_direction = self._average_angles(mutation_directions)

                direction_factor = 1 - self._params.alpha
                avg_mutation_direction_factor = math.sqrt(self._params.alpha * (2 - self._params.alpha) * self._params.population_size)
                direction_vec: AngleList = self._add_angle_lists(
                    self._apply_to_each_angle(direction_vec, lambda a: a * direction_factor),
                    self._apply_to_each_angle(avg_mutation_direction, lambda a: a * avg_mutation_direction_factor)
                )

                genome_length = len(sequence) * 2 # phi und psi for each aa
                sigma = sigma * math.exp((self._distance_of_angle_list(direction_vec) ** 2 - genome_length) / (2 * self._params.tau * genome_length))
                best_offspring = min(population, key=lambda p: p.fitness)

                is_premature_termination = self._is_premature_termination(best_offspring)

                if callback is not None and generation % callback_frequency == 0:
                    callback(generation, best_offspring, sigma, is_premature_termination or (generation >= self._params.generations))
                if is_premature_termination:
                    return best_offspring

        return best_offspring