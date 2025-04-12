import math

import numpy as np

from typing import List, Callable, Tuple
from multiprocessing import get_context
from overrides import overrides

from backend.algorithms.es import ES
from backend.algorithms.params.derandomized_es_params import DerandomizedESParams
from backend.structure.protein import Protein, AngleList

"""
Formulas used marked with >>f. n<< are used from https://doi.org/10.48550/arXiv.1604.00772
"""
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

    @staticmethod
    def _calculate_µ_eff_and_µ_eff_neg(weights: np.array, population_size: int) -> Tuple[float, float]:
        w_pos = weights[:population_size]
        w_neg = weights[population_size:]
        return sum(w_pos) ** 2 / sum(w_pos ** 2), sum(w_neg) ** 2 / sum(w_neg ** 2)

    @staticmethod
    def _calculate_initial_weights(children_size: int) -> np.ndarray:
        return np.array([math.log((children_size + 1) / 2) - math.log(i + 1)                                   # f. 49
                         for i in range(children_size)])

    @staticmethod
    def _calculate_weights(w: np.array, population_size: int, µ_eff: float, µ_eff_neg: float, c_1: float, c_µ :float, n: int) -> np.array:
        α_µ_neg = 1 + c_1 / c_µ                                                                                # f. 50
        print('α_µ_neg: ', α_µ_neg)
        α_µ_eff_neg = 1 + (2 * µ_eff_neg) / (µ_eff + 2)                                                        # f. 51
        α_pos_def_neg = (1 - c_1 - c_µ) / (n * c_µ)                                                            # f. 52

        sum_pos = sum(w[:population_size])
        sum_neg = sum(abs(x) for x in w[population_size:])

        return np.array([
            x / sum_pos if i < population_size else abs(x) * min(α_µ_neg, α_µ_eff_neg, α_pos_def_neg) / sum_neg                      # f. 53
            for i, x in enumerate(w)
        ])

    @overrides
    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None]=None, callback_frequency: int=1) -> Protein:
        generation:     int           = 0
        sigma:          float         = self._params.sigma
        population:     List[Protein] = self._create_initial_population(sequence)
        best_offspring: Protein       = min(population, key=lambda p: p.fitness)

        n = len(sequence) * 3
        w = self._calculate_initial_weights(self._params.children_size)                                        # f. 49
        α_cov = 2
        µ_eff, µ_eff_neg = self._calculate_µ_eff_and_µ_eff_neg(w, self._params.population_size)
        c_1 = α_cov / ((n + 1.3) ** 2 + µ_eff)                                                                 # f. 57
        c_µ = min(1 - c_1, α_cov * (0.25 + µ_eff + 1 / µ_eff - 2) / ((n + 2) ** 2 + α_cov * µ_eff / 2))        # f. 58

        w = self._calculate_weights(w, self._params.population_size, µ_eff, µ_eff_neg, c_1, c_µ, n)
        c_σ = (µ_eff + 2) / (n + µ_eff + 5)                                                                    # f. 54.5
        d_σ = 1 + 2 * max(0.0, math.sqrt((µ_eff - 1) / (n + 1)) - 1) + c_σ                                     # f. 55
        c_c = (4 + µ_eff / n) / (n + 4 + 2 * µ_eff / n)                                                        # f. 56
        c_m = 1 # learning rate for the mean
        chiN = n ** 0.5 * (1 - 1 / (4 * n) + 1 / (21 * n ** 2)) # expectation of ||N(0,I)|| == norm(randn(N,1))
        p_σ = 0
        p_c = 0
        C = np.identity(n)
        I = np.identity(n)

        print('sum to µ: ', sum(w[:self._params.population_size]))
        print('sum from µ: ', sum(w[self._params.population_size:]))
        σ = self._params.sigma
        m = ...
        return

        with get_context("spawn").Pool() as pool:
            while generation < self._params.generations:
                generation += 1

                z = [np.random.multivariate_normal(0, I) for _ in range(self._params.children_size)]     # f. 38
                B, D = np.linalg.eigh(C)
                y = [z_k * B * D for z_k in z]                                                                 # f. 39
                x = [m + σ * y_k for m, y_k in y]                                                              # f. 40

                y_w = ...
                m = m + c_m * σ * y_w

                p_σ = (1 - c_σ) * p_σ + math.sqrt(c_σ * (2 - c_σ) * µ_eff) * np.diag(1.0 / np.sqrt(C)) * y_w   # f. 43
                σ = σ * math.exp(c_σ / d_σ * (np.linalg.norm(p_σ) / chiN - 1))

                results = pool.starmap(self._mutate_protein,
                                       [(_, sigma) for _ in range(self._params.children_size)])
                ....extend(results)
                best_offspring = min(population, key=lambda p: p.fitness)

                is_premature_termination = self._is_premature_termination(best_offspring)

                if callback is not None and generation % callback_frequency == 0:
                    callback(generation, best_offspring, sigma, is_premature_termination or (generation >= self._params.generations))
                if is_premature_termination:
                    return best_offspring

        return best_offspring

if __name__ == '__main__':
    DerandomizedES(DerandomizedESParams()).run("PEPTIDE")