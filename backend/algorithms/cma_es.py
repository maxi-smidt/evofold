import numpy as np

from multiprocessing import get_context
from overrides import overrides
from scipy.linalg import fractional_matrix_power
from typing import Callable, Tuple, override, List, Optional

from backend.algorithms.es import ES
from backend.algorithms.params.derandomized_es_params import DerandomizedESParams
from backend.structure.protein import Protein

"""
The CMA-ES is implemented based on https://doi.org/10.48550/arXiv.1604.00772
"""
class CmaES(ES):
    def __init__(self, params: DerandomizedESParams):
        super().__init__(params)

    def _mutate_protein(self, angles: np.array, sequence: str) -> Protein:
        return Protein(sequence, self._params.force_field, flat_angles=angles)

    @override
    def _finalize_generation(self, children: Optional[List[Protein]], callback: ES.Callback, sigma: float | None) -> bool:
        start_strategy = self._has_reached_stagnation() or self._has_reached_sigma(sigma) or self._has_reached_fitness()
        will_terminate = start_strategy and self._params.premature_strategy == 'terminate'

        if callback is not None and self._generation % self._callback_frequency == 0:
            callback(self._generation, self._local_best_offspring, sigma,
                     will_terminate or (self._generation >= self._params.generations))

        return start_strategy and self._params.premature_strategy != 'none'

    @staticmethod
    def _calculate_µ_eff_and_µ_eff_neg(weights: np.array, population_size: int) -> Tuple[float, float]:
        w_pos = weights[:population_size]
        w_neg = weights[population_size:]
        return sum(w_pos) ** 2 / sum(w_pos ** 2), sum(w_neg) ** 2 / sum(w_neg ** 2)

    @staticmethod
    def _calculate_initial_weights(children_size: int) -> np.ndarray:
        return np.array([np.log((children_size + 1) / 2) - np.log(i + 1)                                         # f. 49
                         for i in range(children_size)])

    @staticmethod
    def _calculate_weights(w: np.array, µ_eff: float, µ_eff_neg: float, c_1: float, c_µ :float, n: int) -> np.array:
        α_µ_neg = 1 + c_1 / c_µ                                                                                  # f. 50
        α_µ_eff_neg = 1 + (2 * µ_eff_neg) / (µ_eff + 2)                                                          # f. 51
        α_pos_def_neg = (1 - c_1 - c_µ) / (n * c_µ)                                                              # f. 52

        sum_pos = sum(wi for wi in w if wi >= 0)
        sum_neg = sum(abs(wi) for wi in w if wi < 0)

        return np.array([
            (wi / sum_pos) if wi >= 0 else (min(α_µ_neg, α_µ_eff_neg, α_pos_def_neg) / -sum_neg)                  # f. 53
            for wi in w
        ])

    @staticmethod
    def _adapt_weights(w: np.array, y: np.array, n: int, C_inv_sqrt: np.array) -> np.ndarray:
        return np.array([w_i * (1 if w_i >= 0 else (n / np.linalg.norm(C_inv_sqrt @ y_i) ** 2))                  # f. 46
                for w_i, y_i in zip(w, y)])

    @staticmethod
    def eval_h_σ(p_σ: np.array, c_σ: float, generation: int, n: int, chiN: float) -> bool:
        return np.linalg.norm(p_σ) / np.sqrt(1 - (1 - c_σ) ** (2 * (generation + 1))) < (1.4 + 2 / (n + 1) * chiN)

    @staticmethod
    def _adapt_C(C: np.array, δ_h_σ: float, c_1: float, c_µ: float, w: np.array, w_adapted: np.array, p_c: np.array,
                 y: np.array):
        part_1 = 1 + c_1 * δ_h_σ - c_1 - c_µ * sum(w)
        part_2 = c_µ * sum(w_i * y_i * y_i.T for w_i, y_i in zip(w_adapted, y))
        return part_1 * C + c_1 * p_c * p_c.T + part_2                                                           # f. 47

    @overrides
    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None]=None, callback_frequency: int=1) -> Protein:
        generation = 0
        µ = self._params.population_size
        n = len(sequence) * 2
        α_cov = 2
        w = self._calculate_initial_weights(self._params.children_size)                                          # f. 49
        µ_eff, µ_eff_neg = self._calculate_µ_eff_and_µ_eff_neg(w, self._params.population_size)
        c_1 = α_cov / ((n + 1.3) ** 2 + µ_eff)                                                                   # f. 57
        c_µ = min(1 - c_1, α_cov * (0.25 + µ_eff + 1 / µ_eff - 2) / ((n + 2) ** 2 + α_cov * µ_eff / 2))          # f. 58

        w = self._calculate_weights(w, µ_eff, µ_eff_neg, c_1, c_µ, n)
        c_σ = (µ_eff + 2) / (n + µ_eff + 5)                                                                      # f. 54.5
        d_σ = 1 + 2 * max(0.0, np.sqrt((µ_eff - 1) / (n + 1)) - 1) + c_σ                                         # f. 55
        c_c = (4 + µ_eff / n) / (n + 4 + 2 * µ_eff / n)                                                          # f. 56
        c_m = 1 # learning rate for the mean
        chiN = np.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n ** 2)) # expectation of ||N(0,I)|| == norm(randn(N,1))
        p_σ = np.zeros(n)
        p_c = np.zeros(n)
        C = np.identity(n)
        I = np.identity(n)

        σ = self._params.sigma
        m = np.random.uniform(Protein.ANGLE_MIN, Protein.ANGLE_MAX, size=n)
        self._global_best_offspring: Protein = Protein(sequence, self._params.force_field, flat_angles=m)

        with get_context("spawn").Pool() as pool:
            while generation < self._params.generations:
                generation += 1

                z = [np.random.multivariate_normal(np.zeros(n), I) for _ in range(self._params.children_size)]   # f. 38
                eigenvalues, B = np.linalg.eigh(C)
                eigenvalues = np.clip(eigenvalues, 1e-8, None)
                D = np.diag(np.sqrt(eigenvalues))
                BD = B @ D
                y = [BD @ z_k for z_k in z]                                                                      # f. 39
                x = [np.clip(m + σ * y_k, Protein.ANGLE_MIN, Protein.ANGLE_MAX) for y_k in y]                    # f. 40

                results = pool.starmap(self._mutate_protein, [(angles, sequence) for angles in x])
                evaluated = list(zip(results, y))
                evaluated_sorted = sorted(evaluated, key=lambda pair: pair[0].fitness)
                proteins_sorted, y_sorted = zip(*evaluated_sorted)

                self._local_best_offspring = proteins_sorted[0]

                print(c_σ)
                print(d_σ)
                print(c_σ / d_σ)

                print(np.linalg.norm(p_σ)) # this is the problem for the overflow
                print(chiN)
                print((np.linalg.norm(p_σ) / chiN - 1))

                print((c_σ / d_σ) * (np.linalg.norm(p_σ) / chiN - 1))
                print("======")

                if self._finalize_generation(None, callback, σ):
                    if self._params.premature_strategy == 'terminate': return self._global_best_offspring
                    if self._params.premature_strategy == 'restart':
                        σ = self._params.sigma
                        C = np.identity(n)
                        p_σ = np.zeros(n)
                        p_c = np.zeros(n)
                        m = np.random.uniform(Protein.ANGLE_MIN, Protein.ANGLE_MAX, size=n)
                        w = self._calculate_initial_weights(self._params.children_size)
                        w = self._calculate_weights(w, µ_eff, µ_eff_neg, c_1, c_µ, n)
                else:
                    y_w = sum(w_i * y_i for w_i, y_i in zip(w[:µ], y_sorted[:µ]))                                # f. 41
                    m = m + c_m * σ * y_w                                                                        # f. 42
                    C_inv_sqrt = fractional_matrix_power(C, -0.5)
                    p_σ = (1 - c_σ) * p_σ + np.sqrt(c_σ * (2 - c_σ) * µ_eff) * C_inv_sqrt @ y_w  # f. 43
                    print(np.linalg.norm(p_σ))
                    print(self._local_best_offspring.fitness)
                    σ = σ * np.exp((c_σ / d_σ) * (np.linalg.norm(p_σ) / chiN - 1))                           # f. 44
                    h_σ = 1 if self.eval_h_σ(p_σ, c_σ, generation, n, chiN) else 0
                    p_c = (1 - c_c) * p_c + h_σ * np.sqrt(c_c * (2 - c_c) * µ_eff) * y_w                         # f. 45
                    w_adapted = self._adapt_weights(w, y_sorted, n, C_inv_sqrt)                                  # f. 46
                    δ_h_σ = (1 - h_σ) * c_c * (2 - c_c)
                    C = self._adapt_C(C, δ_h_σ, c_1, c_µ, w, w_adapted, p_c, y)                                  # f. 47
                    C = (C + C.T) / 2
                    µ_eff, µ_eff_neg = self._calculate_µ_eff_and_µ_eff_neg(w, self._params.population_size)
                    c_1 = α_cov / ((n + 1.3) ** 2 + µ_eff)
                    c_µ = min(1 - c_1, α_cov * (0.25 + µ_eff + 1 / µ_eff - 2) / ((n + 2) ** 2 + α_cov * µ_eff / 2))

        return self._global_best_offspring

if __name__ == '__main__':
    x = CmaES(DerandomizedESParams(population_size=4, children_size=30)).run("PEPTIDE")
    print(x.fitness)