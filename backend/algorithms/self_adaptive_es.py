import math
import numpy as np

from typing import List, Tuple
from multiprocessing import get_context
from overrides import overrides

from backend.algorithms.es import ES
from backend.algorithms.params.self_adaptive_es_params import SelfAdaptiveESParams
from backend.structure.protein import Protein


class SelfAdaptiveES(ES):
    def __init__(self, params: SelfAdaptiveESParams):
        super().__init__(params)

    def _mutate_protein(self, protein: Protein) -> Protein:
        angles, sigma = self._self_adaptive_gaussian_mutation(protein)
        return Protein(protein.sequence, self._params.force_field, flat_angles=angles, sigma=sigma)

    @staticmethod
    def _self_adaptive_gaussian_mutation(protein: Protein) -> Tuple[np.array, np.array]:
        if len(protein.sigma) == 1:
            return SelfAdaptiveES._genome_wise_strategy_param(protein)
        return SelfAdaptiveES._gene_wise_strategy_param(protein)

    @staticmethod
    def _gene_wise_strategy_param(protein: Protein) -> Tuple[np.array, np.array]:
        u = np.random.normal(0, 1)
        sigmas = []
        mutated_angles = []
        l = len(protein.sigma)
        for sigma, angle in zip(protein.sigma, protein.angles_flat):
            u_i = np.random.normal(0, 1)
            s = sigma * math.exp(1 / math.sqrt(2 * l) * u + 1 / math.sqrt(2 * math.sqrt(l)) * u_i)
            mutated_angles.append(np.clip(angle + np.random.normal(0, s), Protein.ANGLE_MIN, Protein.ANGLE_MAX))
            sigmas.append(s)
        return np.array(mutated_angles), np.array(sigmas)

    @staticmethod
    def _genome_wise_strategy_param(protein: Protein) -> Tuple[np.array, np.array]:
        sigma = protein.sigma[0] * math.exp((1 / math.sqrt(len(protein.angles_flat))) * np.random.normal(0, 1))
        mutated_angles = np.clip(
            protein.angles_flat + np.random.normal(0, sigma, len(protein.angles_flat)),
            Protein.ANGLE_MIN,
            Protein.ANGLE_MAX
        )
        return mutated_angles, np.array([sigma])

    def _initialize_sigma(self, sequence):
        for protein in self._population:
            if self._params.strategy_param == 'genome-wise':
                protein.sigma = np.array([self._params.sigma])
            else:
                protein.sigma = np.full(len(sequence) * 2, self._params.sigma)

    @overrides
    def run(self, sequence: str, callback: ES.Callback=None, callback_frequency: int=1) -> Protein:
        self._initialize_run(sequence, callback_frequency)

        self._initialize_sigma(sequence)

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
                        self._initialize_sigma(sequence)

        return self._global_best_offspring