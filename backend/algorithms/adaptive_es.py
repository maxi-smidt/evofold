import numpy as np

from typing import List
from multiprocessing import get_context
from overrides import overrides

from backend.algorithms.es import ES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams
from backend.structure.protein import Protein, AngleList


class AdaptiveES(ES):
    def __init__(self, params: AdaptiveESParams):
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

    def _adaptive_adaption(self, sigma, success_rate: float) -> float:
        if np.isclose(success_rate, self._params.theta):
            return sigma
        if success_rate < self._params.theta:
            return sigma / self._params.alpha
        return sigma * self._params.alpha

    def _process_child(self, parent: Protein, sigma: float):
        child = self._mutate_protein(parent, sigma)
        return child, child.fitness < parent.fitness

    @overrides
    def run(self, sequence: str, callback: ES.Callback=None, callback_frequency: int=1) -> Protein:
        self._initialize_run(sequence, callback_frequency)
        successes: float = 0.0
        sigma:     float = self._params.sigma

        with get_context("spawn").Pool() as pool:
            while self._generation < self._params.generations:
                self._generation += 1
                children: List[Protein] = self._population if self._params.plus_selection else []
                parents_to_mutate = [self._population[np.random.randint(self._params.population_size)] for _ in range(self._params.children_size)]

                results = pool.starmap(self._process_child, [(parent, sigma) for parent in parents_to_mutate])
                successes += sum(success for _, success in results)
                children.extend(child for child, _ in results)

                if self._finalize_generation(children, callback, sigma):
                    return self._global_best_offspring

                if self._generation % self._params.mod_frequency == 0:
                    sigma = self._adaptive_adaption(sigma, successes / (self._params.mod_frequency * self._params.children_size))
                    successes = 0

        return self._global_best_offspring
