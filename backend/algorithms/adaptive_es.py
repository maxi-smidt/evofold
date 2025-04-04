import numpy as np

from heapq import nsmallest
from typing import List, Callable
from multiprocessing import get_context

from backend.algorithms.es import ES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams
from backend.structure.protein import Protein, AngleList


class AdaptiveES(ES):
    def __init__(self, params: AdaptiveESParams):
        self._params = params
        self._previous_best = None
        self._previous_best_count = 0

    def _make_selection(self, children: List[Protein]) -> List[Protein]:
        fitness: Callable[[Protein], float] = lambda p: p.fitness
        return nsmallest(self._params.population_size, children, fitness)

    def _create_initial_population(self, sequence) -> List[Protein]:
        return [Protein(sequence, self._params.force_field) for _ in range(self._params.population_size)]

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

    def _adaptive_adaption(self, sigma, success_rate: float) -> float:
        if np.isclose(success_rate, self._params.theta):
            return sigma
        if success_rate < self._params.theta:
            return sigma / self._params.alpha
        return sigma * self._params.alpha

    def _process_child(self, parent: Protein, sigma: float):
        child = self._mutate_protein(parent, sigma)
        return child, child.fitness < parent.fitness

    def run(self, sequence: str, callback: Callable[[int, Protein, float, bool], None] = None, callback_frequency: int = 1) -> Protein:
        generation:     int           = 0
        successes:      float         = 0.0
        sigma:          float         = self._params.sigma
        population:     List[Protein] = self._create_initial_population(sequence)
        best_offspring: Protein       = min(population, key=lambda p: p.fitness)

        with get_context("spawn").Pool() as pool:
            while generation < self._params.generations:
                generation += 1
                children: List[Protein] = population if self._params.plus_selection else []
                parents_to_mutate = [population[np.random.randint(self._params.population_size)] for _ in range(self._params.children_size)]

                results = pool.starmap(self._process_child, [(parent, sigma) for parent in parents_to_mutate])
                successes += sum(success for _, success in results)
                children.extend(child for child, _ in results)

                population = self._make_selection(children)

                if generation % self._params.mod_frequency == 0:
                    sigma = self._adaptive_adaption(sigma, successes / (self._params.mod_frequency * self._params.children_size))
                    successes = 0

                best_offspring = min(population, key=lambda p: p.fitness)

                is_premature_termination = self._is_premature_termination(best_offspring)

                if callback is not None and generation % callback_frequency == 0:
                    callback(generation, best_offspring, sigma, is_premature_termination or (generation >= self._params.generations))
                if is_premature_termination:
                    return best_offspring

        return best_offspring

    def _is_premature_termination(self, best_offspring: Protein) -> bool:
        return self._params.premature_termination is not None and self._reached_premature_termination(best_offspring)

    def _reached_premature_termination(self, best_offspring: Protein) -> bool:
        if self._previous_best is not None:
            if best_offspring.fitness == self._previous_best.fitness:
                self._previous_best_count += 1
            else:
                self._previous_best_count = 1
            if self._previous_best_count == self._params.premature_termination:
                return True
        self._previous_best = best_offspring
        return False

if __name__ == '__main__':
    import cProfile
    cProfile.run("EvolutionStrategy(EvolutionStrategyParams()).run('AAAA')")