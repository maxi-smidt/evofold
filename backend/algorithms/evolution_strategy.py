import numpy as np

from heapq import nsmallest
from typing import List, Callable

from backend.algorithms.evolution_strategy_params import EvolutionStrategyParams
from backend.structure.protein import Protein, AngleList


class EvolutionStrategy:
    def __init__(self, params: EvolutionStrategyParams):
        self._params = params
        self._previous_best = None
        self._previous_best_count = 0

    def _make_selection(self, children: List[Protein]) -> List[Protein]:
        fitness: Callable[[Protein], float] = lambda p: p.fitness
        return nsmallest(self._params.population_size, children, fitness)

    def _create_initial_population(self, sequence) -> List[Protein]:
        return [Protein(sequence) for _ in range(self._params.population_size)]

    @staticmethod
    def _mutate_protein(protein: Protein, sigma: float) -> Protein:
        return Protein(protein.sequence, EvolutionStrategy._gaussian_mutation(protein.angles, sigma))

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

    @staticmethod
    def process_child(args):
        parent, sigma = args
        child = EvolutionStrategy._mutate_protein(parent, sigma)
        return child, child.fitness < parent.fitness

    def run(self, sequence: str, callback: Callable[[int, Protein, float], None] = None, callback_frequency: int = 1) -> Protein:
        generation: int = 0
        s = 0.0
        k = 2
        sigma: float = Protein.ANGLE_MAX * 0.1
        population: List[Protein] = self._create_initial_population(sequence)

        while generation < self._params.generations:
            generation += 1
            children: List[Protein] = population if self._params.plus_selection else []

            for _ in range(self._params.children_size):
                parent = population[np.random.randint(self._params.population_size)]
                child = self._mutate_protein(parent, sigma)

                if child.fitness < parent.fitness:
                    s += 1

                children.append(child)

            population = self._make_selection(children)

            if generation % k == 0:
                sigma = self._adaptive_adaption(sigma, s / (k * self._params.children_size))
                s = 0

            best_offspring = min(population, key=lambda p: p.fitness)

            if self._params.premature_termination is not None and self._reached_premature_termination(best_offspring):
                return best_offspring

            if callback is not None and generation % callback_frequency == 0:
                callback(generation, best_offspring, sigma)

        return min(population, key=lambda p: p.fitness)

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
