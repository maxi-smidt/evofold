import copy
import time

import numpy as np
from heapq import nsmallest
from typing import List, Callable

from openmm.unit import Quantity

from backend.structure.protein import Protein


class EvolutionStrategy1:
    def __init__(self):
        self.generations        = 500   # the amount of generations that will be produced
        self.population_size    = 100   # the size of the population
        self.children_size      = 600   # the amount of children produced by one population
        self.plus_selection     = False # if true (µ+λ) else (µ,λ)
        self.theta              = 0.2   # threshold for adaptive adaption (1 / 5 success rule)
        self.alpha              = 1.224 # modification factor

    def _make_selection(self, children: List[Protein]) -> List[Protein]:
        fitness: Callable[[Protein], Quantity] = lambda p: p.fitness()
        return nsmallest(self.population_size, children, fitness)

    def _create_initial_population(self, sequence) -> List[Protein]:
        base_protein = Protein(sequence)
        return [copy.deepcopy(base_protein) for _ in range(self.population_size)]

    def _mutate_protein(self, protein: Protein, sigma: float) -> Protein:
        ...

    def _adaptive_adaption(self, sigma, success_rate: float) -> float:
        if np.isclose(success_rate, self.theta):
            return sigma
        if success_rate < self.theta:
            return sigma / self.alpha
        return sigma * self.alpha

    def run(self, sequence: str, callback: Callable[[int, Protein], None] = None, callback_frequency: int = 5) -> Protein:
        generation: int = 0
        s = 0.0
        k = 0
        sigma: float = 0
        population: List[Protein] = self._create_initial_population(sequence)

        while generation < self.generations:
            children: List[Protein] = population if self.plus_selection else []

            for _ in range(self.children_size):
                parent = population[np.random.randint(self.population_size)]
                child = self._mutate_protein(parent, sigma)

                if child.fitness() > parent.fitness():
                    s += 1

                children.append(child)

            population = self._make_selection(children)

            if generation % k == 0:
                sigma = self._adaptive_adaption(sigma, s / (k * self.children_size))
                s = 0

            if callback is not None and generation % callback_frequency == 0:
                callback(generation, min(population, key=lambda p: p.fitness()))

        return min(population, key=lambda p: p.fitness())


def main():
    # EvolutionStrategy1().run("PEPTIDE")
    p = Protein("ARNDCQEGHILKMFPSTWYV")
    print(p.fitness())
    print(p.to_cif())

# current fitness for ARNDCQEGHILKMFPSTWYV: 9.79537141081154e+19 kJ/mol
# current fitness for ARNDCQEGHILKMFPSTWYV: 11734092.587478638 kJ/mol
# current fitness for ARNDCQEGHILKMFPSTWYV: 11562952.821296692 kJ/mol


if __name__ == "__main__":
    start_time = time.time()
    main()
    print(f'{(time.time() - start_time) * 1000:.2f}ms')
