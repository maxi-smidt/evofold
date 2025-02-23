import time
import numpy as np

from heapq import nsmallest
from typing import List, Callable

from backend.algorithms.evolution_strategy_params import EvolutionStrategyParams
from backend.structure.protein import Protein, AngleList


class EvolutionStrategy:
    def __init__(self, params: EvolutionStrategyParams):
        self._params = params

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
            children: List[Protein] = population if self._params.plus_selection else []
            start_time = time.time()

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

            if callback is not None and generation % callback_frequency == 0:
                callback(generation, min(population, key=lambda p: p.fitness), sigma)

            generation += 1
        return min(population, key=lambda p: p.fitness)

def main():
    # p = EvolutionStrategy1().run("AAAA", lambda i, x: print(f";{x.fitness};{x.angles}"))
    p = Protein("AAA", [(120, 30, 0), (-140, 60, 0), (50, 120, 0)])
    print(f"fitness: {p.fitness}")
    """
    print(f"angles: {p.angles}")
    parser = PDB.MMCIFParser()
    structure = parser.get_structure("protein", "test_1902.cif")

    model = structure[0]

    for chain in model:
        polypeptides = PDB.PPBuilder().build_peptides(chain)

        for poly in polypeptides:
            phi_psi_angles = poly.get_phi_psi_list()  # List of (phi, psi) tuples
            for aa in phi_psi_angles:
                if aa[0] and aa[1]:
                    print(f"phi: {np.rad2deg(aa[0]):.2f}, psi: {np.rad2deg(aa[1]):.2f}")
                elif not aa[0]:
                    print(f"phi: None, psi: {np.rad2deg(aa[1]):.2f}")
                else:
                    print(f"phi: {np.rad2deg(aa[0]):.2f}, psi: None")

    # print(p.fitness)
    # print(p.to_cif())
    """

# profiler = Profiler()
# profiler.start()

if __name__ == "__main__":
    start_time = time.time()
    for _ in range(100):
        main()
    print(f'{(time.time() - start_time) * 1000:.2f}ms')

# profiler.stop()
# print(profiler.output_text(unicode=True, color=True))