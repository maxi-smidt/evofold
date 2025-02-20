import time
import numpy as np

from concurrent.futures import ProcessPoolExecutor
from heapq import nsmallest
from typing import List, Callable

from backend.structure.protein import Protein, AngleList

from Bio import PDB

class EvolutionStrategy1:
    def __init__(self):
        self.generations        = 500   # the amount of generations that will be produced
        self.population_size    = 100   # the size of the population
        self.children_size      = 600   # the amount of children produced by one population
        self.plus_selection     = False # if true (µ+λ) else (µ,λ)
        self.theta              = 0.2   # threshold for adaptive adaption (1 / 5 success rule)
        self.alpha              = 1.224 # modification factor

    def _make_selection(self, children: List[Protein]) -> List[Protein]:
        fitness: Callable[[Protein], float] = lambda p: p.fitness
        return nsmallest(self.population_size, children, fitness)

    def _create_initial_population(self, sequence) -> List[Protein]:
        return [Protein(sequence) for _ in range(self.population_size)]

    @staticmethod
    def _mutate_protein(protein: Protein, sigma: float) -> Protein:
        return Protein(protein.sequence, EvolutionStrategy1._gaussian_mutation(protein.angles, sigma))

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
        if np.isclose(success_rate, self.theta):
            return sigma
        if success_rate < self.theta:
            return sigma / self.alpha
        return sigma * self.alpha

    @staticmethod
    def process_child(args):
        parent, sigma = args
        child = EvolutionStrategy1._mutate_protein(parent, sigma)
        return child, child.fitness > parent.fitness

    def run(self, sequence: str, callback: Callable[[int, Protein], None] = None, callback_frequency: int = 5) -> Protein:
        generation: int = 0
        s = 0.0
        k = 2
        sigma: float = 0
        population: List[Protein] = self._create_initial_population(sequence)

        while generation < self.generations:
            children: List[Protein] = population if self.plus_selection else []
            start_time = time.time()

            for _ in range(self.children_size):
                parent = population[np.random.randint(self.population_size)]
                child = self._mutate_protein(parent, sigma)

                if child.fitness > parent.fitness:
                    s += 1

                children.append(child)

            population = self._make_selection(children)

            if generation % k == 0:
                sigma = self._adaptive_adaption(sigma, s / (k * self.children_size))
                s = 0

            if callback is not None and generation % callback_frequency == 0:
                callback(generation, min(population, key=lambda p: p.fitness))

            print(f'Generation {generation}: {(time.time() - start_time) * 1000:.2f}ms')
            generation += 1
        return min(population, key=lambda p: p.fitness)


def main():
    # p = EvolutionStrategy1().run("PEPTIDE")
    p = Protein("AAA", [(120, 30, 0), (-140, 60, 0), (50, 120, 0)])
    p.to_cif("test_1902.cif")
    print(f"fitness: {p.fitness}")
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


if __name__ == "__main__":
    for _ in range(1):
        start_time = time.time()
        main()
        print(f'{(time.time() - start_time) * 1000:.2f}ms')