import numpy as np

from backend.algorithms import self_adaptive_es
from backend.algorithms.params.self_adaptive_es_params import SelfAdaptiveESParams
from backend.structure.protein import Protein

AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


for n in [10, 50, 100, 200, 500, 1000, 1500]:
    seq = ''.join([np.random.choice(AAs) for _ in range(n)])
    print(f"[n={n}] {seq}")
    for ff in ['amber', 'charmm']:
        print(f"  {ff}:")
        Protein.TIME_STRUCTURE = []
        Protein.TIME_CIF = []
        Protein.TIME_FITNESS = []

        params = SelfAdaptiveESParams(generations=5, force_field=ff, population_size=100, children_size=600, sigma=36)
        es = self_adaptive_es.SelfAdaptiveES(params)
        es.run(seq)
        t_structure = np.array(Protein.TIME_STRUCTURE)
        t_cif = np.array(Protein.TIME_CIF)
        t_fitness = np.array(Protein.TIME_FITNESS)

        print("f,min,max,mean,median,std,samples")
        print(f"structure, {t_structure.min()}, {t_structure.max()}, {t_structure.mean()}, {np.median(t_structure)}, {t_structure.std()}, {t_structure.shape[0]}")
        print(f"cif, {t_cif.min()}, {t_cif.max()}, {t_cif.mean()}, {np.median(t_cif)}, {t_cif.std()}, {t_cif.shape[0]}")
        print(f"fitness, {t_fitness.min()}, {t_fitness.max()}, {t_fitness.mean()}, {np.median(t_fitness)}, {t_fitness.std()}, {t_fitness.shape[0]}")
        print()

#
#
#
#
#
#
#