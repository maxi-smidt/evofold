import numpy as np
import gc
from backend.algorithms import self_adaptive_es
from backend.algorithms.params.self_adaptive_es_params import SelfAdaptiveESParams
from backend.structure.protein import Protein

AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def log(msg: str):
    with open('log.txt', 'a') as f:
        f.write(msg + '\n')
    print(msg)

for n in [10, 50, 100, 200, 500, 1000, 1500]:
    seq = ''.join([np.random.choice(AAs) for _ in range(n)])
    log(f"[n={n}] {seq}")
    for ff in ['amber', 'charmm']:
        gc.collect()
        log(f"  {ff}:")
        Protein.TIME_STRUCTURE = []
        Protein.TIME_CIF = []
        Protein.TIME_FITNESS = []

        params = SelfAdaptiveESParams(generations=5, force_field=ff, population_size=100, children_size=600, sigma=36)
        es = self_adaptive_es.SelfAdaptiveES(params)
        es.run(seq)
        t_structure = np.array(Protein.TIME_STRUCTURE)
        t_cif = np.array(Protein.TIME_CIF)
        t_fitness = np.array(Protein.TIME_FITNESS)
        log(f"{np.median(t_structure)},{np.median(t_cif)},{np.median(t_fitness)}")
        gc.collect()