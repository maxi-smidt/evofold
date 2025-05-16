import time
from functools import partial

import numpy as np

from backend.algorithms.params.self_adaptive_es_params import SelfAdaptiveESParams
from backend.algorithms.self_adaptive_es import SelfAdaptiveES
from backend.structure.protein import Protein

AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

t = time.time()

def log(msg: str):
    print(msg)

def generate_sequence(n: int):
    return ''.join([np.random.choice(AAs) for _ in range(n)])

def save_time(generation: int, protein: Protein, fitness: float, is_best: bool, l: list):
    global t
    l.append(time.time() - t)
    t = time.time()

def measure_es(method: str, params: SelfAdaptiveESParams, sequence: str):
    global t
    for i in range(20):
        t = time.time()
        run_times = []
        cb = partial(save_time, l=run_times)
        es = SelfAdaptiveES(params)
        es.run(sequence, callback=cb)
        log(f"{method};{len(sequence)};{np.mean(run_times)}")

def run_measurement(sequence: str):
    log(f"seq: {sequence}")
    log(f"starting amber n={len(sequence)}")
    measure_es('amber', SelfAdaptiveESParams(), sequence)
    log(f"starting charmm n={len(sequence)}")
    measure_es('charmm', SelfAdaptiveESParams(force_field='charmm'), sequence)

def main():
    for n in [10, 30, 75, 100, 150]:
        s = generate_sequence(n)
        run_measurement(s)


if __name__ == '__main__':
    main()