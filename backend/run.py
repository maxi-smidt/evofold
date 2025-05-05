import time
from functools import partial

import numpy as np

from backend.algorithms.params.self_adaptive_es_params import SelfAdaptiveESParams
from backend.algorithms.self_adaptive_es import SelfAdaptiveES
from backend.structure.protein import Protein

AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

t = time.time()

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
        print(f"{method};{len(sequence)};{np.mean(run_times)}")
        with open("log.txt", "a") as f:
            f.write(f"{method};{len(sequence)};{np.mean(run_times)}\n")

def run_measurement(sequence: str):
    print(f"seq: {sequence}")
    print(f"starting amber n={len(sequence)}")
    measure_es('amber', SelfAdaptiveESParams(), sequence)
    print(f"starting charmm n={len(sequence)}")
    measure_es('charmm', SelfAdaptiveESParams(force_field='charmm'), sequence)


def main():
    for n in [10, 100, 1000]:
        s = generate_sequence(n)
        run_measurement(s)



if __name__ == '__main__':
    main()

"""
seq: GVCNNGNYGF
starting amber n=10
amber;10;1.0839974474906922
amber;10;1.275071988105774
amber;10;1.2327506279945373
amber;10;1.281261808872223
amber;10;1.3382270622253418
amber;10;1.375825891494751
amber;10;1.3036393237113952
amber;10;1.3307637071609497
amber;10;1.4193973565101623
amber;10;1.4190781831741333
amber;10;1.3250740122795106
amber;10;2.274530487060547
amber;10;2.474404036998749
amber;10;1.73793776512146
amber;10;1.358688669204712
amber;10;1.423206751346588
amber;10;1.368365659713745
amber;10;1.3712623524665832
amber;10;1.310444815158844
amber;10;1.4669630217552185
starting charmm n=10
charmm;10;1.597320418357849
charmm;10;1.6318775272369386
charmm;10;1.4677457928657531
charmm;10;1.5169169640541076
charmm;10;1.565239109992981
charmm;10;1.640558431148529
charmm;10;1.5441135358810425
charmm;10;1.5207390809059143
charmm;10;1.5411707830429078
charmm;10;1.5330018758773805
charmm;10;1.43572402715683
charmm;10;1.4349070620536803
charmm;10;1.5778377676010131
charmm;10;1.4313464283943176
charmm;10;1.491726167201996
charmm;10;1.5287493777275085
charmm;10;1.5198726153373718
charmm;10;1.4468026781082153
charmm;10;1.499310553073883
charmm;10;1.4845602440834045
seq: IRLLYDAQWMILGQGNWELAAVYAPFKDIDTHRWMCMQAQLGNWDIIIGESCQAEMVDSSAVCSRNNINVHPVWIWVTRYWAINCGEQPNMPWNMCEDEF
starting amber n=100
amber;100;20.304173810482027
amber;100;18.925528120994567
amber;100;16.770832476615904
amber;100;16.993850195407866
amber;100;16.688202965259553
amber;100;16.922659816741945
amber;100;16.856764874458314
amber;100;16.864094879627228
amber;100;16.86622601032257
amber;100;16.987021453380585
amber;100;16.997552857398986
amber;100;16.93378183364868
amber;100;17.188699266910554
amber;100;17.148275270462037
amber;100;16.65438227891922
amber;100;16.882444179058076
amber;100;16.449246118068697
amber;100;16.646972377300262
amber;100;16.691976654529572
amber;100;16.81734754085541
starting charmm n=100
"""