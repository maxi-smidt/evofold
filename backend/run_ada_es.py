from functools import partial

import optuna

from backend.algorithms.adaptive_es import AdaptiveES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams

def run_adaptive_es(sequence: str, params: AdaptiveESParams) -> float:
    es = AdaptiveES(params)
    best_protein = es.run(sequence)
    return best_protein.fitness  # assuming lower fitness is better

def print_callback(study, trial):
    print(f"Trial {trial.number} finished.")
    print(f"  Value: {trial.value}")
    print(f"  Params: {trial.params}")
    print(f"  Intermediate values: {trial.intermediate_values}")

def objective(trial: optuna.trial.Trial, force_field) -> float:
    # 1. Suggest hyperparameters
    sigma = trial.suggest_float('sigma', 0.01, 360, log=True)
    theta = trial.suggest_float('theta', 0.1, 0.5)
    alpha = trial.suggest_float('alpha', 1.1, 1.5)
    population_size = trial.suggest_int('population_size', 1, 300)
    children_size = trial.suggest_int('children_size', population_size, 1000)
    mod_frequency = trial.suggest_int('mod_frequency', 1, 10)
    plus_selection = trial.suggest_categorical('plus_selection', [True, False])

    print(f"[Trial {trial.number}] Testing: sigma={sigma:.3f}, pop={population_size}, plus={plus_selection}")

    # 2. Assemble AdaptiveESParams
    params = AdaptiveESParams(
        generations=150,  # fixed
        population_size=population_size,
        children_size=children_size,
        plus_selection=plus_selection,
        force_field=force_field,
        sigma=sigma,
        theta=theta,
        alpha=alpha,
        mod_frequency=mod_frequency
    )

    # 3. Run the AdaptiveES
    sequence = "PEPTIDE"  # or any sequence you want to optimize on
    fitness = run_adaptive_es(sequence, params)

    print(f"[Trial {trial.number}] Final fitness: {fitness:.4f}")

    # 4. Report the objective
    return -fitness  # maximize (if lower fitness is better)



if __name__ == "__main__":
    # 5. Create and run Optuna study
    for ff in ['amber', 'charmm']:
        print(f"Running with force field: {ff}")
        objective_with_ff = partial(objective, force_field=ff)
        study = optuna.create_study(direction="maximize")  # because we minimize fitness, we maximize (-fitness)
        study.optimize(objective_with_ff, n_trials=1000, n_jobs=10)

        # 6. Best hyperparameters
        print("Best trial:")
        print(study.best_trial.params)
        print(study.best_trial.value)