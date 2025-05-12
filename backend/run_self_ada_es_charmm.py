import optuna

from backend.algorithms.params.self_adaptive_es_params import SelfAdaptiveESParams
from backend.algorithms.self_adaptive_es import SelfAdaptiveES


def run_adaptive_es(sequence: str, params: SelfAdaptiveESParams) -> float:
    es = SelfAdaptiveES(params)
    best_protein = es.run(sequence)
    return best_protein.fitness  # assuming lower fitness is better

def objective(trial: optuna.trial.Trial) -> float:
    # 1. Suggest hyperparameters
    sigma = trial.suggest_float('sigma', 0.01, 360, log=True)
    population_size = trial.suggest_int('population_size', 1, 300)
    children_size = trial.suggest_int('children_size', population_size, 1000)
    plus_selection = trial.suggest_categorical('plus_selection', [True, False])
    strategy_param = trial.suggest_categorical('strategy_param', ['gene-wise', 'genome-wise'])

    # 2. Assemble AdaptiveESParams
    params = SelfAdaptiveESParams(
        generations=150,  # fixed
        population_size=population_size,
        children_size=children_size,
        plus_selection=plus_selection,
        force_field='charmm',
        sigma=sigma,
        strategy_param=strategy_param,
        premature_strategy='terminate',
        premature_stagnation=30,
        premature_sigma=0
    )

    # 3. Run the AdaptiveES
    sequence = "PEPTIDE"  # or any sequence you want to optimize on
    fitness = run_adaptive_es(sequence, params)

    # 4. Report the objective
    return -fitness  # maximize (if lower fitness is better)


if __name__ == "__main__":
    # 5. Create and run Optuna study
    print(f"Running self adaptive with charmm:")
    study = optuna.create_study(direction="maximize")  # because we minimize fitness, we maximize (-fitness)
    study.optimize(objective, n_trials=1000, n_jobs=-1)

    # 6. Best hyperparameters
    print("Best trial:")
    print(study.best_trial.params)
    print(study.best_trial.value)