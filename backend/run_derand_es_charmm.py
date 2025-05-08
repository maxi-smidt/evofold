import optuna

from backend.algorithms.derandomized_es import DerandomizedES
from backend.algorithms.params.derandomized_es_params import DerandomizedESParams


def run_adaptive_es(sequence: str, params: DerandomizedESParams) -> float:
    es = DerandomizedES(params)
    best_protein = es.run(sequence)
    return best_protein.fitness  # assuming lower fitness is better

def objective(trial: optuna.trial.Trial) -> float:
    # 1. Suggest hyperparameters
    sigma = trial.suggest_float('sigma', 0.01, 360, log=True)
    population_size = trial.suggest_int('population_size', 1, 300)
    children_size = trial.suggest_int('children_size', population_size, 1000)
    plus_selection = trial.suggest_categorical('plus_selection', [True, False])

    crossover = trial.suggest_categorical('crossover', ['global-arithmetic', 'global-uniform'])
    tau = trial.suggest_float('tau', 1, 14)
    alpha = trial.suggest_float('alpha', 0.01, 1)

    print(f"[Trial {trial.number}] Testing: sigma={sigma:.3f}, pop={population_size}, plus={plus_selection}")

    # 2. Assemble AdaptiveESParams
    params = DerandomizedESParams(
        generations=150,  # fixed
        population_size=population_size,
        children_size=children_size,
        plus_selection=plus_selection,
        force_field='charmm',
        sigma=sigma,
        alpha=alpha,
        crossover=crossover,
        tau=tau,
    )

    # 3. Run the AdaptiveES
    sequence = "PEPTIDE"  # or any sequence you want to optimize on
    fitness = run_adaptive_es(sequence, params)

    print(f"[Trial {trial.number}] Final fitness: {fitness:.4f}")

    # 4. Report the objective
    return -fitness  # maximize (if lower fitness is better)


if __name__ == "__main__":
    # 5. Create and run Optuna study
    print(f"Running derandomized with charmm:")
    study = optuna.create_study(direction="maximize")  # because we minimize fitness, we maximize (-fitness)
    study.optimize(objective, n_trials=1000, n_jobs=-1)

    # 6. Best hyperparameters
    print("Best trial:")
    print(study.best_trial.params)
    print(study.best_trial.value)