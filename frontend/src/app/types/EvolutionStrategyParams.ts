export interface EvolutionStrategyParams {
  generations: number;
  populationSize: number;
  childrenSize: number;
  plusSelection: boolean;
  forceField: 'charmm' | 'amber';
  sigma: number;
  prematureStrategy: 'none' | 'terminate' | 'restart';
  prematureStagnation: null | number;
  prematureSigma: null | number;
  prematureFitness: null | number;
}

export interface AdaptiveEvolutionStrategyParams {
  theta: number;
  alpha: number;
  modFrequency: number;
}

export interface DerandomizedEvolutionStrategyParams {
  alpha: number;
  tau: number;
}
