export interface EvolutionStrategyParams {
  generations: number;
  populationSize: number;
  childrenSize: number;
  plusSelection: boolean;
  forceField: 'charmm' | 'amber';
  sigma: number;
  prematureTermination: null | number;
}

export interface AdaptiveEvolutionStrategyParams extends EvolutionStrategyParams {
  theta: number;
  alpha: number;
  modFrequency: number;
}
