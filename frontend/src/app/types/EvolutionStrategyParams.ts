export interface EvolutionStrategyParams {
  generations: number;
  populationSize: number;
  childrenSize: number;
  plusSelection: boolean;
  forceField: 'charmm' | 'amber';
  sigma: number;
  theta: number;
  alpha: number;
  prematureTermination: null | number;
}
