export interface EvolutionStrategyParams {
  generations: number;
  populationSize: number;
  childrenSize: number;
  plusSelection: boolean;
  theta: number;
  alpha: number;
  prematureTermination: null | number;
}
