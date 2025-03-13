export interface SimulationData {
  generation: number;
  fitness: number;
  sequence: string;
  cifFile: string;
  sigma: number;
}

export interface ResultEntries {
  [localStorageKey: string]: SimulationData;
}
