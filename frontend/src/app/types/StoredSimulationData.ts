export interface ReceivedSimulationData extends StoredSimulationData {
  atomPositions: [number, string, [number, string, number, number, number][]][];
}

export interface StoredSimulationData {
  generation: number;
  fitness: number;
  sequence: string;
  sigma: number;
  angles: [number, number][];
  isLast: boolean;
}

export interface ResultEntries {
  [localStorageKey: string]: StoredSimulationData;
}
