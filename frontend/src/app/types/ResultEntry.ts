export interface ResultEntries {
  [localStorageKey: string]: ResultEntry;
}

export interface ResultEntry {
  generation: number;
  fitness: string;
  sequence: string;
}
