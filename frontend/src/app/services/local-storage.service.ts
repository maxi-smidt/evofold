import {Injectable} from '@angular/core';
import * as LZString from 'lz-string';

@Injectable({
  providedIn: 'root'
})
export class LocalStorageService {

  constructor() {
  }

  set(key: string, value: any) {
    const jsonString = JSON.stringify(value);
    const compressed = LZString.compress(jsonString);
    localStorage.setItem(key, compressed);
  }

  get(key: string): any | null {
    const compressed = localStorage.getItem(key);
    if (compressed) {
      const decompressed = LZString.decompress(compressed);
      return decompressed ? JSON.parse(decompressed) : null;
    }
    return null;
  }

  clearAll() {
    localStorage.clear();
  }
}
