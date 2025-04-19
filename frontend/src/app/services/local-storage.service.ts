import {Injectable} from '@angular/core';
import * as LZString from 'lz-string';

@Injectable({
  providedIn: 'root'
})
export class LocalStorageService {

  constructor() {
  }

  public set(key: string, value: any) {
    const jsonString = JSON.stringify(value);
    const compressed = LZString.compress(jsonString);
    localStorage.setItem(key, compressed);
  }

  public get(key: string): any | null {
    const compressed = localStorage.getItem(key);
    if (compressed) {
      const decompressed = LZString.decompress(compressed);
      return decompressed ? JSON.parse(decompressed) : null;
    }
    return null;
  }

  public delete(key: string) {
    localStorage.removeItem(key);
  }

  public clearAll() {
    localStorage.clear();
  }

  public has(key: string): boolean {
    return localStorage.getItem(key) !== null;
  }
}
