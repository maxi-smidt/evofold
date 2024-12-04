import { Injectable } from '@angular/core';

@Injectable({
  providedIn: 'root'
})
export class LocalStorageService {

  constructor() { }

  set(key: string, value: any) {
    localStorage.setItem(key, JSON.stringify(value));
  }

  get(key: string): string | null {
    return localStorage.getItem(key);
  }

  clearAll() {
    localStorage.clear();
  }
}
