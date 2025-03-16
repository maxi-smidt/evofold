import { Injectable } from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {environment} from '../../environments/environment';
import {AnalysisData} from '../types/AnalysisData';
import {Observable} from 'rxjs';

@Injectable({
  providedIn: 'root'
})
export class AnalysisService {
  constructor(private http: HttpClient) { }

  public analyze(cifString: string): Observable<AnalysisData> {
    return this.http.post<AnalysisData>(`${environment.apiUrl}/analyze`, {cifFile: cifString});
  }
}
