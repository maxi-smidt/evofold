import { Injectable } from '@angular/core';
import {HttpClient} from '@angular/common/http';
import {environment} from '../../environments/environment';
import {AnalysisData} from '../types/AnalysisData';
import {Observable} from 'rxjs';
import {AlignmentData} from '../types/AlignmentData';

@Injectable({
  providedIn: 'root'
})
export class AnalysisService {
  constructor(private http: HttpClient) { }

  public analyze(cifString: string): Observable<AnalysisData> {
    return this.http.post<AnalysisData>(`${environment.apiUrl}/analyze`, { cifFile1: cifString });
  }

  public align(cifString1: string, cifString2: string): Observable<AlignmentData> {
    return this.http.post<AlignmentData>(`${environment.apiUrl}/align`, { cifFile1: cifString1, cifFile2: cifString2 });
  }
}
