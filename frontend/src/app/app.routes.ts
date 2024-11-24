import { Routes } from '@angular/router';
import { HomeComponent } from './features/home/home.component';
import {StructureViewerComponent} from './features/structure-viewer/structure-viewer.component';

export const routes: Routes = [
  { path: '', component: HomeComponent },
  { path: 'structure', component: StructureViewerComponent },
  { path: '**', redirectTo: '', pathMatch: 'full' },
];
