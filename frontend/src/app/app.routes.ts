import {Routes} from '@angular/router';
import {HomeComponent} from './features/home/home.component';
import {StructureViewerComponent} from './features/structure-viewer/structure-viewer.component';
import {LayoutComponent} from './core/layout/layout.component';
import {AnalysisComponent} from './features/analysis/analysis.component';
import {AlignmentComponent} from './features/alignment/alignment.component';

export const routes: Routes = [
  {
    path: '', component: LayoutComponent, children: [
      { path: '', component: HomeComponent },
      { path: 'structure', component: StructureViewerComponent },
      { path: 'analysis', component: AnalysisComponent },
      { path: 'alignment', component: AlignmentComponent }
    ]
  },
  { path: '**', redirectTo: '', pathMatch: 'full' },
];
