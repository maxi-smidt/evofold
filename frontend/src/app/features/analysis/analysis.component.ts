import {Component} from '@angular/core';
import {NgxStructureViewerComponent, Settings, Source} from 'ngx-structure-viewer';
import {FileSelectEvent, FileUpload} from 'primeng/fileupload';
import {ButtonModule} from 'primeng/button';
import {AnalysisService} from '../../services/analysis.service';
import {AnalysisData} from '../../types/AnalysisData';

@Component({
  selector: 'app-analysis',
  imports: [
    NgxStructureViewerComponent,
    FileUpload,
    ButtonModule,
  ],
  templateUrl: './analysis.component.html',
  styleUrl: './analysis.component.css'
})
export class AnalysisComponent {
  protected source1: Source | undefined;
  protected source2: Source | undefined;
  protected settings: Partial<Settings> = {
    'background-color': '#2b3035ff',
    'backbone-color': '#6ea8fecc',
    'interaction-color': '#ff0000ff',
    'interaction-size': 1,
    prefer_label_asym_id: true,
  };
  protected sequence: string | undefined;
  protected angles: string | undefined;
  protected file: File | undefined;

  constructor(private analysisService: AnalysisService) {
  }

  protected async onSelect(event: FileSelectEvent) {
    this.file = event.currentFiles[0];
    this.source1 = {
      type: 'local' as const,
      format: 'mmcif' as const,
      label: 'EvoFold',
      binary: false,
      data: await this.file.text()
    };
  }

  protected async onUpload() {
    this.analysisService.analyze(await this.file!.text()).subscribe({
      next: (data: AnalysisData) => {
        this.source2 = {
          type: 'local' as const,
          format: 'mmcif' as const,
          label: 'EvoFold',
          binary: false,
          data: data.cifFile
        };
        this.sequence = data.sequence;
        this.angles = data.angles;
      }
    });
  }

  protected onAlign() {

  }
}
