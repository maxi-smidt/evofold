import {Component} from '@angular/core';
import {NgxStructureViewerComponent, Settings, Source} from 'ngx-structure-viewer';
import {FileSelectEvent, FileUpload} from 'primeng/fileupload';
import {ButtonModule} from 'primeng/button';
import {AnalysisService} from '../../services/analysis.service';
import {AnalysisData} from '../../types/AnalysisData';
import {AlignmentData} from '../../types/AlignmentData';
import {Tooltip} from 'primeng/tooltip';
import {Router} from '@angular/router';

@Component({
  selector: 'app-analysis',
  imports: [
    NgxStructureViewerComponent,
    FileUpload,
    ButtonModule,
    Tooltip,
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

  constructor(private analysisService: AnalysisService, private router: Router) {
  }

  protected async onSelect(event: FileSelectEvent, isSecond: boolean = false): Promise<void> {
    this.file = event.currentFiles[0];
    const sourceData = {
      type: 'local' as const,
      format: 'mmcif' as const,
      label: 'EvoFold',
      binary: false,
      data: await this.file.text()
    };

    if (!isSecond) {
      this.source1 = sourceData;
    } else {
      this.source2 = sourceData;
    }
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
    this.analysisService.align((this.source1 as any).data, (this.source2 as any).data).subscribe({
      next: (data: AlignmentData) => {
        this.router.navigate(['alignment'], {state: {alignmentData: data}}).then();
      }
    })
  }
}
