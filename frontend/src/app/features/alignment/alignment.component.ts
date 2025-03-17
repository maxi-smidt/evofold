import {Component, OnInit} from '@angular/core';
import {AlignmentData} from '../../types/AlignmentData';
import {Locus, NgxStructureViewerComponent, Settings, Source} from 'ngx-structure-viewer';
import {DecimalPipe} from '@angular/common';

@Component({
  selector: 'app-alignment',
  imports: [
    NgxStructureViewerComponent,
    DecimalPipe,
  ],
  templateUrl: './alignment.component.html',
  styleUrl: './alignment.component.css'
})
export class AlignmentComponent implements OnInit {
  protected settings: Partial<Settings> = {
    'background-color': '#2b3035ff',
    'backbone-color': '#6ea8fecc',
    'interaction-color': '#ff0000ff',
    'interaction-size': 1,
    prefer_label_asym_id: true,
  };
  protected loci: Locus[] = [
    {chain: 'A', color: '#ff0000ff'},
    {chain: 'B', color: '#0d6efd'},
  ];
  protected source: Source | undefined;
  protected alignmentData: AlignmentData | undefined;

  public ngOnInit() {
    this.alignmentData = history.state.alignmentData;

    this.source = {
      type: 'local' as const,
      format: 'mmcif' as const,
      label: 'EvoFold',
      binary: false,
      data: this.alignmentData!.cifFile
    }
  }
}
