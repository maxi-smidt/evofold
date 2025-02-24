import {Component, OnDestroy, OnInit} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';
import {SimulationService} from '../../services/simulation.service';
import {Subject, Subscription} from 'rxjs';
import {NgxStructureViewerComponent, Settings, Source} from 'ngx-structure-viewer';
import {ResultEntries, SimulationData} from '../../types/SimulationData';
import {LocalStorageService} from '../../services/local-storage.service';
import {v4 as uuidv4} from 'uuid';
import {environment} from '../../../environments/environment';
import {Button} from 'primeng/button';

@Component({
  selector: 'app-structure-viewer',
  imports: [
    NgxStructureViewerComponent,
    Button,
  ],
  templateUrl: './structure-viewer.component.html',
  styleUrl: './structure-viewer.component.css'
})
export class StructureViewerComponent implements OnInit, OnDestroy {
  sequence: string = '';
  private subject: Subject<MessageEvent> | undefined;
  private subscription: Subscription | undefined;
  protected settings: Partial<Settings> = {
    'background-color': '#2b3035ff',
    'backbone-color': '#6ea8fecc',
    'interaction-color': '#ff0000ff',
    'interaction-size': 1,
    prefer_label_asym_id: true,
  };
  protected source: Source | undefined;
  protected results: ResultEntries = {};
  protected currentResultEntry: SimulationData | undefined;

  constructor(private route: ActivatedRoute,
              private router: Router,
              private simulationService: SimulationService,
              private localStorageService: LocalStorageService) {
  }

  ngOnInit() {
    const seq = this.route.snapshot.paramMap.get('sequence');
    this.localStorageService.clearAll();
    if (!seq) {
      this.router.navigate(['/']).then();
      return;
    }

    this.subject = this.simulationService.connect(environment.apiUrl);
    this.subscription = this.subject.subscribe({
      next: (msg: MessageEvent) => this.handleMessage(msg),
      error: (err: any) => console.error('WebSocket error:', err),
      complete: () => console.log('WebSocket connection closed')
    });
    this.sequence = seq;
    this.subject.next(this.sequence as unknown as MessageEvent);
  }

  ngOnDestroy() {
    if (this.subscription) {
      this.subscription.unsubscribe();
    }
    this.simulationService.disconnect();
  }

  handleMessage(msg: MessageEvent): void {
    const data = JSON.parse(msg.data) as SimulationData;
    const key = uuidv4();
    this.localStorageService.

    set(key, data.cifFile);
    this.results[key] = data;

    if (this.source === undefined) {
      this.onEntryClick(key);
    }
  }

  onEntryClick(localStorageKey: string) {
    const cif = this.localStorageService.get(localStorageKey);
    if (cif === null) return;

    this.source = {
      type: 'local' as const,
      format: 'mmcif' as const,
      label: 'EvoFold',
      binary: false,
      data: JSON.parse(cif)
    };

    this.currentResultEntry = this.results[localStorageKey];
  }

  protected readonly Object = Object;

  onDownloadClick(cif: string) {
    const newBlob = new Blob([cif], {type: "text/plain"});
    const data = window.URL.createObjectURL(newBlob);
    const link = document.createElement("a");
    link.href = data;
    link.download = `EvoFold_${new Date().toISOString()}.cif`;
    link.click();
  }
}
