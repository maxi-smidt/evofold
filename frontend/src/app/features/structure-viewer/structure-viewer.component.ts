import {Component, OnDestroy, OnInit} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';
import {SimulationService} from '../../services/simulation.service';
import {Subscription} from 'rxjs';
import {NgxStructureViewerComponent, Settings, Source} from 'ngx-structure-viewer';
import {ReceivedSimulationData, ResultEntries, StoredSimulationData} from '../../types/StoredSimulationData';
import {LocalStorageService} from '../../services/local-storage.service';
import {v4 as uuidv4} from 'uuid';
import {Button} from 'primeng/button';
import {EvolutionStrategyParams} from '../../types/EvolutionStrategyParams';
import {ProgressSpinner} from 'primeng/progressspinner';

@Component({
  selector: 'app-structure-viewer',
  imports: [
    NgxStructureViewerComponent,
    Button,
    ProgressSpinner,
  ],
  templateUrl: './structure-viewer.component.html',
  styleUrl: './structure-viewer.component.css'
})
export class StructureViewerComponent implements OnInit, OnDestroy {
  protected isLoading: boolean = false;
  protected sequence: string = '';
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
  protected currentResultEntry: StoredSimulationData | undefined;

  constructor(private route: ActivatedRoute,
              private router: Router,
              private simulationService: SimulationService,
              private localStorageService: LocalStorageService) {
  }

  ngOnInit() {
    this.isLoading = true;
    const sequence: string = history.state.sequence;
    const params: EvolutionStrategyParams = history.state.params;
    this.localStorageService.clearAll();

    this.subscription = this.simulationService.getMessages().subscribe({
      next: (msg: ReceivedSimulationData) => this.handleMessage(msg),
      error: (err: any) => {
        console.error('WebSocket error:', err);
        this.isLoading = false;
      },
      complete: () => console.log('WebSocket connection closed')
    });

    this.simulationService.sendMessage({params, sequence} as unknown as MessageEvent)
    this.sequence = sequence;
  }

  ngOnDestroy() {
    if (this.subscription) {
      this.subscription.unsubscribe();
    }
    this.localStorageService.clearAll();
  }

  handleMessage(simulationData: ReceivedSimulationData): void {
    this.isLoading = false;
    const key = uuidv4();
    this.localStorageService.set(key, simulationData.atomPositions);
    this.results[key] = simulationData;
    if (this.source === undefined) {
      this.onEntryClick(key, simulationData.sequence);
    }
  }

  onEntryClick(localStorageKey: string, sequence: string): void {
    const atomPositions = this.localStorageService.get(localStorageKey);
    if (atomPositions === null) return;

    this.source = {
      type: 'local' as const,
      format: 'mmcif' as const,
      label: 'EvoFold',
      binary: false,
      data: this.simulationService.convertAtomPositionsToCif(atomPositions, sequence)
    };

    this.currentResultEntry = this.results[localStorageKey];
  }

  protected readonly Object = Object;

  onDownloadClick() {
    const newBlob = new Blob([(this.source as any).data], {type: "text/plain"});
    const data = window.URL.createObjectURL(newBlob);
    const link = document.createElement("a");
    link.href = data;
    link.download = `EvoFold_${new Date().toISOString()}.cif`;
    link.click();
  }
}
