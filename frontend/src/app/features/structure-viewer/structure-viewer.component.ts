import {Component, OnDestroy, OnInit} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';
import {SimulationService} from '../../services/simulation.service';
import {Subscription} from 'rxjs';
import {NgxStructureViewerComponent, Settings, Source} from 'ngx-structure-viewer';
import {ReceivedSimulationData, ResultEntries, StoredSimulationData} from '../../types/StoredSimulationData';
import {LocalStorageService} from '../../services/local-storage.service';
import {v4 as uuidv4} from 'uuid';
import {Button} from 'primeng/button';
import {AdaptiveEvolutionStrategyParams, EvolutionStrategyParams} from '../../types/EvolutionStrategyParams';
import {ProgressSpinner} from 'primeng/progressspinner';
import {RamachandranPlotComponent} from '../../shared/ramachandran-plot/ramachandran-plot.component';
import {Dialog} from 'primeng/dialog';
import {Tooltip} from 'primeng/tooltip';
import {FitnessPlotComponent} from '../../shared/fitness-plot/fitness-plot.component';
import {Slider, SliderSlideEndEvent} from 'primeng/slider';
import {FormsModule} from '@angular/forms';

@Component({
  selector: 'app-structure-viewer',
  imports: [
    NgxStructureViewerComponent,
    Button,
    ProgressSpinner,
    RamachandranPlotComponent,
    Dialog,
    Tooltip,
    FitnessPlotComponent,
    Slider,
    FormsModule,
  ],
  templateUrl: './structure-viewer.component.html',
  styleUrl: './structure-viewer.component.css'
})
export class StructureViewerComponent implements OnInit, OnDestroy {
  protected receivedLast: boolean = false;
  protected isLoading: boolean = false;
  protected sequence: string = '';
  protected isRamachandranVisible: boolean = false;
  protected isSummaryVisible: boolean = false;
  protected generationRange: [number, number] = [0, 0];
  protected fixedGenerationRange: [number, number] = [0, 0];

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
  protected resultsArray: StoredSimulationData[] = [];
  protected selectedResultsArray: StoredSimulationData[] = [];
  protected currentResultEntry: StoredSimulationData | undefined;

  constructor(private route: ActivatedRoute,
              private router: Router,
              private simulationService: SimulationService,
              private localStorageService: LocalStorageService) {
  }

  get obj() {
    return Object;
  }

  ngOnInit() {
    this.isLoading = true;
    const sequence: string = history.state.sequence;
    const method: string = history.state.method;
    const params: EvolutionStrategyParams | AdaptiveEvolutionStrategyParams = history.state.params;
    this.localStorageService.clearAll();

    this.subscription = this.simulationService.getMessages().subscribe({
      next: (msg: ReceivedSimulationData) => this.handleMessage(msg),
      error: (err: any) => {
        console.error('WebSocket error:', err);
        this.isLoading = false;
      },
      complete: () => console.log('WebSocket connection closed')
    });

    this.simulationService.sendMessage({params, sequence, method} as unknown as MessageEvent)
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
    this.resultsArray = [...this.resultsArray, simulationData];
    if (simulationData.isLast) {
      this.receivedLast = simulationData.isLast;
      this.generationRange = [1, simulationData.generation];
      this.fixedGenerationRange = [1, simulationData.generation];
      this.selectedResultsArray = this.resultsArray;
    }
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

  protected onSlideEnd(event: SliderSlideEndEvent) {
    if (event.values) {
      this.selectedResultsArray = this.resultsArray.slice(event.values[0] - 1, event.values[1]);
    }
  }

  protected onDownloadClick() {
    const newBlob = new Blob([(this.source as any).data], {type: "text/plain"});
    const data = window.URL.createObjectURL(newBlob);
    const link = document.createElement("a");
    link.href = data;
    link.download = `EvoFold_${new Date().toISOString()}.cif`;
    link.click();
  }
}
