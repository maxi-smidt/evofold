import {Component, effect, OnDestroy, OnInit, signal, WritableSignal} from '@angular/core';
import {ActivatedRoute, Router} from '@angular/router';
import {SimulationService} from '../../services/simulation.service';
import {Subscription} from 'rxjs';
import {NgxStructureViewerComponent, Settings, Source} from 'ngx-structure-viewer';
import {ReceivedSimulationData, StoredSimulationData} from '../../types/StoredSimulationData';
import {LocalStorageService} from '../../services/local-storage.service';
import {v4 as uuidv4} from 'uuid';
import {Button} from 'primeng/button';
import {AdaptiveEvolutionStrategyParams, EvolutionStrategyParams} from '../../types/EvolutionStrategyParams';
import {ProgressSpinner} from 'primeng/progressspinner';
import {RamachandranPlotComponent} from '../../shared/ramachandran-plot/ramachandran-plot.component';
import {Dialog, DialogModule} from 'primeng/dialog';
import {Tooltip} from 'primeng/tooltip';
import {FitnessPlotComponent} from '../../shared/fitness-plot/fitness-plot.component';
import {Slider} from 'primeng/slider';
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
    DialogModule,
  ],
  templateUrl: './structure-viewer.component.html',
  styleUrl: './structure-viewer.component.css'
})
export class StructureViewerComponent implements OnInit, OnDestroy {
  protected receivedLast: boolean = false;
  protected isLoading: boolean = false;
  protected sequence: string = '';
  protected isPlotsVisible: boolean = false;
  protected generationRange: WritableSignal<[number, number]> = signal([1, 1]);

  private subscription: Subscription | undefined;
  protected settings: Partial<Settings> = {
    'background-color': '#2b3035ff',
    'backbone-color': '#6ea8fecc',
    'interaction-color': '#ff0000ff',
    'interaction-size': 1,
    prefer_label_asym_id: true,
  };
  protected source: Source | undefined;
  protected results: [string, StoredSimulationData][] = [];
  protected selectedResults: [string, StoredSimulationData][] = [];
  protected currentResultEntry: StoredSimulationData | undefined;
  protected csvHeader: string = '';

  constructor(private route: ActivatedRoute,
              private router: Router,
              private simulationService: SimulationService,
              private localStorageService: LocalStorageService) {
    effect(() => {
      this.selectedResults = this.results.slice(this.generationRange()[0] - 1, this.generationRange()[1]);
      console.log(this.selectedResults);
      console.log(this.generationRange());
    });
  }

  public ngOnInit() {
    this.isLoading = true;
    const sequence: string = history.state.sequence;
    const method: string = history.state.method;
    const params: EvolutionStrategyParams | AdaptiveEvolutionStrategyParams = history.state.params;
    this.localStorageService.clearAll();

    this.csvHeader = `sep=;\nsequence: ${sequence}\nmethod: ${method}\nparams: ${JSON.stringify(params)}\n`;

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

  public ngOnDestroy() {
    if (this.subscription) {
      this.subscription.unsubscribe();
    }
    this.localStorageService.clearAll();
  }

  protected handleMessage(simulationData: ReceivedSimulationData): void {
    this.isLoading = false;
    const key = uuidv4();
    this.localStorageService.set(key, simulationData.atomPositions);
    this.results = [...this.results, [key, simulationData]];
    this.receivedLast = simulationData.isLast;
    this.generationRange.set([1, simulationData.generation]);
    if (this.source === undefined) {
      this.onEntryClick(0, simulationData.sequence);
    }
  }

  protected onEntryClick(index: number, sequence: string): void {
    const localStorageKey = this.results[index][0];
    const atomPositions = this.localStorageService.get(localStorageKey);
    if (atomPositions === null) return;
    this.source = {
      type: 'local' as const,
      format: 'mmcif' as const,
      label: 'EvoFold',
      binary: false,
      data: this.simulationService.convertAtomPositionsToCif(atomPositions, sequence)
    };

    this.currentResultEntry = this.results[index][1];
  }

  protected onDownloadCifClick() {
    this.download(`EvoFold_${new Date().toISOString()}.cif`, (this.source as any).data);
  }

  protected onDownloadCsvClick() {
    let csv = this.csvHeader;
    csv += 'generation;fitness;sigma;angles;cif\n';
    for (const result of this.results) {
      const r = result[1];
      const cif = this.simulationService.convertAtomPositionsToCif(this.localStorageService.get(result[0]), r.sequence);
      csv += `${r.generation};${r.fitness};${r.sigma};${JSON.stringify(r.angles)};"${cif}"\n`;
    }
    this.download(`EvoFold_${new Date().toISOString()}.csv`, csv);
  }

  private download(filename: string, text: string) {
    const newBlob = new Blob([text], {type: "text/plain"});
    const data = window.URL.createObjectURL(newBlob);
    const link = document.createElement("a");
    link.href = data;
    link.download = filename;
    link.click();
  }
}
