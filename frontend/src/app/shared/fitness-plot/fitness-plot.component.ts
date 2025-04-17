import {Component, computed, input, InputSignal, signal, Signal} from '@angular/core';
import {PlotlyModule} from 'angular-plotly.js';
import {StoredSimulationData} from '../../types/StoredSimulationData';
import * as PlotlyJS from 'plotly.js-dist-min';
import {Data, Layout} from 'plotly.js-dist-min';
import {FormsModule} from '@angular/forms';
import {ToggleSwitch} from 'primeng/toggleswitch';

PlotlyModule.plotlyjs = PlotlyJS;

@Component({
  selector: 'app-fitness-plot',
  imports: [
    PlotlyModule,
    FormsModule,
    ToggleSwitch,
  ],
  templateUrl: './fitness-plot.component.html',
})
export class FitnessPlotComponent {
  protected layout: Partial<Layout> = {
    autosize: false,
    xaxis: {
      title: {
        text: 'Generation'
      },
      fixedrange: true,
    },
    yaxis: {
      title: {
        text: 'Fitness kJ/mol'
      },
      fixedrange: true,
    },
    yaxis2: {
      title: {
        text: 'Sigma'
      },
      fixedrange: true,
      overlaying: 'y',
      side: 'right'
    },
    modebar: {
      remove: ['select2d', 'lasso2d']
    },
    legend: {
      "orientation": "h",
    },
    height: 500,
    margin: {
      t: 40,
      r: 40,
      b: 40,
    },
  };
  public simulationData: InputSignal<[string, StoredSimulationData][]> = input.required();
  protected _showSigma = signal(false);
  protected data: Signal<Data[]> = computed(() => {
    const fitnessData: Data[] = [
      {
        x: this.simulationData().map(d => d[1].generation),
        y: this.simulationData().map(d => d[1].fitness),
        yaxis: 'y1',
        type: 'scatter',
        name: 'fitness'
      }
    ];
    if (this._showSigma()) {
      fitnessData.push(
        {
          x: this.simulationData().map(d => d[1].generation),
          y: this.simulationData().map(d => d[1].sigma),
          yaxis: 'y2',
          type: 'scatter',
          name: 'sigma'
        }
      );
    }
    return fitnessData;
  });

  protected get showSigma(): boolean {
    return this._showSigma();
  }

  protected set showSigma(value: boolean) {
    this._showSigma.set(value);
  }
}
