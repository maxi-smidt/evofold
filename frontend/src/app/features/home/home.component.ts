import {Component} from '@angular/core';
import {Button} from 'primeng/button';
import {FormsModule} from '@angular/forms';
import {TextareaModule} from 'primeng/textarea';
import {NgClass} from '@angular/common';
import {MessageModule} from 'primeng/message';
import {Router} from '@angular/router';
import {SliderModule} from 'primeng/slider';
import {AdaptiveEvolutionStrategyParams, DerandomizedEvolutionStrategyParams, EvolutionStrategyParams} from '../../types/EvolutionStrategyParams';
import {SelectButtonChangeEvent, SelectButtonModule} from 'primeng/selectbutton';
import {InputNumberModule} from 'primeng/inputnumber';
import {InputSwitchModule} from 'primeng/inputswitch';
import {TooltipModule} from 'primeng/tooltip';
import {MessageService} from 'primeng/api';
import {ToastModule} from 'primeng/toast';
import {Tab, TabList, TabPanel, TabPanels, Tabs} from 'primeng/tabs';
import {LocalStorageService} from '../../services/local-storage.service';

@Component({
  selector: 'app-home',
  imports: [
    Button,
    FormsModule,
    TextareaModule,
    NgClass,
    MessageModule,
    SliderModule,
    SelectButtonModule,
    InputNumberModule,
    InputSwitchModule,
    TooltipModule,
    ToastModule,
    Tabs,
    TabList,
    Tab,
    TabPanels,
    TabPanel,
  ],
  templateUrl: './home.component.html',
  styleUrl: './home.component.css',
  providers: [MessageService]
})
export class HomeComponent {
  protected aminoAcids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];
  protected sequence: string = '';
  protected errorMessage: string | null = null;
  protected selectionOptions: any[] = [{label: '(µ+λ)-Selection', value: 'plus'}, {label: '(µ,λ)-Selection', value: 'comma'}];
  protected forceFieldOptions: any[] = [{label: 'AMBER', value: 'amber'}, {label: 'CHARMM', value: 'charmm'}];
  protected esOptions: any[] = [{label: 'Adaptive ES', value: 'adaptive'}, {label: 'Self-adaptive ES', value: 'self-adaptive'}, {label: 'Derandomized ES', value: 'derandomized'}];
  protected selectedESOption: 'adaptive' | 'self-adaptive' | 'derandomized' = 'adaptive';
  protected pcOptions: any[] = [{label: 'None', value: 'none'}, {label: 'Termination', value: 'terminate'}, {label: 'Restart', value: 'restart'}];
  private previousPrematureStrategy: string = 'none';

  protected params: EvolutionStrategyParams = {
    generations: 500,
    populationSize: 100,
    childrenSize: 600,
    plusSelection: true,
    forceField: 'amber',
    prematureStrategy: 'none',
    prematureStagnation: null,
    prematureSigma: null,
    prematureFitness: null,
    sigma: 36,
  }

  protected derandomizedParams: DerandomizedEvolutionStrategyParams = {
    tau: 0,
    alpha: 1 / 0,
  }

  protected adaptiveParams: AdaptiveEvolutionStrategyParams = {
    theta: 0.2,
    alpha: 1.225,
    modFrequency: 2
  }

  constructor(private router: Router,
              private messageService: MessageService,
              private localStorageService: LocalStorageService) {
  }


  onInput(event: Event) {
    const inputElement = event.target as HTMLTextAreaElement;
    this.sequence = inputElement.value.toUpperCase();
    inputElement.value = this.sequence;

    this.errorMessage = null;

    const sqrtOfLength = Math.sqrt(this.sequence.length * 2);
    this.derandomizedParams.alpha = 1 / sqrtOfLength;
    this.derandomizedParams.tau = sqrtOfLength;

    if (this.sequence.length < 2) {
      this.errorMessage = "The sequence must be at least two amino acids."
    }

    if ([...this.sequence].some(char => !this.aminoAcids.includes(char))) {
      this.errorMessage = "This sequence should only contain 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', " +
        "'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'."
    }
  }

  onSimulateClick() {
    if (this.validateParams()) {
      this.router.navigate(['structure'],
        {
          state: {
            sequence: this.sequence,
            method: this.selectedESOption,
            params: this.correctParams()
          }
        }
      ).then();
    }
  }

  get plusSelection(): string {
    return this.params.plusSelection ? 'plus' : 'comma';
  }

  set plusSelection(value: string) {
    this.params.plusSelection = value === 'plus';
  }

  get prematureStagnation(): boolean {
    return this.params.prematureStagnation !== null;
  }

  set prematureStagnation(value: boolean) {
    this.params.prematureStagnation = value ? 15 : null;
  }

  get prematureSigma(): boolean {
    return this.params.prematureSigma !== null;
  }

  set prematureSigma(value: boolean) {
    this.params.prematureSigma = value ? 0.01 : null;
  }

  get prematureFitness(): boolean {
    return this.params.prematureFitness !== null;
  }

  set prematureFitness(value: boolean) {
    this.params.prematureFitness = value ? 2000 : null;
  }

  private validateParams(): boolean {
    const errors: string[] = []

    if (this.params.populationSize > this.params.childrenSize) {
      errors.push("The population size must be smaller than the children size.");
    }

    if (this.selectedESOption === 'adaptive' && this.adaptiveParams.alpha <= 1) {
      errors.push("The modification factor alpha must be greater than 1.");
    }

    if (this.selectedESOption === 'adaptive' && this.adaptiveParams.modFrequency > this.params.generations) {
      errors.push("The modification frequency cannot exceed the number of generations.");
    }

    for (const error of errors) {
      this.messageService.add({severity: 'error', summary: 'Error', detail: error});
    }

    return errors.length === 0;
  }

  private correctParams(): any {
    switch (this.selectedESOption) {
      case 'adaptive': return {...this.params, ...this.adaptiveParams};
      case 'self-adaptive': return this.params;
      case 'derandomized': return {...this.params, ...this.derandomizedParams};
      default: throw new Error('Unknown option');
    }
  }

  protected onPrematureStrategyChange(event: SelectButtonChangeEvent) {
    if (event.value === 'none') {
      this.params.prematureStagnation = null;
      this.params.prematureSigma = null;
      this.params.prematureFitness = null;
    } else {
      if (this.previousPrematureStrategy === 'none') {
        this.params.prematureStagnation = 15;
        this.params.prematureSigma = 0.01;
        this.params.prematureFitness = null;
      }
    }
    this.previousPrematureStrategy = event.value;
  }

  protected onSaveParamsClick() {
    this.localStorageService.set('defaultParams', this.params);
    this.localStorageService.set('method', this.selectedESOption);

    if (this.selectedESOption === 'adaptive') {
      this.localStorageService.set('additionalParams', this.adaptiveParams);
    }
    if (this.selectedESOption === 'derandomized') {
      this.localStorageService.set('additionalParams', this.derandomizedParams);
    }

    this.messageService.add({severity: 'success', summary: 'Success', detail: 'Parameters saved successfully.'});
  }

  protected onLoadParamsClick() {
    if (!this.localStorageService.has('defaultParams')) {
      this.messageService.add({severity: 'error', summary: 'Error', detail: 'No parameters have been saved yet.'});
      return;
    }
    this.params = this.localStorageService.get('defaultParams');
    this.selectedESOption = this.localStorageService.get('method');
    if (this.selectedESOption === 'adaptive') {
      this.adaptiveParams = this.localStorageService.get('additionalParams');
    }
    if (this.selectedESOption === 'derandomized') {
      this.derandomizedParams = this.localStorageService.get('additionalParams');
    }
  }
}
