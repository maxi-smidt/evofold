import {Component} from '@angular/core';
import {Button} from 'primeng/button';
import {FormsModule} from '@angular/forms';
import {TextareaModule} from 'primeng/textarea';
import {NgClass} from '@angular/common';
import {MessageModule} from 'primeng/message';
import {Router} from '@angular/router';
import {SliderModule} from 'primeng/slider';
import {AdaptiveEvolutionStrategyParams, EvolutionStrategyParams} from '../../types/EvolutionStrategyParams';
import {SelectButtonModule} from 'primeng/selectbutton';
import {InputNumberModule} from 'primeng/inputnumber';
import {InputSwitchModule} from 'primeng/inputswitch';
import {TooltipModule} from 'primeng/tooltip';
import {MessageService} from 'primeng/api';
import {ToastModule} from 'primeng/toast';

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
  protected esOptions: any[] = [{label: 'Adaptive ES', value: 'adaptive'}, {label: 'Self-adaptive ES', value: 'self-adaptive'}];
  protected selectedESOption: 'adaptive' | 'self-adaptive' = 'adaptive';

  protected params: AdaptiveEvolutionStrategyParams = {
    generations: 500,
    populationSize: 100,
    childrenSize: 600,
    plusSelection: true,
    forceField: 'amber',
    prematureTermination: 10,
    sigma: 36,
    theta: 0.2,
    alpha: 1.225,
    modFrequency: 2
  }

  constructor(private router: Router,
              private messageService: MessageService) {
  }


  onInput(event: Event) {
    const inputElement = event.target as HTMLTextAreaElement;
    this.sequence = inputElement.value.toUpperCase();
    inputElement.value = this.sequence;

    this.errorMessage = null;

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

  get prematureTermination(): boolean {
    return this.params.prematureTermination !== null;
  }

  set prematureTermination(value: boolean) {
    this.params.prematureTermination = value ? 2 : null;
  }

  private validateParams(): boolean {
    const errors: string[] = []

    if (this.params.populationSize > this.params.childrenSize) {
      errors.push("The population size must be smaller than the children size.");
    }

    if (this.params.alpha <= 1) {
      errors.push("The modification factor alpha must be greater than 1.");
    }

    if (this.params.modFrequency > this.params.generations) {
      errors.push("The modification frequency cannot exceed the number of generations.");
    }

    for (const error of errors) {
      this.messageService.add({severity: 'error', summary: 'Error', detail: error});
    }

    return errors.length === 0;
  }

  private correctParams(): EvolutionStrategyParams | AdaptiveEvolutionStrategyParams {
    if (this.selectedESOption === 'self-adaptive') {
      return {
        generations: this.params.generations,
        populationSize: this.params.populationSize,
        childrenSize: this.params.childrenSize,
        plusSelection: this.params.plusSelection,
        forceField: this.params.forceField,
        sigma: this.params.sigma,
        prematureTermination: this.params.prematureTermination,
      };
    }
    return this.params;
  }
}
