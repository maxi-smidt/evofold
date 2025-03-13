import {Component} from '@angular/core';
import {Button} from 'primeng/button';
import {FormsModule} from '@angular/forms';
import {TextareaModule} from 'primeng/textarea';
import {NgClass} from '@angular/common';
import {MessageModule} from 'primeng/message';
import {Router} from '@angular/router';
import {SliderModule} from 'primeng/slider';
import {EvolutionStrategyParams} from '../../types/EvolutionStrategyParams';
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
  aminoAcids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];
  sequence: string = '';
  hasError = false;
  errorMessage = "This sequence should only contain 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', " +
    "'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'."
  stateOptions: any[] = [{label: '(µ+λ)-Selection', value: 'plus'}, {label: '(µ,λ)-Selection', value: 'comma'}];

  params: EvolutionStrategyParams = {
    generations: 500,
    populationSize: 100,
    childrenSize: 600,
    plusSelection: true,
    theta: 0.2,
    alpha: 1.225,
    prematureTermination: 10
  }

  constructor(private router: Router,
              private messageService: MessageService) {
  }


  onInput($event: Event) {
    const inputElement = $event.target as HTMLTextAreaElement;
    this.sequence = inputElement.value.toUpperCase();
    inputElement.value = this.sequence;

    this.hasError = [...this.sequence].some(char => !this.aminoAcids.includes(char));
  }

  onSimulateClick() {
    if (this.validateParams()) {
      this.router.navigate(['structure'],
        {
          state: {
            sequence: this.sequence,
            params: this.params
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

    if (errors.length > 0) {
      for (const error of errors) {
        this.messageService.add({severity: 'error', summary: 'Error', detail: error});
      }

      return false;
    }

    return true;
  }
}
