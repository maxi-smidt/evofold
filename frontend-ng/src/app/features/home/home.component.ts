import {Component} from '@angular/core';
import {Button} from 'primeng/button';
import {FormsModule} from '@angular/forms';
import {InputTextareaModule} from 'primeng/inputtextarea';
import {NgClass} from '@angular/common';
import {MessageModule} from 'primeng/message';
import {Router} from '@angular/router';

@Component({
  selector: 'app-home',
  imports: [
    Button,
    FormsModule,
    InputTextareaModule,
    NgClass,
    MessageModule,
  ],
  templateUrl: './home.component.html',
  styleUrl: './home.component.css'
})
export class HomeComponent {
  aminoAcids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];
  sequence: string = '';
  hasError = false;
  errorMessage = "This sequence should only contain 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', " +
    "'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'."

  constructor(private router: Router) {
  }


  onInput($event: Event) {
    const inputElement = $event.target as HTMLTextAreaElement;
    this.sequence = inputElement.value.toUpperCase();
    inputElement.value = this.sequence;

    this.hasError = [...this.sequence].some(char => !this.aminoAcids.includes(char));
  }

  onSimulateClick() {
    this.router.navigate(['structure', this.sequence]).then();
  }
}
