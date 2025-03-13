import { Component } from '@angular/core';
import {Button} from 'primeng/button';
import {RouterLink} from '@angular/router';

@Component({
  selector: 'app-header',
  imports: [
    Button,
    RouterLink
  ],
  templateUrl: './header.component.html'
})
export class HeaderComponent {
  routeToGithub() {
    window.open('https://github.com/maxi-smidt/evofold', "_blank");
  }
}
