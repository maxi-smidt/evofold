import { Component } from '@angular/core';
import {HeaderComponent} from './header/header.component';
import {RouterOutlet} from '@angular/router';

@Component({
  selector: 'app-layout',
  imports: [
    HeaderComponent,
    RouterOutlet
  ],
  templateUrl: './layout.component.html'
})
export class LayoutComponent {
}
