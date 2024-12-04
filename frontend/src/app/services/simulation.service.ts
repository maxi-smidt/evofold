import {Injectable} from '@angular/core';
import {Observable, Subject} from 'rxjs';

@Injectable({
  providedIn: 'root'
})
export class SimulationService {
  private subject: Subject<MessageEvent> | undefined;
  private ws: WebSocket | undefined;

  constructor() {
  }

  public connect(url: string): Subject<MessageEvent> {
    if (!this.subject) {
      this.subject = this.create(url);
    }
    return this.subject;
  }

  private create(url: string): Subject<MessageEvent> {
    const subject = new Subject<MessageEvent>();

    this.ws = new WebSocket(url);

    const observable = new Observable<MessageEvent>(observer => {
      this.ws!.onmessage = (event) => observer.next(event);
      this.ws!.onerror = (error) => observer.error(error);
      this.ws!.onclose = () => observer.complete();
      return () => this.ws!.close();
    });

    observable.subscribe(subject);

    const sendToWebSocket = (data: any) => {
      let attempts = 0;
      const maxAttempts = 10;
      const interval = 10; // milliseconds

      const trySend = () => {
        if (this.ws!.readyState === WebSocket.OPEN) {
          this.ws!.send(JSON.stringify(data));
        } else if (attempts < maxAttempts) {
          attempts++;
          console.log(`Retrying (${attempts}/${maxAttempts})...`);
          setTimeout(trySend, interval);
        } else {
          console.error("Failed to send data: WebSocket is not open.");
        }
      };

      trySend();
    };

    subject.subscribe({
      next: sendToWebSocket
    });

    return subject;
  }

  disconnect(): void {
    if (this.ws) {
      this.ws.close();
      this.ws = undefined;
      this.subject = undefined;
    }
  }
}
