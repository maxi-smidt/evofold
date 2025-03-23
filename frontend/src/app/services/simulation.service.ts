import {Injectable} from '@angular/core';
import {Observable} from 'rxjs';
import {webSocket, WebSocketSubject} from 'rxjs/webSocket';
import {environment} from '../../environments/environment';

@Injectable({
  providedIn: 'root'
})
export class SimulationService {
  private socket: WebSocketSubject<any> | null = null;

  private createWebSocket() {
    if (!this.socket || this.socket.closed) {
      this.socket = webSocket({url: environment.wsUrl});
    }
  }

  public sendMessage(message: any) {
    this.socket!.next(message);
  }

  public getMessages(): Observable<any> {
    this.createWebSocket();
    return this.socket!.asObservable();
  }

  public closeConnection() {
    this.socket?.complete();
  }

  public convertAtomPositionsToCif(atomPositions: string[][], sequence: string): string {
    let cifStr = `data_EvoFold_${sequence}\n`;
    cifStr += '#\n';
    cifStr += this.cifPositionsToStr(atomPositions);
    return cifStr;
  }

  private cifPositionsToStr(atomPositions: string[][]): string {
    const atomHeader = [
      'loop_', '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol', '_atom_site.label_atom_id',
      '_atom_site.label_alt_id', '_atom_site.label_comp_id', '_atom_site.label_asym_id',
      '_atom_site.label_entity_id', '_atom_site.label_seq_id', '_atom_site.pdbx_PDB_ins_code',
      '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
      '_atom_site.B_iso_or_equiv', '_atom_site.auth_seq_id', '_atom_site.auth_asym_id',
      '_atom_site.pdbx_PDB_model_num'
    ];

    let positionsString = atomHeader.join('\n');
    for (const aa of atomPositions) {
      for (const atom of aa[2]) {
        const data = [
          'ATOM',       // group_PDB (here always ATOM)
          `${atom[0]}`, // id
          atom[1][0],   // type_symbol (e.g. N, O, C)
          atom[1],      // label_atom_id (e.g. CA, O, H1)
          '.',          // label_alt_id (always .)
          aa[1],        // label_comp_id
          'C',          // label_asym_id (always C)
          '3',          // label_entity_id (always 3)
          `${aa[0]}`,   // label_seq_id
          '?',          // pdbx_PDB_ins_code (always ?)
          `${atom[2]}`, // Cartn_x
          `${atom[3]}`, // Cartn_y
          `${atom[4]}`, // Cartn_z
          '1.00',       // occupancy (always 1.00)
          '0.00',       // B_iso_or_equiv (always 0.00)
          `${aa[0]}`,   // auth_seq_id (like label_seq_id)
          'A',          // auth_asym_id (always A)
          '1'           // pdbx_PDB_model_num (1 because we always have only one chain)
        ];
        positionsString += '\n' + data.join('\t');
      }
    }
    return positionsString;
  }
}
