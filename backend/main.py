import asyncio
import json
import re
import traceback
import numpy as np

from typing import Optional, List
from Bio import PDB
from Bio.PDB import Superimposer
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from fastapi import FastAPI, Response
from fastapi.middleware.cors import CORSMiddleware
from io import StringIO
from pydantic import BaseModel
from starlette.websockets import WebSocket

from backend.algorithms.adaptive_es import AdaptiveES
from backend.algorithms.params.adaptive_es_params import AdaptiveESParams
from backend.algorithms.params.es_params import ESParams
from backend.algorithms.self_adaptive_es import SelfAdaptiveES
from backend.structure.protein import Protein


app = FastAPI(root_path="/api")

origins = [
    "http://localhost:4200",
    "http://localhost",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def get_health():
    return "Server is up and running!"

@app.websocket("/simulate")
async def simulate(websocket: WebSocket):
    await websocket.accept()
    try:
        message = await websocket.receive_json()
        sequence = message["sequence"]
        method = message["method"]
        params = json.loads(re.sub(r'(?<!^)(?=[A-Z])', '_', json.dumps(message["params"])).lower())

        esp = AdaptiveESParams(**params) if method == "adaptive" else ESParams(**params)
        es = AdaptiveES(esp) if method == "adaptive" else SelfAdaptiveES(esp)

        async def send_data(generation: int, protein: Protein, sigma: float, is_last: bool) -> None:
            data = {
                "generation": generation,
                "fitness": round(protein.fitness, 2),
                "atomPositions": protein.atom_positions,
                "sequence": sequence,
                "sigma": round(sigma, 2),
                "angles": [(angle[0], angle[1]) for angle in protein.angles],
                "isLast": is_last,
            }
            await websocket.send_json(data)

        def sync_callback(generation: int, protein: Protein, sigma: float, is_last: bool) -> None:
            future = asyncio.run_coroutine_threadsafe(send_data(generation, protein, sigma, is_last), loop)
            future.result()

        loop = asyncio.get_running_loop()
        await loop.run_in_executor(None, lambda: es.run(sequence, sync_callback))

    except Exception as e:
        print(f"Error: {e}")
        traceback.print_exc()

class AnalyzeModel(BaseModel):
    cifFile1: str
    cifFile2: Optional[str] = None

@app.post("/analyze")
def perform_analysis(am: AnalyzeModel):
    cif_io = StringIO(am.cifFile1)
    parser = PDB.MMCIFParser(QUIET=True)
    chain = parser.get_structure("s", cif_io)[0].child_list[0]

    phi_psi_angles = None
    sequence = None

    polypeptides = PDB.PPBuilder().build_peptides(chain)
    for poly in polypeptides:
        phi_psi_angles = poly.get_phi_psi_list()
        sequence = str(poly.get_sequence())

    if not phi_psi_angles:
        raise Exception("No chain found")

    p = Protein(sequence, angles=[(round(np.rad2deg(phi), 2) if phi else None, round(np.rad2deg(psi), 2) if psi else None, 180) for phi, psi in phi_psi_angles])

    angles = [(float(phi) if isinstance(phi, np.float64) else phi,
               float(psi) if isinstance(psi, np.float64) else psi,
               float(omega) if isinstance(omega, np.float64) else omega) for phi, psi, omega in p.angles]

    data = {
        'angles': str(angles),
        'sequence': p.sequence,
        'cifFile': p.cif_str,
    }
    return Response(content=json.dumps(data), media_type="application/json")

@app.post("/align")
def align(am: AnalyzeModel):
    structure1 = get_structure_from_cif(am.cifFile1, "s1")
    structure2 = get_structure_from_cif(am.cifFile2, "s2")

    model1 = structure1[0]
    model2 = structure2[0]

    atoms1 = [atom for residue in model1.get_residues() for atom in residue.get_atoms()]
    atoms2 = [atom for residue in model2.get_residues() for atom in residue.get_atoms()]

    if len(atoms1) != len(atoms2):
        raise Exception("Different number of atoms")

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    sup.apply(structure2.get_atoms())

    data = {
        'sequence1': str(PDB.PPBuilder().build_peptides(structure1)[0].get_sequence()),
        'sequence2': str(PDB.PPBuilder().build_peptides(structure2)[0].get_sequence()),
        'angles1': '',
        'angles2': '',
        'cifFile': get_cif_from_structure(merge_structures(structure1, structure2)),
        'rmsd': sup.rms
    }
    return Response(content=json.dumps(data), media_type="application/json")

def get_cif_from_structure(structure: Structure) -> str:
    cif_io = StringIO()
    io_writer = PDB.MMCIFIO()
    io_writer.set_structure(structure)
    io_writer.save(cif_io)
    return cif_io.getvalue()

def get_structure_from_cif(cif_file: str, structure_id: str) -> Structure:
    parser = PDB.MMCIFParser(QUIET=True)
    cif_io = StringIO(cif_file)
    return parser.get_structure(structure_id, cif_io)

def merge_structures(structure1: Structure, structure2: Structure) -> Structure:
    merged_structure = Structure("merged_structure")
    merged_model = Model(0)
    add_chains(merged_model, [structure1, structure2], ['A', 'B'])
    merged_structure.add(merged_model)
    return merged_structure

def add_chains(model, structures: List[Structure], ids: List[str]) -> None:
    for structure, i in zip(structures, ids):
        for chain in structure.get_chains():
            chain.id = i
            model.add(chain)