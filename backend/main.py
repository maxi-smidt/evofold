import asyncio
import json
import re
import traceback
from io import StringIO

import numpy as np
from Bio import PDB
from fastapi import FastAPI, Response
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from starlette.websockets import WebSocket

from backend.algorithms.evolution_strategy import EvolutionStrategy
from backend.algorithms.evolution_strategy_params import EvolutionStrategyParams
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
        params = json.loads(re.sub(r'(?<!^)(?=[A-Z])', '_', json.dumps(message["params"])).lower())

        esp = EvolutionStrategyParams(**params)
        es = EvolutionStrategy(esp)

        async def send_data(generation: int, protein: Protein, sigma: float) -> None:
            data = {
                "generation": generation,
                "fitness": round(protein.fitness, 2),
                "cifFile": protein.cif_str,
                "sequence": sequence,
                "sigma": round(sigma, 2),
            }
            await websocket.send_json(data)

        def sync_callback(generation: int, protein: Protein, sigma: float):
            future = asyncio.run_coroutine_threadsafe(send_data(generation, protein, sigma), loop)
            future.result()

        loop = asyncio.get_running_loop()
        await loop.run_in_executor(None, lambda: es.run(sequence, sync_callback))

    except Exception as e:
        print(f"Error: {e}")
        traceback.print_exc()

class AnalyzeModel(BaseModel):
    cifFile: str

@app.post("/analyze")
def perform_analysis(am: AnalyzeModel):
    cif_io = StringIO(am.cifFile)
    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("s", cif_io)
    model = structure[0]

    phi_psi_angles = None
    sequence = None
    for chain in model:
        polypeptides = PDB.PPBuilder().build_peptides(chain)
        for poly in polypeptides:
            phi_psi_angles = poly.get_phi_psi_list()
            sequence = str(poly.get_sequence())

    if not phi_psi_angles:
        raise Exception("No chain found")

    angles = []
    for phi, psi in phi_psi_angles:
        angles.append((round(np.rad2deg(phi), 2) if phi else None, round(np.rad2deg(psi), 2) if psi else None, 180))
    p = Protein(sequence, angles)

    data = {
        'angles': p.angles,
        'sequence': p.sequence,
        'cifFile': p.cif_str,
    }
    return Response(content=json.dumps(data), media_type="application/json")