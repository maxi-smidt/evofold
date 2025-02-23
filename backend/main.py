import asyncio
import json

from Bio.Pathway.Rep.MultiGraph import df_search
from fastapi import FastAPI
from starlette.websockets import WebSocket

from backend.algorithms.evolution_strategy import EvolutionStrategy
from backend.algorithms.evolution_strategy_params import EvolutionStrategyParams
from backend.structure.protein import Protein


app = FastAPI(root_path="/api")


@app.get("/health")
def get_health():
    return "Server is up and running!"

@app.websocket("/simulate")
async def simulate(websocket: WebSocket):
    await websocket.accept()
    try:
        sequence = json.loads(await websocket.receive_text())
        esp = EvolutionStrategyParams()
        es = EvolutionStrategy(esp)

        async def send_data(generation: int, protein: Protein, sigma: float) -> None:
            data = {
                "generation": generation,
                "fitness": protein.fitness,
                "cifFile": protein.cif_str,
                "sequence": sequence
            }
            await websocket.send_json(data)

        def sync_callback(generation: int, protein: Protein, sigma: float):
            future = asyncio.run_coroutine_threadsafe(send_data(generation, protein, sigma), loop)
            future.result()

        loop = asyncio.get_running_loop()
        await loop.run_in_executor(None, lambda: es.run(sequence, sync_callback))

    except Exception as e:
        print(f"Connection closed: {e}")