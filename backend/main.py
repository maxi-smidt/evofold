import asyncio
import base64
import json

from fastapi import FastAPI
from starlette.websockets import WebSocket

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
        print("entered")
        print(sequence)
        generation = 0
        while generation < 20:
            structure = Protein(sequence)
            data = {
                "generation": generation,
                "fitness": str(structure.fitness()),
                "cifFile": structure.to_cif(),
            }
            generation += 1
            print(data)
            await websocket.send_json(data)
            await asyncio.sleep(20)
    except Exception as e:
        print(f"Connection closed: {e}")