'use client';

import {SimulationData} from "@/app/types";
import {useEffect, useState} from "react";
import {Button} from "primereact/button";
import {useParams, useRouter} from "next/navigation";
import MolstarViewer from "@/app/structure/[sequence]/MolstarViewer";

export default function Structure() {
  const params = useParams<{ sequence: string }>();
  const router = useRouter();
  const [results, setResults] = useState<Array<SimulationData>>([]);
  const [selectedResult, setSelectedResult] = useState<SimulationData | null>(null);

  useEffect(() => {
    if (!params.sequence) router.push('/');

    const ws = new WebSocket("ws://localhost:8000/api/simulate");

    ws.onopen = () => {
      ws.send(JSON.stringify(params.sequence));
    };

    ws.onmessage = (event) => handleMessage(event);

    ws.onclose = () => {
      console.log("WebSocket disconnected");
    };

    return () => {
      ws.close();
    };
  }, [])

  const handleMessage = (message: MessageEvent) => {
    const data = JSON.parse(message.data) as SimulationData;
    let len = 0;
    setResults(prevResults => {
      len = prevResults.length;
      return [...prevResults, data]
    });

    if (!len) {
      onResultClick(data)
    }
  }

  const onResultClick = (data: SimulationData) => {
    setSelectedResult(data);
  }

  return (
    <div className="flex h-full">
      <div className="flex flex-column bg-gray-200 w-2 p-1 gap-1 overflow-y-scroll">
        {results.map((data) => {
          return (
            <div key={data.generation} className="flex border border-1 border-round-lg align-items-center hover:bg-gray-300
                  justify-content-center p-3 select-none cursor-pointer" onClick={() => onResultClick(data)}>
              <div className="text-gray-500">Generation {data.generation}</div>
            </div>)
        })}
      </div>

      <div className="flex-column-reverse w-full flex-column">
        <div className="flex h-11rem bg-gray-600">
          {
            selectedResult &&
              <div className="p-2 flex items-start justify-between w-full">
                  <table className="border-separate border-spacing-y-2">
                      <tbody>
                      <tr>
                          <td className="text-lg text-gray-200">Generation:</td>
                          <td className="text-lg text-gray-200">{selectedResult.generation}</td>
                      </tr>
                      <tr>
                          <td className="text-lg text-gray-200">Fitness:</td>
                          <td className="text-lg text-gray-200">{selectedResult.fitness}</td>
                      </tr>
                      <tr>
                          <td className="text-lg text-gray-200">Sigma:</td>
                          <td className="text-lg text-gray-200">{selectedResult.sigma}</td>
                      </tr>
                      <tr>
                          <td className="text-lg text-gray-200">Sequence:</td>
                          <td className="text-lg text-gray-200">{selectedResult.sequence}</td>
                      </tr>
                      </tbody>
                  </table>
                  <div className="ml-auto">
                    <Button outlined severity="secondary" icon="pi pi-download" onClick={() => onDownloadClick(selectedResult!.cifFile)}/>
                  </div>
              </div>
          }
        </div>
        <div className="h-full">
          { selectedResult && <MolstarViewer cifFile={selectedResult.cifFile} showCartoon={false} />
          }
        </div>
      </div>
    </div>
  );

  function onDownloadClick(cif: string) {
    const newBlob = new Blob([cif], {type: "text/plain"});
    const data = window.URL.createObjectURL(newBlob);
    const link = document.createElement("a");
    link.href = data;
    link.download = `EvoFold_${new Date().toISOString()}.cif`;
    link.click();
  }
}
