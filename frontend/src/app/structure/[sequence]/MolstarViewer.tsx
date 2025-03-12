"use client";


import {PluginContext} from "molstar/lib/mol-plugin/context";
import {DefaultPluginSpec} from "molstar/lib/mol-plugin/spec";
import {useEffect} from "react";

export default function MolstarViewer({cifFile, showCartoon}: { cifFile: string, showCartoon: boolean }) {

  useEffect(() => {
    createRootViewer().then(async (plugin) => {
      const fileData = await plugin.builders.data.rawData({ data: cifFile, label: 'Custom CIF' });
      const trajectory = await plugin.builders.structure.parseTrajectory(fileData, "mmcif");
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);

      await plugin.builders.structure.representation.addRepresentation(
        structure,
        {
          type: "cartoon",
          typeParams: {aspectRatio: 1, sizeFactor: 0.5},
          color: "sequence-id",
        },
        {tag: "cartoon-style"}
      );

      await plugin.builders.structure.representation.addRepresentation(
        structure,
        {
          type: "ball-and-stick",
          color: "element-symbol",
          colorParams: {
            carbonColor: {name: "element-symbol", params: {}}
          },
        },
        {tag: "ballstick-style"}
      );
    });
  }, [cifFile, showCartoon]);

  return (
    <div id="app" className="flex w-full flex-column">
      <canvas id="canvas" style={{width: "100%", height: "100%"}}></canvas>
    </div>
  );

  async function createRootViewer() {
    const viewport = document.getElementById("app") as HTMLDivElement;
    const canvas = document.getElementById("canvas") as HTMLCanvasElement;

    const plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();

    if (!plugin.initViewer(canvas, viewport)) {
      viewport.innerHTML = "Failed to init Mol*";
      throw new Error("init failed");
    }
    // @ts-expect-error just
    window["molstar"] = plugin;

    return plugin;
  }
}
