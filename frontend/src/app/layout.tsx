import type {Metadata} from "next";
import React from "react";
import Link from "next/link";
import {PrimeReactProvider} from "primereact/api";
import "./globals.css";
import "primereact/resources/themes/lara-light-cyan/theme.css";
import {Button} from "primereact/button";


export const metadata: Metadata = {
  title: "EvoFold",
  description: "Compute protein structures using evolutionary algorithms",
};

export default function RootLayout({children}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
    <body className="overflow-y-hidden">
    <PrimeReactProvider>
      <div className="flex flex-column bg-gray-100">
        <div className="flex p-3 bg-gray-50 justify-content-between">
          <Link href="/" className="no-underline">
            <div className="flex gap-3 align-items-center select-none">
              <div className="text-5xl text-gray-400 font-light">EvoFold</div>
              <div className="text-xl text-gray-300 border-solid border-2 border-round-xl p-1">BETA</div>
            </div>
          </Link>

          <div className="flex align-items-center mr-3">
            <Button severity="secondary" outlined icon="pi pi-github" />
          </div>
        </div>

        <div style={{height: `calc(100vh - 75px)`}}>
          {children}
        </div>
      </div>
    </PrimeReactProvider>
    </body>
    </html>
  );
}
