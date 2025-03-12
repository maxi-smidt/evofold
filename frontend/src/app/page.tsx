'use client';

import {InputTextarea} from "primereact/inputtextarea";
import {Button} from "primereact/button";
import {Message} from "primereact/message";
import {ChangeEvent, useState} from "react";
import {useRouter} from "next/navigation";

export default function Home() {
  const aminoAcids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];
  const errorMessage: string = "The sequence can only contain 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', " +
    "'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'.";
  const [sequence, setSequence] = useState<string>('');
  const [hasError, setHasError] = useState<boolean>();
  const router = useRouter();

  const onInput = (e: ChangeEvent<HTMLTextAreaElement>) => {
    const inputElement = e.target as HTMLTextAreaElement;
    const modifiedSequence = inputElement.value.toUpperCase();
    inputElement.value = modifiedSequence;
    setSequence(modifiedSequence);

    setHasError([...modifiedSequence].some(char => !aminoAcids.includes(char)));
  }

  return (
    <div className="flex flex-column h-full justify-content-center align-items-center gap-2">
      <div className="flex gap-3 align-items-center">
        <InputTextarea
          placeholder="Enter your sequence..."
          value={sequence}
          onChange={(e) => onInput(e)}
          rows={5}
          cols={92}/>

        <Button label="Simulate" severity="secondary" outlined disabled={hasError || sequence.trim().length === 0}
                onClick={() => router.push(`/structure/${sequence}`)}/>
      </div>
      { hasError &&
        <Message severity="error" text={errorMessage} />
      }
    </div>
  );
}
