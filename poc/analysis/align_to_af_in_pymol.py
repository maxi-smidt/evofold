import csv
import tempfile
from typing import List

from pymol import cmd


def load_cif_string(cif_str, object_name):
    with tempfile.NamedTemporaryFile(delete=False, suffix=".cif") as temp_file:
        temp_file.write(cif_str.encode())
        temp_file.close()
        cmd.load(temp_file.name, object_name)

def align_with_af(files: List[str], af_file: str):
    cmd.load(f"{af_file}.cif", "af_structure")
    for fn in files:
        with open(f'best_qual/{fn}.csv', mode='r') as file:
            reader = csv.reader(file)
            next(reader, None) # skip sequence
            next(reader, None) # skip header
            print(fn, end='')
            for i, row in enumerate(reader):
                cif = row[3]
                load_cif_string(cif, f"{fn}_{i}")
                cmd.remove("hydro")
                rmsd = cmd.align(f"{fn}_{i}", "af_structure")[0]
                print(rmsd, end=', ')
            print("\n=======")

# ['comma-amber', 'plus-amber', 'comma-charmm']
align_with_af(['plus-charmm'], 'best_qual/fold_peptide_v2_model_0')