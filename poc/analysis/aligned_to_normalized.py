import os
import numpy as np

def main():
    input_path = 'aa_cifs/aligned'
    pdb_files = [f for f in os.listdir(input_path) if f.endswith('.cif')]

    for pdb_file in sorted(pdb_files):
        current_atoms = []
        with open(os.path.join(input_path, pdb_file), 'r') as f:
            for line in f.read().splitlines()[::-1]:
                if not line.startswith('HETATM'):
                    break
                data = line.strip().split()
                current_atoms.append([data[3], np.array([float(data[10]), float(data[11]), float(data[12])])])

        n = next(atom for atom in current_atoms if atom[0] == 'N')
        normalized_atoms = []

        for atom in current_atoms:
            a = atom[1] - n[1]
            normalized_atoms.append([atom[0], a])

        print(pdb_file)
        for atom in sorted(normalized_atoms, key=lambda x: x[0]):
            print(f"{atom[0]:<4} {atom[1][0]:>6.3f}, {atom[1][1]:>6.3f}, {atom[1][2]:>6.3f}")
        print()


if __name__ == '__main__':
    main()