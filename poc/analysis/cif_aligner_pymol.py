import os
from pymol import cmd


def align_and_save(input_dir):
    pdb_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.cif')])

    reference_structure = pdb_files[0]
    reference_path = os.path.join(input_dir, reference_structure)
    ref_aa = os.path.splitext(reference_structure)[0]

    print(f"Loading reference structure: {reference_structure}")
    cmd.load(reference_path)

    for pdb_file in pdb_files[1:]:
        aa = os.path.splitext(pdb_file)[0]

        print(f"Loading and aligning: {pdb_file}")
        target_path = os.path.join(input_dir, pdb_file)
        cmd.load(target_path)

        if not aa == 'GLY':
            cmd.pair_fit(f'{aa}///0/N',  f'{ref_aa}///0/N',
                         f'{aa}///0/CA', f'{ref_aa}///0/CA',
                         f'{aa}///0/C',  f'{ref_aa}///0/C',
                         f'{aa}///0/CB', f'{ref_aa}///0/CB')
        else:
            cmd.pair_fit(f'{aa}///0/N',  f'{ref_aa}///0/N',
                         f'{aa}///0/CA', f'{ref_aa}///0/CA',
                         f'{aa}///0/C',  f'{ref_aa}///0/C')

        aligned_filename = os.path.join(input_dir, f"aligned/aligned_{pdb_file}")
        cmd.save(aligned_filename)
        print(f"Saved aligned structure to: {aligned_filename}")

        cmd.delete(aa)
    aligned_filename = os.path.join(input_dir, f"aligned/aligned_{reference_structure}")
    cmd.save(aligned_filename)
    cmd.delete(ref_aa)
    print("All files have been aligned and saved.")


align_and_save('/Users/maximiliansmidt/PycharmProjects/evofold/poc/analysis/aa_cifs')
