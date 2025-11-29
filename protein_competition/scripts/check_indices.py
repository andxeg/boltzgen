import sys
import os


def check_structure_indices(cif_path):
    print(f"Checking file: {cif_path}")
    if not os.path.exists(cif_path):
        print("Error: File not found.")
        return

    # Critical Residues for Nipah G (Author Numbering)
    # We want to find what 'label_seq_id' matches these 'auth_seq_id's
    target_residues = {
        '242': 'ARG (Rim)',
        '504': 'TRP (Trp Gate)',  # <--- MOST IMPORTANT
        '505': 'GLU (Loop)',
        '506': 'GLY (Loop)',
        '533': 'GLU (Salt Bridge)',
        '558': 'ALA (Pocket)',
        '559': 'GLN (Pocket)',
        '580': 'ILE (Floor)',
        '581': 'TYR (Floor)',
        '588': 'ILE (Floor)'
    }

    found_mapping = {}

    with open(cif_path, 'r') as f:
        lines = f.readlines()

    # 1. Parse Headers to find column positions
    # We look for the loop_ that contains _atom_site.label_seq_id
    headers = {}
    in_loop = False
    data_start = 0

    for i, line in enumerate(lines):
        line = line.strip()
        if line == "loop_":
            in_loop = True
            headers = {}  # Reset headers for new loop
            continue

        if line.startswith("_atom_site."):
            headers[line] = len(headers)
        elif in_loop and not line.startswith("_"):
            # We found the data block
            if "_atom_site.id" in headers:  # Confirm it's the atom site loop
                data_start = i
                break

    if "_atom_site.label_seq_id" not in headers:
        print("CRITICAL ERROR: Could not find _atom_site definitions in CIF.")
        return

    # Get column indices
    idx_auth = headers.get('_atom_site.auth_seq_id')
    idx_label = headers.get('_atom_site.label_seq_id')
    idx_res = headers.get('_atom_site.label_comp_id')

    print(f"{'AUTH':<10} | {'RES':<5} | {'LABEL (Use in YAML)':<20} | {'ROLE'}")
    print("-" * 60)

    # 2. Scan Data
    for i in range(data_start, len(lines)):
        parts = lines[i].split()
        if len(parts) < len(headers): continue
        if parts[0] == '#': break  # End of loop structure

        curr_auth = parts[idx_auth]
        curr_label = parts[idx_label]
        curr_res = parts[idx_res]

        if curr_auth in target_residues and curr_auth not in found_mapping:
            role = target_residues[curr_auth]
            print(f"{curr_auth:<10} | {curr_res:<5} | {curr_label:<20} | {role}")
            found_mapping[curr_auth] = curr_label

            # Stop if we found everything
            if len(found_mapping) == len(target_residues):
                break

    # 3. Output the Config Line
    print("-" * 60)
    if found_mapping:
        sorted_labels = [found_mapping[k] for k in target_residues.keys() if k in found_mapping]
        print("\nSUCCESS! Paste this into your YAML constraints:")
        print(f"binding: {','.join(sorted_labels)}")
    else:
        print("FAILED: No matching residues found. Check if Chain A exists in the file.")


if __name__ == "__main__":
    check_structure_indices(
        "/lambda/nfs/Nipah-hackathon/andrei/boltzgen/protein_competition/pdb/2VSM_target_only_A.cif"
    )
