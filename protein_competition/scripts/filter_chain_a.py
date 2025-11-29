#!/usr/bin/env python3
"""
Comprehensive mmCIF filter to keep only chain A.
Removes all other chains, water (HOH), ligands, and associated metadata.

This script handles:
1. _atom_site section - keeps only chain A atoms
2. _struct_conn section - removes connections involving other chains
3. _struct_asym section - removes other chain definitions
4. All other sections - kept as-is

Usage:
    python filter_chain_a.py
"""

import sys
from pathlib import Path


def filter_mmcif_chain_a(input_path, output_path, keep_chain='A'):
    """
    Filter mmCIF file to keep only specified chain.

    Args:
        input_path: Path to input mmCIF file
        output_path: Path to output filtered mmCIF file
        keep_chain: Chain ID to keep (default: 'A')
    """
    print(f"Reading: {input_path}")

    with open(input_path, 'r') as f:
        lines = f.readlines()

    output_lines = []

    # State tracking
    in_atom_site_data = False
    in_struct_conn_data = False
    in_struct_asym_data = False
    atom_site_header_line = None
    struct_conn_header_line = None
    struct_asym_header_line = None

    for i, line in enumerate(lines):
        # ===== Track section headers =====

        # _atom_site section (atom coordinates)
        if line.startswith('_atom_site.'):
            output_lines.append(line)
            if atom_site_header_line is None:
                atom_site_header_line = i
            continue

        # _struct_conn section (disulfide bonds, covalent bonds)
        if line.startswith('_struct_conn.'):
            output_lines.append(line)
            if struct_conn_header_line is None:
                struct_conn_header_line = i
            continue

        # _struct_asym section (chain/asymmetric unit definitions)
        if line.startswith('_struct_asym.'):
            output_lines.append(line)
            if struct_asym_header_line is None:
                struct_asym_header_line = i
            continue

        # ===== Detect when we enter data sections =====

        if atom_site_header_line is not None and not in_atom_site_data:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                in_atom_site_data = True

        if struct_conn_header_line is not None and not in_struct_conn_data:
            # Connection data starts with disulf, covale, metalc, etc.
            if (line.startswith('disulf') or line.startswith('covale') or
                line.startswith('metalc') or line.startswith('hydrog')):
                in_struct_conn_data = True

        if struct_asym_header_line is not None and not in_struct_asym_data:
            # struct_asym data: non-header, non-comment lines after headers
            if not line.startswith('_') and not line.startswith('#') and line.strip():
                parts = line.split()
                # Check if this looks like chain ID (A, B, C, etc.)
                if len(parts) >= 1 and len(parts[0]) == 1 and parts[0].isalpha():
                    in_struct_asym_data = True

        # ===== Detect when we leave data sections =====

        if in_atom_site_data:
            if line.startswith('#') or (line.startswith('_') and not line.startswith('_atom_site')):
                in_atom_site_data = False
                atom_site_header_line = None
                output_lines.append(line)
                continue

        if in_struct_conn_data:
            if line.startswith('#') or line.startswith('_'):
                in_struct_conn_data = False
                struct_conn_header_line = None
                output_lines.append(line)
                continue

        if in_struct_asym_data:
            if line.startswith('#') or line.startswith('_'):
                in_struct_asym_data = False
                struct_asym_header_line = None
                output_lines.append(line)
                continue

        # ===== Process data in each section =====

        if in_atom_site_data:
            # Filter atom records - keep only specified chain
            if not line.strip():
                output_lines.append(line)
                continue

            parts = line.split()
            if len(parts) < 7:
                output_lines.append(line)
                continue

            # Field 7 (index 6) is label_asym_id (chain identifier)
            label_asym_id = parts[6]

            if label_asym_id == keep_chain:
                output_lines.append(line)
            # else: skip this atom (different chain)

        elif in_struct_conn_data:
            # Filter structural connections - keep only if both partners are in keep_chain
            if not line.strip():
                output_lines.append(line)
                continue

            parts = line.split()
            if len(parts) >= 14:
                # Field 5 (index 4): ptnr1_label_asym_id
                # Field 13 (index 12): ptnr2_label_asym_id
                ptnr1_asym = parts[4]
                ptnr2_asym = parts[12]

                # Keep only if both partners are in the keep_chain
                if ptnr1_asym == keep_chain and ptnr2_asym == keep_chain:
                    output_lines.append(line)
                # else: skip this connection (involves removed chains)
            else:
                # Malformed line, keep it to avoid breaking file structure
                output_lines.append(line)

        elif in_struct_asym_data:
            # Filter chain definitions - keep only specified chain
            if not line.strip():
                output_lines.append(line)
                continue

            parts = line.split()
            if len(parts) >= 1:
                asym_id = parts[0]

                if asym_id == keep_chain:
                    output_lines.append(line)
                # else: skip this chain definition
            else:
                output_lines.append(line)

        else:
            # Not in any filtered section - keep all lines
            output_lines.append(line)

    # Write filtered output
    print(f"Writing: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    print(f"✓ Successfully filtered mmCIF file")
    print(f"✓ Kept only chain {keep_chain}")
    print(f"✓ Removed all other chains, water (HOH), and ligands")

    return output_path


def main():
    """Main entry point."""
    input_file = Path('/lambda/nfs/Nipah-hackathon/andrei/boltzgen/protein_competition/pdb/2VSM.cif')
    output_file = Path('/lambda/nfs/Nipah-hackathon/andrei/boltzgen/protein_competition/pdb/2VSM_target_only_A.cif')

    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    filter_mmcif_chain_a(input_file, output_file, keep_chain='A')

    # Verify output
    print("\n=== Verification ===")
    import subprocess

    # Count atoms in original
    result_orig = subprocess.run(
        ['grep', '-c', '^ATOM\\|^HETATM', str(input_file)],
        capture_output=True, text=True
    )
    orig_atoms = result_orig.stdout.strip()

    # Count atoms in filtered
    result_filt = subprocess.run(
        ['grep', '-c', '^ATOM\\|^HETATM', str(output_file)],
        capture_output=True, text=True
    )
    filt_atoms = result_filt.stdout.strip()

    print(f"Original atoms: {orig_atoms}")
    print(f"Filtered atoms: {filt_atoms}")

    # Check which chains remain
    result_chains = subprocess.run(
        f"grep '^ATOM\\|^HETATM' {output_file} | awk '{{print $7}}' | sort -u",
        shell=True, capture_output=True, text=True
    )
    chains = result_chains.stdout.strip()
    print(f"Chains in filtered file: {chains}")


if __name__ == '__main__':
    main()
