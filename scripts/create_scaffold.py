import sys


def create_scaffold(input_pdb, output_pdb, chain_mappings):
    """Extracts and renames chains from input_pdb to output_pdb."""

    chains_to_extract = set(chain_mappings.keys())

    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM') or line.startswith('TER'):
                # Extract original chain ID from column 22 (index 21)
                if len(line) > 21:
                    original_chain_id = line[21]

                    if original_chain_id in chains_to_extract:
                        new_chain_id = chain_mappings[original_chain_id]
                        # Replace original chain ID with new ID
                        new_line = line[:21] + new_chain_id + line[22:]
                        f_out.write(new_line)

            elif line.startswith('END'):
                f_out.write(line)
                break

    print(f"✓ Created {output_pdb}. Chains {list(chains_to_extract)} renamed to {list(chain_mappings.values())}.")


if __name__ == "__main__":
    # Usage: python create_scaffold.py <input.pdb> <output.pdb> <map1> <map2> ...
    # Example: python create_scaffold.py 8u4q_full.pdb 8u4q_fab.pdb B:H A:L

    if len(sys.argv) < 4:
        print("Usage: python create_scaffold.py <input.pdb> <output.pdb> <old:new_map> [old:new_map...]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    mapping_args = sys.argv[3:]

    mappings = {}
    try:
        for arg in mapping_args:
            original, new = arg.split(':')
            mappings[original] = new
    except:
        print("Error: Mappings must be in 'ORIG:NEW' format (e.g., 'B:H')")
        sys.exit(1)

    create_scaffold(input_file, output_file, mappings)
