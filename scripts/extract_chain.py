import sys

pdb_file = sys.argv[1]
chain_id = sys.argv[2]
output_file = sys.argv[3]

with open(pdb_file, 'r') as f, open(output_file, 'w') as out:
    for line in f:
        if line.startswith('ATOM') and line[21] == chain_id:
            out.write(line)
        elif line.startswith('TER') and len(line) > 21 and line[21] == chain_id:
            out.write(line)
        elif line.startswith('END'):
            out.write(line)
            break

print(f"✓ Extracted chain {chain_id} from {pdb_file} to {output_file}")
