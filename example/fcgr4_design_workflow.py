#!/usr/bin/env python3
"""
FCGR4 Design Workflow Script
Performs structure-based design of binders against FCGR4 using 3AY4 binders as scaffolds.
"""

import sys
from pathlib import Path
from typing import List, Dict, Set
from collections import defaultdict
import logging

try:
    import yaml
    import gemmi
    from Bio import PDB
    from Bio.PDB import Superimposer, PDBIO, Select
    from Bio.SeqUtils import seq1
    import numpy as np
except ImportError as e:
    print(f"Error: Missing required package. {e}")
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class ChainSelector(Select):
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids


class StructureProcessor:
    def __init__(self, interface_cutoff: float = 4.5):
        self.interface_cutoff = interface_cutoff
        self.parser = PDB.MMCIFParser(QUIET=True)
        self.pdb_parser = PDB.PDBParser(QUIET=True)
        self.io = PDBIO()
        self.hydrophobic = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}
    
    def _get_parser(self, path: Path):
        return self.parser if path.suffix == '.cif' else self.pdb_parser
    
    def extract_chains_from_cif(self, cif_path: Path, chain_ids: List[str], output_path: Path):
        """Extract chains from CIF and save as PDB"""
        logger.info(f"Extracting chains {chain_ids} from {cif_path}")
        doc = gemmi.read_structure(str(cif_path))
        new_doc = gemmi.Structure()
        new_doc.name = doc.name
        
        for model in doc:
            new_model = gemmi.Model(model.name)
            for chain in model:
                if chain.name in chain_ids:
                    new_chain = new_model.add_chain(chain.name)
                    for residue in chain:
                        new_chain.add_residue(residue.clone())
            if len(new_model) > 0:
                new_doc.add_model(new_model)
        
        new_doc.write_pdb(str(output_path))
        logger.info(f"Saved to {output_path}")
        return output_path
    
    def calculate_interface_residues(self, structure_path: Path, chain1: str, chain2: str) -> Dict:
        """Find interface residues between two chains"""
        logger.info(f"Analyzing interface between chain {chain1} and chain {chain2}")
        structure = self._get_parser(structure_path).get_structure("complex", str(structure_path))
        model = structure[0]
        chain1_obj, chain2_obj = model[chain1], model[chain2]
        
        interface = {chain1: set(), chain2: set()}
        for res1 in chain1_obj:
            if not PDB.is_aa(res1, standard=True):
                continue
            atoms1 = [a for a in res1 if a.element != 'H']
            for res2 in chain2_obj:
                if not PDB.is_aa(res2, standard=True):
                    continue
                atoms2 = [a for a in res2 if a.element != 'H']
                min_dist = min(atom1 - atom2 for atom1 in atoms1 for atom2 in atoms2)
                if min_dist < self.interface_cutoff:
                    interface[chain1].add(res1.id[1])
                    interface[chain2].add(res2.id[1])
        
        interface[chain1] = sorted(interface[chain1])
        interface[chain2] = sorted(interface[chain2])
        logger.info(f"Found {len(interface[chain1])} and {len(interface[chain2])} interface residues")
        return {'interface_residues': interface}
    
    def align_structures(self, ref_path: Path, mobile_path: Path, ref_chain: str, mobile_chain: str):
        """Align structures and return RMSD"""
        logger.info(f"Aligning {mobile_path} to {ref_path}")
        ref_struct = self._get_parser(ref_path).get_structure("ref", str(ref_path))
        mobile_struct = self._get_parser(mobile_path).get_structure("mobile", str(mobile_path))
        
        ref_res = {r.id[1]: r for r in ref_struct[0][ref_chain] if PDB.is_aa(r, standard=True)}
        mobile_res = {r.id[1]: r for r in mobile_struct[0][mobile_chain] if PDB.is_aa(r, standard=True)}
        common = sorted(set(ref_res.keys()) & set(mobile_res.keys()))
        
        ref_atoms = [ref_res[r]['CA'] for r in common if 'CA' in ref_res[r] and 'CA' in mobile_res[r]]
        mobile_atoms = [mobile_res[r]['CA'] for r in common if 'CA' in ref_res[r] and 'CA' in mobile_res[r]]
        
        if len(ref_atoms) < 3:
            raise ValueError("Not enough common atoms for alignment")
        
        superimposer = Superimposer()
        superimposer.set_atoms(ref_atoms, mobile_atoms)
        superimposer.apply(mobile_struct[0].get_atoms())
        logger.info(f"Alignment RMSD: {superimposer.rms:.3f} Å")
        return superimposer.rotran, superimposer.rms
    
    def compute_interface_rmsd(self, struct1_path: Path, struct2_path: Path, 
                              chain1: str, chain2: str, interface_residues: Set[int]) -> Dict:
        """Compute RMSD for interface residues"""
        logger.info(f"Computing interface RMSD for {len(interface_residues)} residues")
        struct1 = self._get_parser(struct1_path).get_structure("s1", str(struct1_path))
        struct2 = self._get_parser(struct2_path).get_structure("s2", str(struct2_path))
        
        chain1_obj, chain2_obj = struct1[0][chain1], struct2[0][chain2]
        atoms1, atoms2 = [], []
        residue_rmsd = {}
        
        for res_num in interface_residues:
            if res_num in chain1_obj and res_num in chain2_obj:
                r1, r2 = chain1_obj[res_num], chain2_obj[res_num]
                if PDB.is_aa(r1, standard=True) and PDB.is_aa(r2, standard=True):
                    if 'CA' in r1 and 'CA' in r2:
                        atoms1.append(r1['CA'])
                        atoms2.append(r2['CA'])
        
        if len(atoms1) < 3:
            return {'rmsd': float('inf'), 'residue_rmsd': {}}
        
        coords1 = np.array([a.coord for a in atoms1])
        coords2 = np.array([a.coord for a in atoms2])
        coords1 -= coords1.mean(axis=0)
        coords2 -= coords2.mean(axis=0)
        
        diff = coords1 - coords2
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        
        for i, res_num in enumerate(interface_residues):
            if i < len(atoms1):
                residue_rmsd[res_num] = np.sqrt(np.sum((coords1[i] - coords2[i])**2))
        
        logger.info(f"Interface RMSD: {rmsd:.3f} Å")
        return {'rmsd': rmsd, 'residue_rmsd': residue_rmsd}
    
    def identify_hydrophobic_patches(self, structure_path: Path, chain: str, 
                                     interface_residues: Set[int]) -> Set[int]:
        """Find hydrophobic residues in interface"""
        structure = self._get_parser(structure_path).get_structure("struct", str(structure_path))
        chain_obj = structure[0][chain]
        hydrophobic = {r.id[1] for r in chain_obj 
                      if r.id[1] in interface_residues and PDB.is_aa(r, standard=True)
                      and r.get_resname() in self.hydrophobic}
        logger.info(f"Found {len(hydrophobic)} hydrophobic interface residues")
        return hydrophobic
    
    def extract_sequence(self, structure_path: Path, chain_id: str) -> str:
        """Extract amino acid sequence from structure"""
        structure = self._get_parser(structure_path).get_structure("struct", str(structure_path))
        chain = structure[0][chain_id]
        sequence = ''.join(seq1(r.get_resname()) for r in chain 
                          if PDB.is_aa(r, standard=True) and r.get_resname())
        return sequence
    
    def extract_sequence_with_mapping(self, structure_path: Path, chain_id: str) -> tuple:
        """Extract amino acid sequence and mapping from PDB residue numbers to sequence positions.
        
        Returns:
            (sequence, pdb_to_seq_map) where pdb_to_seq_map maps PDB residue number -> sequence position (1-indexed)
        """
        structure = self._get_parser(structure_path).get_structure("struct", str(structure_path))
        chain = structure[0][chain_id]
        sequence = []
        pdb_to_seq = {}
        seq_pos = 1  # 1-indexed sequence position
        
        for residue in chain:
            if PDB.is_aa(residue, standard=True) and residue.get_resname():
                pdb_res_num = residue.id[1]
                aa = seq1(residue.get_resname())
                sequence.append(aa)
                pdb_to_seq[pdb_res_num] = seq_pos
                seq_pos += 1
        
        return ''.join(sequence), pdb_to_seq


def expand_residue_set(residues: Set[int], window: int = 5) -> Set[int]:
    """Expand residue set with window"""
    expanded = set(residues)
    for res in residues:
        expanded.update(range(res - window, res + window + 1))
    return expanded


def format_residue_ranges(residues: List[int]) -> str:
    """Format sorted residue list as ranges (e.g., '1..5,10,15..20')"""
    if not residues:
        return ""
    ranges = []
    start = end = residues[0]
    for res in residues[1:]:
        if res == end + 1:
            end = res
        else:
            ranges.append(f"{start}..{end}" if start != end else str(start))
            start = end = res
    ranges.append(f"{start}..{end}" if start != end else str(start))
    return ",".join(ranges)


def interleave_design_ranges_in_sequence(sequence: str, design_residues: Set[int]) -> str:
    """Interleave design ranges with fixed amino acids in sequence.
    
    Args:
        sequence: Full amino acid sequence (1-indexed positions)
        design_residues: Set of residue positions (1-indexed) that should be designed
    
    Returns:
        String with design ranges interspersed with fixed amino acids
        Example: "AB3..5FG8..10KLMNOP" where 3..5 and 8..10 are design regions
    """
    if not design_residues:
        return sequence
    
    # Convert to 0-indexed for sequence access
    design_set = {pos - 1 for pos in design_residues if 1 <= pos <= len(sequence)}
    
    if not design_set:
        return sequence
    
    result_parts = []
    i = 0
    seq_len = len(sequence)
    
    while i < seq_len:
        if i in design_set:
            # Start of a design region - find the end
            start = i
            while i < seq_len and i in design_set:
                i += 1
            end = i - 1
            
            # Format the design range (convert back to 1-indexed)
            if start == end:
                result_parts.append(str(start + 1))
            else:
                result_parts.append(f"{start + 1}..{end + 1}")
        else:
            # Fixed region - collect consecutive fixed amino acids
            fixed_start = i
            while i < seq_len and i not in design_set:
                i += 1
            fixed_end = i
            result_parts.append(sequence[fixed_start:fixed_end])
    
    return "".join(result_parts)


def create_design_spec(fcgr4_binding_site: Set[int], binder_interface: Dict[str, Set[int]],
                      binder_sequences: Dict[str, str]) -> Dict:
    """Create design specification YAML"""
    expanded_fcgr4 = sorted(expand_residue_set(fcgr4_binding_site, window=7))
    expanded_binder = {chain: sorted(expand_residue_set(residues, window=7))
                       for chain, residues in binder_interface.items()}
    
    entities = [{
        'file': {
            'path': 'fc4_model_0.pdb',
            'include': [{'chain': {'id': 'A'}}],
            'binding_types': [{'chain': {'id': 'A', 'binding': format_residue_ranges(expanded_fcgr4)}}],
            'structure_groups': [{'group': {'visibility': 1, 'id': 'A'}}]
        }
    }]
    
    # Create a single protein entity with ID 'G' combining chains A and B
    # Use chain A as the primary scaffold (or combine both if needed)
    if 'A' in binder_sequences:
        # Get design residues for chain A (convert list to set for interleaving)
        design_residues_list = expanded_binder.get('A', [])
        design_residues = set(design_residues_list)
        sequence = binder_sequences['A']
        
        # Interleave design ranges with fixed amino acids in sequence
        sequence_with_design = interleave_design_ranges_in_sequence(sequence, design_residues)
        
        # Create protein entity without binding_types
        protein = {'protein': {'id': 'G', 'sequence': sequence_with_design}}
        entities.append(protein)
    
    return {'entities': entities, 'constraints': []}


def main():
    base_dir = Path(__file__).parent
    cif_path = base_dir / "3AY4.cif"
    fcgr4_path = base_dir / "fc4_model_0.pdb"
    work_dir = base_dir / "fcgr4_workflow_output"
    work_dir.mkdir(exist_ok=True)
    
    logger.info("=" * 80)
    logger.info("FCGR4 Design Workflow")
    logger.info("=" * 80)
    
    processor = StructureProcessor(interface_cutoff=4.5)
    
    # Step 1: Extract chains
    logger.info("\n=== Step 1: Separating chains ===")
    fcgr3_path = work_dir / "fcgr3_chain_C.pdb"
    processor.extract_chains_from_cif(cif_path, ['C'], fcgr3_path)
    processor.extract_chains_from_cif(cif_path, ['A', 'B'], work_dir / "binder_chains_AB.pdb")
    
    # Step 2: Align structures
    logger.info("\n=== Step 2: Aligning structures ===")
    fcgr4_struct = processor.pdb_parser.get_structure("fcgr4", str(fcgr4_path))
    fcgr3_struct = processor.pdb_parser.get_structure("fcgr3", str(fcgr3_path))
    
    ref_res = {r.id[1]: r for r in fcgr4_struct[0]['A'] if PDB.is_aa(r, standard=True)}
    mobile_res = {r.id[1]: r for r in fcgr3_struct[0]['C'] if PDB.is_aa(r, standard=True)}
    common = sorted(set(ref_res.keys()) & set(mobile_res.keys()))
    
    ref_atoms = [ref_res[r]['CA'] for r in common if 'CA' in ref_res[r] and 'CA' in mobile_res[r]]
    mobile_atoms = [mobile_res[r]['CA'] for r in common if 'CA' in ref_res[r] and 'CA' in mobile_res[r]]
    
    superimposer = Superimposer()
    superimposer.set_atoms(ref_atoms, mobile_atoms)
    superimposer.apply(fcgr3_struct[0].get_atoms())
    rmsd = superimposer.rms
    logger.info(f"Overall alignment RMSD: {rmsd:.3f} Å")
    
    aligned_fcgr3_path = work_dir / "fcgr3_aligned.pdb"
    processor.io.set_structure(fcgr3_struct)
    processor.io.save(str(aligned_fcgr3_path), ChainSelector(['C']))
    
    # Step 3: Identify interface residues
    logger.info("\n=== Step 3: Identifying hotspot residues ===")
    interface_ca = processor.calculate_interface_residues(cif_path, 'C', 'A')
    interface_cb = processor.calculate_interface_residues(cif_path, 'C', 'B')
    
    fcgr3_interface = set(interface_ca['interface_residues']['C']) | set(interface_cb['interface_residues']['C'])
    binder_a = set(interface_ca['interface_residues']['A'])
    binder_b = set(interface_cb['interface_residues']['B'])
    
    logger.info(f"FCGR3 interface: {sorted(fcgr3_interface)}")
    logger.info(f"Binder A interface: {sorted(binder_a)}")
    logger.info(f"Binder B interface: {sorted(binder_b)}")
    
    # Step 4: Compute interface RMSD
    logger.info("\n=== Step 4: Computing interface RMSD ===")
    interface_rmsd = processor.compute_interface_rmsd(
        aligned_fcgr3_path, fcgr4_path, 'C', 'A', fcgr3_interface)
    
    # Step 5: Identify binding site
    logger.info("\n=== Step 5: Identifying low RMSD and hydrophobic patches ===")
    low_rmsd = {r for r, val in interface_rmsd['residue_rmsd'].items() if val < 4.0}
    hydrophobic = processor.identify_hydrophobic_patches(fcgr4_path, 'A', low_rmsd)
    
    fcgr4_binding_site = low_rmsd & hydrophobic
    if not fcgr4_binding_site:
        fcgr4_binding_site = low_rmsd if low_rmsd else fcgr3_interface
        logger.warning("Using fallback binding site residues")
    
    logger.info(f"FCGR4 binding site: {sorted(fcgr4_binding_site)}")
    
    # Step 6: Create design specification
    logger.info("\n=== Step 6: Creating design specification ===")
    logger.info("Extracting binder sequences with residue mappings...")
    binder_sequences = {}
    binder_mappings = {}
    for chain in ['A', 'B']:
        seq, mapping = processor.extract_sequence_with_mapping(cif_path, chain)
        binder_sequences[chain] = seq
        binder_mappings[chain] = mapping
        logger.info(f"Chain {chain}: {len(seq)} residues ({seq[:50]}...)")
    
    # Map PDB residue numbers to sequence positions for design residues
    binder_a_seq_positions = {binder_mappings['A'][pdb_num] for pdb_num in binder_a 
                              if pdb_num in binder_mappings['A']}
    binder_b_seq_positions = {binder_mappings['B'][pdb_num] for pdb_num in binder_b 
                              if pdb_num in binder_mappings['B']}
    
    design_spec = create_design_spec(
        fcgr4_binding_site,
        {'A': binder_a_seq_positions, 'B': binder_b_seq_positions},
        binder_sequences
    )
    
    output_yaml = work_dir / "Fcgr4_binder_design.yaml"
    with open(output_yaml, 'w') as f:
        yaml.dump(design_spec, f, default_flow_style=False, sort_keys=False)
    
    logger.info(f"Design specification saved to {output_yaml}")
    logger.info("\n" + "=" * 80)
    logger.info("Workflow Summary")
    logger.info("=" * 80)
    logger.info(f"FCGR3 interface residues: {len(fcgr3_interface)}")
    logger.info(f"FCGR4 binding site residues: {len(fcgr4_binding_site)}")
    logger.info(f"Binder A design residues: {len(binder_a)}")
    logger.info(f"Binder B design residues: {len(binder_b)}")
    logger.info(f"Overall alignment RMSD: {rmsd:.3f} Å")
    logger.info(f"Interface RMSD: {interface_rmsd['rmsd']:.3f} Å")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
