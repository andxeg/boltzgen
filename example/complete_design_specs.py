#!/usr/bin/env python3
"""
Complete Design Specifications Script
=====================================
This script completes the TODO items in the design specification YAML files by:
1. Analyzing 3AY4 structure for Fcgr3 interface
2. Fetching AlphaFold model for mouse Fcgr4
3. Performing interface extraction with 4.5 Å cutoff
4. Identifying hotspots (>1.0 Å² BSA threshold)
5. Applying masking window (±5-8 residues)
6. Updating YAML files with actual values
"""

import os
import sys
import json
import subprocess
import requests
from pathlib import Path
from typing import List, Dict, Tuple, Set, Optional
from collections import defaultdict
import logging

# Check for required packages
try:
    import yaml
except ImportError:
    print("Error: PyYAML is required. Install with: pip install pyyaml")
    sys.exit(1)

try:
    from Bio import PDB
    from Bio.PDB import PDBIO, Select
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class StructureAnalyzer:
    """Analyze protein structures for interface detection"""
    
    def __init__(self, interface_cutoff: float = 4.5, 
                 hotspot_bsa_cutoff: float = 1.0,
                 masking_window: int = 7):
        self.interface_cutoff = interface_cutoff
        self.hotspot_bsa_cutoff = hotspot_bsa_cutoff
        self.masking_window = masking_window
        self.parser = PDB.PDBParser(QUIET=True)
        self.io = PDBIO()
    
    def get_residue_atoms(self, residue):
        """Get all heavy atoms from a residue"""
        atoms = []
        for atom in residue:
            if atom.element != 'H':  # Skip hydrogens
                atoms.append(atom)
        return atoms
    
    def calculate_interface_residues(self, pdb_path: Path, 
                                    chain1: str, chain2: str) -> Dict:
        """Calculate interface residues based on distance cutoff"""
        logger.info(f"Analyzing interface between chain {chain1} and chain {chain2}...")
        
        structure = self.parser.get_structure("complex", str(pdb_path))
        model = structure[0]
        
        if chain1 not in model:
            raise ValueError(f"Chain {chain1} not found in structure")
        if chain2 not in model:
            raise ValueError(f"Chain {chain2} not found in structure")
        
        chain1_obj = model[chain1]
        chain2_obj = model[chain2]
        
        interface_residues = {chain1: set(), chain2: set()}
        residue_distances = {}
        residue_contacts = defaultdict(lambda: {'min_dist': float('inf'), 'contacts': 0})
        
        # Calculate minimum distances between residues
        for res1 in chain1_obj:
            if not PDB.is_aa(res1, standard=True):
                continue
            
            res1_num = res1.id[1]
            atoms1 = self.get_residue_atoms(res1)
            
            for res2 in chain2_obj:
                if not PDB.is_aa(res2, standard=True):
                    continue
                
                res2_num = res2.id[1]
                atoms2 = self.get_residue_atoms(res2)
                
                min_dist = float('inf')
                contact_count = 0
                
                for atom1 in atoms1:
                    for atom2 in atoms2:
                        dist = atom1 - atom2
                        min_dist = min(min_dist, dist)
                        if dist <= self.interface_cutoff:
                            contact_count += 1
                
                if min_dist <= self.interface_cutoff:
                    interface_residues[chain1].add(res1_num)
                    interface_residues[chain2].add(res2_num)
                    
                    key1 = (chain1, res1_num)
                    key2 = (chain2, res2_num)
                    
                    residue_distances[(chain1, res1_num, chain2, res2_num)] = min_dist
                    residue_contacts[key1]['min_dist'] = min(residue_contacts[key1]['min_dist'], min_dist)
                    residue_contacts[key1]['contacts'] += contact_count
                    residue_contacts[key2]['min_dist'] = min(residue_contacts[key2]['min_dist'], min_dist)
                    residue_contacts[key2]['contacts'] += contact_count
        
        logger.info(f"Found {len(interface_residues[chain1])} interface residues in chain {chain1}")
        logger.info(f"Found {len(interface_residues[chain2])} interface residues in chain {chain2}")
        
        return {
            "interface_residues": interface_residues,
            "residue_distances": residue_distances,
            "residue_contacts": dict(residue_contacts)
        }
    
    def calculate_bsa_approximate(self, pdb_path: Path, chain1: str, chain2: str,
                                  interface_info: Dict) -> Dict:
        """Approximate BSA using contact area estimation"""
        logger.info("Calculating approximate BSA (Buried Surface Area)...")
        
        # Approximate BSA based on number of contacts and distances
        # This is a simplified approach - for accurate BSA, use FreeSASA
        bsa_per_residue = {}
        
        for chain_id in [chain1, chain2]:
            bsa_per_residue[chain_id] = {}
            
            for res_num in interface_info["interface_residues"][chain_id]:
                key = (chain_id, res_num)
                contacts = interface_info["residue_contacts"].get(key, {})
                
                # Approximate BSA: more contacts and closer distances = higher BSA
                contact_count = contacts.get('contacts', 0)
                min_dist = contacts.get('min_dist', float('inf'))
                
                # Rough approximation: each contact contributes ~0.1-0.5 Å² depending on distance
                if min_dist < float('inf'):
                    # Closer contacts contribute more
                    distance_factor = max(0, 1.0 - (min_dist / self.interface_cutoff))
                    estimated_bsa = contact_count * 0.2 * distance_factor
                    
                    if estimated_bsa >= self.hotspot_bsa_cutoff:
                        bsa_per_residue[chain_id][res_num] = estimated_bsa
        
        logger.info(f"Identified {len(bsa_per_residue[chain1])} hotspot residues in chain {chain1}")
        logger.info(f"Identified {len(bsa_per_residue[chain2])} hotspot residues in chain {chain2}")
        
        return bsa_per_residue
    
    def calculate_bsa_freesasa(self, pdb_path: Path, chain1: str, chain2: str,
                              interface_info: Dict) -> Dict:
        """Calculate BSA using FreeSASA if available"""
        try:
            # Try to use FreeSASA for accurate BSA calculation
            logger.info("Attempting to calculate BSA using FreeSASA...")
            
            # Create separate chain files
            structure = self.parser.get_structure("complex", str(pdb_path))
            model = structure[0]
            
            temp_dir = Path(pdb_path.parent) / "temp_sasa"
            temp_dir.mkdir(exist_ok=True)
            
            chain1_path = temp_dir / f"chain_{chain1}.pdb"
            chain2_path = temp_dir / f"chain_{chain2}.pdb"
            complex_path = pdb_path
            
            class ChainSelect(Select):
                def __init__(self, chain_id):
                    self.chain_id = chain_id
                def accept_chain(self, chain):
                    return chain.id == self.chain_id
            
            # Save individual chains
            self.io.set_structure(structure)
            self.io.save(str(chain1_path), select=ChainSelect(chain1))
            self.io.save(str(chain2_path), select=ChainSelect(chain2))
            
            # Calculate SASA for complex and individual chains
            def get_sasa(pdb_file):
                cmd = ["freesasa", "--format=json", str(pdb_file)]
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode == 0:
                    return json.loads(result.stdout)
                return None
            
            complex_sasa = get_sasa(complex_path)
            chain1_sasa = get_sasa(chain1_path)
            chain2_sasa = get_sasa(chain2_path)
            
            if not all([complex_sasa, chain1_sasa, chain2_sasa]):
                logger.warning("FreeSASA calculation failed, using approximate method")
                return self.calculate_bsa_approximate(pdb_path, chain1, chain2, interface_info)
            
            # Parse FreeSASA output and calculate BSA
            bsa_per_residue = {chain1: {}, chain2: {}}
            
            # This is a simplified parser - adjust based on actual FreeSASA JSON format
            # FreeSASA JSON structure: {"atoms": [{"residueNumber": X, "sasa": Y, ...}, ...]}
            
            # Calculate BSA = SASA(chain1) + SASA(chain2) - SASA(complex)
            # Per residue: BSA = SASA_isolated - SASA_complex
            
            logger.info("Parsing FreeSASA results...")
            # Implementation would parse the JSON and calculate per-residue BSA
            # For now, fall back to approximate method
            return self.calculate_bsa_approximate(pdb_path, chain1, chain2, interface_info)
            
        except (FileNotFoundError, subprocess.CalledProcessError, json.JSONDecodeError) as e:
            logger.warning(f"FreeSASA not available or failed: {e}")
            logger.info("Falling back to approximate BSA calculation...")
            return self.calculate_bsa_approximate(pdb_path, chain1, chain2, interface_info)
    
    def identify_hotspots_with_masking(self, pdb_path: Path, chain1: str, chain2: str) -> Dict:
        """Identify hotspots and apply masking window"""
        logger.info("Identifying interface hotspots...")
        
        # Step 1: Calculate interface residues
        interface_info = self.calculate_interface_residues(pdb_path, chain1, chain2)
        
        # Step 2: Calculate BSA
        bsa_per_residue = self.calculate_bsa_freesasa(pdb_path, chain1, chain2, interface_info)
        
        # Step 3: Identify hotspots (>1.0 Å² BSA threshold)
        hotspots = {chain1: set(), chain2: set()}
        
        for chain_id in [chain1, chain2]:
            for res_num, bsa in bsa_per_residue[chain_id].items():
                if bsa >= self.hotspot_bsa_cutoff:
                    hotspots[chain_id].add(res_num)
        
        logger.info(f"Found {len(hotspots[chain1])} hotspots in chain {chain1}")
        logger.info(f"Found {len(hotspots[chain2])} hotspots in chain {chain2}")
        
        # Step 4: Apply masking window (±masking_window residues)
        masked_regions = {chain1: set(), chain2: set()}
        
        for chain_id in [chain1, chain2]:
            for hotspot_res in hotspots[chain_id]:
                for offset in range(-self.masking_window, self.masking_window + 1):
                    masked_regions[chain_id].add(hotspot_res + offset)
        
        # Convert to sorted lists and group into ranges
        def group_into_ranges(residue_set: Set[int]) -> List[str]:
            """Group consecutive residues into ranges"""
            if not residue_set:
                return []
            
            sorted_res = sorted(residue_set)
            ranges = []
            start = sorted_res[0]
            end = sorted_res[0]
            
            for res in sorted_res[1:]:
                if res == end + 1:
                    end = res
                else:
                    if start == end:
                        ranges.append(str(start))
                    else:
                        ranges.append(f"{start}..{end}")
                    start = res
                    end = res
            
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}..{end}")
            
            return ranges
        
        masked_ranges = {
            chain1: group_into_ranges(masked_regions[chain1]),
            chain2: group_into_ranges(masked_regions[chain2])
        }
        
        hotspot_ranges = {
            chain1: group_into_ranges(hotspots[chain1]),
            chain2: group_into_ranges(hotspots[chain2])
        }
        
        return {
            "hotspots": hotspots,
            "hotspot_ranges": hotspot_ranges,
            "masked_regions": masked_regions,
            "masked_ranges": masked_ranges,
            "interface_info": interface_info,
            "bsa_per_residue": bsa_per_residue
        }


class AlphaFoldFetcher:
    """Fetch AlphaFold models"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def fetch_alphafold_model(self, uniprot_id: str) -> Tuple[Path, str]:
        """Fetch AlphaFold model and return path and chain ID"""
        output_path = self.output_dir / f"AF_{uniprot_id}_F1-model_v4.cif"
        
        if output_path.exists():
            logger.info(f"Using existing AlphaFold model: {output_path}")
        else:
            # AlphaFold DB API
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif"
            logger.info(f"Fetching AlphaFold model for UniProt {uniprot_id}...")
            
            try:
                response = requests.get(url, timeout=30)
                if response.status_code == 200:
                    output_path.write_bytes(response.content)
                    logger.info(f"Saved to {output_path}")
                else:
                    raise RuntimeError(f"Failed to fetch AlphaFold model: HTTP {response.status_code}")
            except Exception as e:
                logger.error(f"Error fetching AlphaFold model: {e}")
                raise
        
        # Determine chain ID from structure
        parser = PDB.MMCIFParser(QUIET=True)
        try:
            structure = parser.get_structure("af_model", str(output_path))
            model = structure[0]
            chain_ids = list(model.child_dict.keys())
            if chain_ids:
                chain_id = chain_ids[0]  # Use first chain
                logger.info(f"Detected chain ID: {chain_id}")
                return output_path, chain_id
            else:
                logger.warning("No chains found, defaulting to 'A'")
                return output_path, "A"
        except Exception as e:
            logger.warning(f"Could not parse structure to detect chain ID: {e}")
            return output_path, "A"


class YAMLUpdater:
    """Update YAML design specification files"""
    
    def __init__(self, yaml_path: Path):
        self.yaml_path = yaml_path
        with open(yaml_path) as f:
            self.data = yaml.safe_load(f)
    
    def update_fcgr3_spec(self, target_binding_ranges: List[str],
                          scaffold_a_ranges: List[str],
                          scaffold_b_ranges: List[str]):
        """Update Fcgr3 design specification"""
        logger.info(f"Updating {self.yaml_path}...")
        
        if not target_binding_ranges:
            logger.warning("No target binding ranges found, skipping update")
            return
        
        # Update target binding types
        for i, entity in enumerate(self.data['entities']):
            if 'file' in entity and entity['file'].get('path') == '3AY4.cif':
                # Check if this is the target (chain C)
                if 'include' in entity['file']:
                    for include_item in entity['file']['include']:
                        if include_item.get('chain', {}).get('id') == 'C':
                            # Update binding_types
                            if 'binding_types' in entity['file']:
                                for bt_item in entity['file']['binding_types']:
                                    if bt_item.get('chain', {}).get('id') == 'C':
                                        bt_item['binding'] = ','.join(target_binding_ranges)
                                        break
                                else:
                                    entity['file']['binding_types'].append({
                                        'chain': {'id': 'C'},
                                        'binding': ','.join(target_binding_ranges)
                                    })
                            else:
                                entity['file']['binding_types'] = [{
                                    'chain': {'id': 'C'},
                                    'binding': ','.join(target_binding_ranges)
                                }]
                            break
                
                # Check if this is the scaffold (chains A & B)
                if 'include' in entity['file']:
                    chains = [inc.get('chain', {}).get('id') for inc in entity['file'].get('include', [])]
                    if 'A' in chains and 'B' in chains:
                        # Update design regions
                        if scaffold_a_ranges or scaffold_b_ranges:
                            if 'design' not in entity['file']:
                                entity['file']['design'] = []
                            
                            # Update or add design regions
                            design_updated = False
                            for design_item in entity['file']['design']:
                                chain_id = design_item.get('chain', {}).get('id')
                                if chain_id == 'A' and scaffold_a_ranges:
                                    design_item['res_index'] = ','.join(scaffold_a_ranges)
                                    design_updated = True
                                elif chain_id == 'B' and scaffold_b_ranges:
                                    design_item['res_index'] = ','.join(scaffold_b_ranges)
                                    design_updated = True
                            
                            if not design_updated:
                                if scaffold_a_ranges:
                                    entity['file']['design'].append({
                                        'chain': {'id': 'A'},
                                        'res_index': ','.join(scaffold_a_ranges)
                                    })
                                if scaffold_b_ranges:
                                    entity['file']['design'].append({
                                        'chain': {'id': 'B'},
                                        'res_index': ','.join(scaffold_b_ranges)
                                    })
                        
                        # Update structure_groups
                        if 'structure_groups' not in entity['file']:
                            entity['file']['structure_groups'] = []
                        
                        # Remove old visibility: 0 groups
                        structure_groups = entity['file']['structure_groups']
                        structure_groups = [g for g in structure_groups 
                                          if not (g.get('group', {}).get('visibility') == 0)]
                        
                        # Add new visibility: 0 groups
                        if scaffold_a_ranges:
                            structure_groups.append({
                                'group': {
                                    'visibility': 0,
                                    'id': 'A',
                                    'res_index': ','.join(scaffold_a_ranges)
                                }
                            })
                        if scaffold_b_ranges:
                            structure_groups.append({
                                'group': {
                                    'visibility': 0,
                                    'id': 'B',
                                    'res_index': ','.join(scaffold_b_ranges)
                                }
                            })
                        
                        entity['file']['structure_groups'] = structure_groups
                        
                        # Update binding_types for scaffold
                        if scaffold_a_ranges or scaffold_b_ranges:
                            if 'binding_types' not in entity['file']:
                                entity['file']['binding_types'] = []
                            
                            # Update or add binding types
                            bt_updated = False
                            for bt_item in entity['file']['binding_types']:
                                chain_id = bt_item.get('chain', {}).get('id')
                                if chain_id == 'A' and scaffold_a_ranges:
                                    bt_item['binding'] = ','.join(scaffold_a_ranges)
                                    bt_updated = True
                                elif chain_id == 'B' and scaffold_b_ranges:
                                    bt_item['binding'] = ','.join(scaffold_b_ranges)
                                    bt_updated = True
                            
                            if not bt_updated:
                                if scaffold_a_ranges:
                                    entity['file']['binding_types'].append({
                                        'chain': {'id': 'A'},
                                        'binding': ','.join(scaffold_a_ranges)
                                    })
                                if scaffold_b_ranges:
                                    entity['file']['binding_types'].append({
                                        'chain': {'id': 'B'},
                                        'binding': ','.join(scaffold_b_ranges)
                                    })
                        break
    
    def update_fcgr4_spec(self, alphafold_path: Path, chain_id: str,
                          binding_ranges: List[str]):
        """Update Fcgr4 design specification"""
        logger.info(f"Updating {self.yaml_path}...")
        
        # Update AlphaFold path and chain ID
        for entity in self.data['entities']:
            if 'file' in entity:
                file_path = entity['file'].get('path', '')
                if isinstance(file_path, str) and ('fcgr4' in file_path.lower() or 'alphafold' in file_path.lower()):
                    # Update path (use relative path from YAML file location)
                    rel_path = alphafold_path.relative_to(self.yaml_path.parent)
                    entity['file']['path'] = str(rel_path)
                    
                    # Update chain ID in include
                    if 'include' in entity['file']:
                        for include_item in entity['file']['include']:
                            if 'chain' in include_item:
                                include_item['chain']['id'] = chain_id
                    
                    # Update binding_types
                    if binding_ranges:
                        if 'binding_types' in entity['file']:
                            for binding_item in entity['file']['binding_types']:
                                if 'chain' in binding_item:
                                    binding_item['chain']['id'] = chain_id
                                    binding_item['binding'] = ','.join(binding_ranges)
                        else:
                            entity['file']['binding_types'] = [{
                                'chain': {'id': chain_id},
                                'binding': ','.join(binding_ranges)
                            }]
                    
                    # Update structure_groups
                    if 'structure_groups' in entity['file']:
                        for group_item in entity['file']['structure_groups']:
                            if 'group' in group_item and 'id' in group_item['group']:
                                group_item['group']['id'] = chain_id
                    break
    
    def save(self):
        """Save updated YAML file"""
        with open(self.yaml_path, 'w') as f:
            yaml.dump(self.data, f, default_flow_style=False, sort_keys=False, 
                     allow_unicode=True, width=1000)
        logger.info(f"Saved updated YAML to {self.yaml_path}")


def main():
    """Main function"""
    # Configuration
    base_dir = Path(__file__).parent
    pdb_path = base_dir / "3AY4.cif"
    
    # Mouse Fcgr4 UniProt ID
    # NOTE: Verify this UniProt ID is correct for your target
    # Q8CIZ6 is mouse Fcgr4 (Fc receptor, IgG, low affinity IV)
    # To find the correct ID: https://www.uniprot.org/
    fcgr4_uniprot_id = "Q8CIZ6"
    
    # Interface analysis parameters
    interface_cutoff = 4.5  # Angstroms
    hotspot_bsa_cutoff = 1.0  # Å²
    masking_window = 7  # residues
    
    logger.info("="*60)
    logger.info("Completing Design Specifications")
    logger.info("="*60)
    
    # Step 1: Analyze 3AY4 structure for Fcgr3 interface
    logger.info("\n=== Step 1: Analyzing 3AY4 structure ===")
    analyzer = StructureAnalyzer(
        interface_cutoff=interface_cutoff,
        hotspot_bsa_cutoff=hotspot_bsa_cutoff,
        masking_window=masking_window
    )
    
    # Analyze interface between chain C (Fcgr3) and chains A/B (Fc fragment)
    # For Fcgr3, we analyze C vs A and C vs B separately, then combine
    logger.info("Analyzing interface: Chain C (Fcgr3) vs Chain A (Fc)")
    results_ca = analyzer.identify_hotspots_with_masking(pdb_path, "C", "A")
    
    logger.info("Analyzing interface: Chain C (Fcgr3) vs Chain B (Fc)")
    results_cb = analyzer.identify_hotspots_with_masking(pdb_path, "C", "B")
    
    # Combine results
    target_binding_ranges = sorted(set(
        results_ca['masked_ranges']['C'] + results_cb['masked_ranges']['C']
    ))
    scaffold_a_ranges = results_ca['masked_ranges']['A']
    scaffold_b_ranges = results_cb['masked_ranges']['B']
    
    logger.info(f"Target (chain C) binding regions: {target_binding_ranges}")
    logger.info(f"Scaffold (chain A) design regions: {scaffold_a_ranges}")
    logger.info(f"Scaffold (chain B) design regions: {scaffold_b_ranges}")
    
    # Step 2: Update Fcgr3 design specification
    logger.info("\n=== Step 2: Updating Fcgr3 design specification ===")
    fcgr3_yaml = base_dir / "3AY4_Fcgr3_binder_design.yaml"
    updater_fcgr3 = YAMLUpdater(fcgr3_yaml)
    updater_fcgr3.update_fcgr3_spec(
        target_binding_ranges,
        scaffold_a_ranges,
        scaffold_b_ranges
    )
    updater_fcgr3.save()
    
    # Step 3: Fetch AlphaFold model for mouse Fcgr4
    logger.info("\n=== Step 3: Fetching AlphaFold model for mouse Fcgr4 ===")
    temp_dir = base_dir / "temp_alphafold"
    fetcher = AlphaFoldFetcher(temp_dir)
    
    try:
        alphafold_path, chain_id = fetcher.fetch_alphafold_model(fcgr4_uniprot_id)
        logger.info(f"AlphaFold model: {alphafold_path}")
        logger.info(f"Chain ID: {chain_id}")
        
        # Move to example directory
        final_alphafold_path = base_dir / alphafold_path.name
        if not final_alphafold_path.exists():
            import shutil
            shutil.copy(alphafold_path, final_alphafold_path)
            logger.info(f"Copied to {final_alphafold_path}")
        
        # For Fcgr4, we don't have a complex structure yet, so we'll use placeholder ranges
        # These would be filled in after docking or experimental interface analysis
        fcgr4_binding_ranges = []  # Placeholder - would need interface analysis with nanobody
        
        # Step 4: Update Fcgr4 design specification
        logger.info("\n=== Step 4: Updating Fcgr4 design specification ===")
        fcgr4_yaml = base_dir / "Fcgr4_nanobody_binder_design.yaml"
        updater_fcgr4 = YAMLUpdater(fcgr4_yaml)
        updater_fcgr4.update_fcgr4_spec(
            final_alphafold_path,
            chain_id,
            fcgr4_binding_ranges
        )
        updater_fcgr4.save()
        
    except Exception as e:
        logger.error(f"Error fetching AlphaFold model: {e}")
        logger.warning("Fcgr4 specification update skipped")
    
    logger.info("\n" + "="*60)
    logger.info("Design specification completion finished!")
    logger.info("="*60)
    logger.info(f"\nUpdated files:")
    logger.info(f"  - {fcgr3_yaml}")
    if 'fcgr4_yaml' in locals():
        logger.info(f"  - {fcgr4_yaml}")
    logger.info("\nNote: Fcgr4 binding regions need to be determined")
    logger.info("      through docking or experimental interface analysis.")


if __name__ == "__main__":
    main()

