#!/usr/bin/env python3
"""
Score Boltzgen candidates using Boltz-2 predictions + IPSAE.

This script:
1. Extracts sequences from Boltzgen CIF files
2. Creates YAML input for Boltz-2
3. Runs Boltz-2 predictions to get full PAE matrices
4. Runs IPSAE on Boltz predictions to get ipSAE and pDockQ scores
5. Compares with Boltzgen native metrics

Usage:
    # Test on 1 candidate
    python score_with_boltz_ipsae.py --test-one

    # Run on top N candidates
    python score_with_boltz_ipsae.py --top-n 5

    # Run on all candidates
    python score_with_boltz_ipsae.py --all
"""

import os
import sys
import re
import subprocess
import tempfile
import shutil
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
from tqdm import tqdm

# --- CONFIGURATION ---
IPSAE_SCRIPT = "/lambda/nfs/Nipah-hackathon/andrei/IPSAE/ipsae.py"
RESULT_DIR = "/lambda/nfs/Nipah-hackathon/andrei/boltzgen/protein_competition/results/20designs_5budget"
BOLTZGEN_RANKINGS = "./final_rankings_boltzgen.csv"

# Boltz settings
PAE_CUTOFF = 10
DIST_CUTOFF = 10
BOLTZ_DIFFUSION_SAMPLES = 1
BOLTZ_TIMEOUT = 600  # 10 minutes per prediction


def extract_sequences_from_cif(cif_path: str) -> Dict[str, str]:
    """
    Extract protein sequences from CIF file.

    Returns:
        Dict mapping chain_id -> sequence
    """
    sequences = {}
    current_chain = None

    with open(cif_path, 'r') as f:
        in_poly_section = False

        for line in f:
            # Find the entity_poly section which has sequences
            if line.startswith('_entity_poly.'):
                in_poly_section = True
                continue

            if in_poly_section:
                if line.startswith('#') or line.strip() == '':
                    in_poly_section = False
                    continue

                # Parse lines like: "1 polypeptide(L) ? MKVL..."
                parts = line.split(None, 3)
                if len(parts) >= 4 and 'polypeptide' in parts[1]:
                    entity_id = parts[0]
                    sequence = parts[3].strip()

                    # Map entity to chain (simple: 1->A, 2->B, etc.)
                    # For more complex cases, we'd parse _struct_asym
                    chain_id = chr(ord('A') + int(entity_id) - 1)
                    sequences[chain_id] = sequence

    if not sequences:
        # Fallback: try to parse from _entity_poly section more carefully
        with open(cif_path, 'r') as f:
            content = f.read()
            # Look for sequences in different format
            seq_pattern = r'polypeptide\(L\)\s+\?\s+([A-Z]+)'
            matches = re.findall(seq_pattern, content)
            for i, seq in enumerate(matches):
                chain_id = chr(ord('A') + i)
                sequences[chain_id] = seq

    return sequences


def create_boltz_yaml(sequences: Dict[str, str], output_path: str) -> bool:
    """
    Create YAML input file for Boltz-2.

    Args:
        sequences: Dict mapping chain_id -> sequence
        output_path: Where to save YAML file

    Returns:
        True if successful
    """
    try:
        yaml_content = "version: 1\nsequences:\n"

        for chain_id, sequence in sorted(sequences.items()):
            yaml_content += f"  - protein:\n"
            yaml_content += f"      id: {chain_id}\n"
            yaml_content += f"      sequence: {sequence}\n"

        with open(output_path, 'w') as f:
            f.write(yaml_content)

        return True
    except Exception as e:
        print(f"  ❌ Error creating YAML: {e}")
        return False


def run_boltz_prediction(yaml_path: str, output_dir: str, design_id: str) -> Optional[Tuple[str, str]]:
    """
    Run Boltz-2 prediction.

    Returns:
        (pae_npz_path, cif_path) or None if failed
    """
    try:
        cmd = [
            "boltz", "predict",
            str(yaml_path),
            "--out_dir", str(output_dir),
            "--diffusion_samples", str(BOLTZ_DIFFUSION_SAMPLES),
            "--num_workers", "1",
            "--accelerator", "gpu",
            "--devices", "1",
            "--use_msa_server"  # Generate MSAs automatically
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=BOLTZ_TIMEOUT
        )

        if result.returncode != 0:
            print(f"  ❌ Boltz failed: {result.stderr[:200]}")
            return None

        # Find output files
        # Boltz structure: output_dir/boltz_results_<job_name>/predictions/<job_name>/
        pred_base = Path(output_dir)

        # Find boltz_results directory
        boltz_result_dirs = list(pred_base.glob("boltz_results_*"))
        if not boltz_result_dirs:
            print(f"  ❌ No boltz_results directory found in {pred_base}")
            return None

        boltz_dir = boltz_result_dirs[0]
        pred_subdir = boltz_dir / "predictions"

        if not pred_subdir.exists():
            print(f"  ❌ No predictions subdirectory in {boltz_dir}")
            return None

        # Find the job directory inside predictions
        job_dirs = [d for d in pred_subdir.iterdir() if d.is_dir()]
        if not job_dirs:
            print(f"  ❌ No job directory found in {pred_subdir}")
            return None

        pred_dir = job_dirs[0]

        # Find PAE and CIF files
        pae_files = list(pred_dir.glob("pae_*_model_0.npz"))
        cif_files = list(pred_dir.glob("*_model_0.cif"))
        # Exclude pae, plddt, pde files from cif search
        cif_files = [f for f in cif_files if not any(x in f.name for x in ['pae_', 'plddt_', 'pde_'])]

        if not pae_files or not cif_files:
            print(f"  ❌ Missing output files in {pred_dir}")
            print(f"     PAE files: {pae_files}")
            print(f"     CIF files: {cif_files}")
            return None

        return str(pae_files[0]), str(cif_files[0])

    except subprocess.TimeoutExpired:
        print(f"  ❌ Boltz prediction timeout ({BOLTZ_TIMEOUT}s)")
        return None
    except Exception as e:
        print(f"  ❌ Boltz prediction error: {e}")
        return None


def parse_ipsae_output(output_text: str) -> Dict[str, float]:
    """
    Extract ipSAE and pDockQ scores from IPSAE output file.

    Returns:
        Dict with 'ipSAE' and 'pDockQ' keys
    """
    scores = {'ipSAE': 0.0, 'pDockQ': 0.0}

    # Look for the "max" row which contains the best scores
    # Format: Chn1 Chn2 PAE Dist Type ipSAE ... pDockQ ...
    # Example: A    B    10  10  max  0.788510 ... 0.3463 ...
    lines = output_text.split('\n')
    for line in lines:
        if 'max' in line and len(line.split()) >= 11:
            parts = line.split()
            try:
                # ipSAE is column 5 (0-indexed)
                scores['ipSAE'] = float(parts[5])
                # pDockQ is column 10 (0-indexed)
                scores['pDockQ'] = float(parts[10])
                break  # We found the max row, we're done
            except (ValueError, IndexError):
                pass

    return scores


def run_ipsae(pae_path: str, cif_path: str) -> Optional[Dict[str, float]]:
    """
    Run IPSAE scoring.

    Returns:
        Dict with scores or None if failed
    """
    try:
        cmd = [
            "python", IPSAE_SCRIPT,
            pae_path,
            cif_path,
            str(PAE_CUTOFF),
            str(DIST_CUTOFF)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        if result.returncode != 0:
            print(f"  ❌ IPSAE failed: {result.stderr[:200]}")
            return None

        # IPSAE writes output to a file, not stdout
        # Find the output file: <cif_basename>_<pae>_<dist>.txt
        cif_dir = Path(cif_path).parent
        cif_base = Path(cif_path).stem  # e.g., "job_A_medium_00_model_0"
        output_file = cif_dir / f"{cif_base}_{PAE_CUTOFF}_{DIST_CUTOFF}.txt"

        if not output_file.exists():
            print(f"  ❌ IPSAE output file not found: {output_file}")
            return None

        # Read and parse the output file
        with open(output_file, 'r') as f:
            output_text = f.read()

        scores = parse_ipsae_output(output_text)
        return scores

    except Exception as e:
        print(f"  ❌ IPSAE error: {e}")
        return None


def score_single_candidate(
    design_id: str,
    cif_path: str,
    temp_dir: str
) -> Optional[Dict]:
    """
    Score a single candidate through the full pipeline.

    Returns:
        Dict with all scores or None if failed
    """
    print(f"\n{'='*80}")
    print(f"Scoring: {design_id}")
    print(f"{'='*80}")

    # Step 1: Extract sequences
    print("[1/4] Extracting sequences from CIF...")
    sequences = extract_sequences_from_cif(cif_path)

    if not sequences:
        print(f"  ❌ Failed to extract sequences")
        return None

    print(f"  ✓ Found {len(sequences)} chains: {list(sequences.keys())}")
    for chain, seq in sequences.items():
        print(f"    Chain {chain}: {len(seq)} residues")

    # Step 2: Create Boltz YAML
    print("[2/4] Creating Boltz input YAML...")
    yaml_path = Path(temp_dir) / "inputs" / f"{design_id}.yaml"
    yaml_path.parent.mkdir(parents=True, exist_ok=True)

    if not create_boltz_yaml(sequences, str(yaml_path)):
        return None

    print(f"  ✓ Created {yaml_path}")

    # Step 3: Run Boltz prediction
    print("[3/4] Running Boltz-2 prediction (this may take 3-5 minutes)...")
    boltz_output_dir = Path(temp_dir) / "predictions"
    boltz_output_dir.mkdir(parents=True, exist_ok=True)

    result = run_boltz_prediction(str(yaml_path), str(boltz_output_dir), design_id)

    if result is None:
        return None

    pae_path, pred_cif_path = result
    print(f"  ✓ Boltz prediction complete")
    print(f"    PAE: {Path(pae_path).name}")
    print(f"    CIF: {Path(pred_cif_path).name}")

    # Step 4: Run IPSAE
    print("[4/4] Running IPSAE scoring...")
    ipsae_scores = run_ipsae(pae_path, pred_cif_path)

    if ipsae_scores is None:
        return None

    print(f"  ✓ IPSAE complete")
    print(f"    ipSAE:  {ipsae_scores['ipSAE']:.4f}")
    print(f"    pDockQ: {ipsae_scores['pDockQ']:.4f}")

    return {
        'design_id': design_id,
        'boltz_ipSAE': ipsae_scores['ipSAE'],
        'boltz_pDockQ': ipsae_scores['pDockQ'],
        'original_cif': cif_path,
        'boltz_cif': pred_cif_path,
        'boltz_pae': pae_path
    }


def load_boltzgen_rankings() -> pd.DataFrame:
    """Load Boltzgen rankings to get top candidates."""
    if not os.path.exists(BOLTZGEN_RANKINGS):
        print(f"Warning: {BOLTZGEN_RANKINGS} not found")
        return pd.DataFrame()

    return pd.read_csv(BOLTZGEN_RANKINGS)


def get_cif_path(design_id: str) -> Optional[str]:
    """Get CIF path for a design ID."""
    cif_path = Path(RESULT_DIR) / "intermediate_designs_inverse_folded" / "refold_cif" / f"{design_id}.cif"

    if cif_path.exists():
        return str(cif_path)

    return None


def main():
    parser = argparse.ArgumentParser(description="Score candidates with Boltz-2 + IPSAE")
    parser.add_argument('--test-one', action='store_true', help='Test on one candidate')
    parser.add_argument('--top-n', type=int, help='Score top N candidates from Boltzgen rankings')
    parser.add_argument('--all', action='store_true', help='Score all candidates')
    parser.add_argument('--design-id', type=str, help='Score specific design by ID')

    args = parser.parse_args()

    print("="*80)
    print("Boltz-2 + IPSAE Scoring Pipeline")
    print("="*80)

    # Check IPSAE script exists
    if not os.path.exists(IPSAE_SCRIPT):
        print(f"❌ IPSAE script not found: {IPSAE_SCRIPT}")
        sys.exit(1)

    # Determine which candidates to score
    candidates_to_score = []

    if args.design_id:
        cif_path = get_cif_path(args.design_id)
        if cif_path:
            candidates_to_score = [(args.design_id, cif_path)]
        else:
            print(f"❌ CIF not found for {args.design_id}")
            sys.exit(1)

    elif args.test_one:
        print("\n🧪 Test Mode: Scoring 1 candidate")
        # Get top candidate from Boltzgen rankings
        rankings = load_boltzgen_rankings()
        if rankings.empty:
            print("❌ No rankings found, using job_A_medium_00")
            design_id = "job_A_medium_00"
        else:
            design_id = rankings.iloc[0]['id']

        cif_path = get_cif_path(design_id)
        if cif_path:
            candidates_to_score = [(design_id, cif_path)]
        else:
            print(f"❌ CIF not found for {design_id}")
            sys.exit(1)

    elif args.top_n:
        print(f"\n📊 Scoring top {args.top_n} candidates")
        rankings = load_boltzgen_rankings()
        if rankings.empty:
            print("❌ No rankings found")
            sys.exit(1)

        for i in range(min(args.top_n, len(rankings))):
            design_id = rankings.iloc[i]['id']
            cif_path = get_cif_path(design_id)
            if cif_path:
                candidates_to_score.append((design_id, cif_path))

    elif args.all:
        print("\n📊 Scoring ALL candidates")
        rankings = load_boltzgen_rankings()
        for _, row in rankings.iterrows():
            design_id = row['id']
            cif_path = get_cif_path(design_id)
            if cif_path:
                candidates_to_score.append((design_id, cif_path))

    else:
        print("❌ Please specify --test-one, --top-n N, --all, or --design-id ID")
        sys.exit(1)

    print(f"\n📝 Will score {len(candidates_to_score)} candidates")
    print(f"⏱️  Estimated time: {len(candidates_to_score) * 4} minutes\n")

    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix="boltz_ipsae_")
    print(f"📁 Temp directory: {temp_dir}\n")

    results = []

    try:
        # Score each candidate
        for design_id, cif_path in candidates_to_score:
            result = score_single_candidate(design_id, cif_path, temp_dir)
            if result:
                results.append(result)
                print(f"\n✅ Successfully scored {design_id}")
            else:
                print(f"\n❌ Failed to score {design_id}")

        # Save results
        if results:
            output_csv = "boltz_ipsae_scores.csv"
            df = pd.DataFrame(results)
            df.to_csv(output_csv, index=False)

            print(f"\n{'='*80}")
            print(f"✅ SUCCESS! Scored {len(results)}/{len(candidates_to_score)} candidates")
            print(f"{'='*80}")
            print(f"\n📄 Results saved to: {os.path.abspath(output_csv)}")

            # Display results
            print("\n📊 Results:")
            print(df[['design_id', 'boltz_ipSAE', 'boltz_pDockQ']].to_string(index=False))

            # Compare with Boltzgen if available
            if os.path.exists(BOLTZGEN_RANKINGS):
                print(f"\n🔄 Comparing with Boltzgen metrics...")
                boltzgen_df = pd.read_csv(BOLTZGEN_RANKINGS)

                comparison = df.merge(
                    boltzgen_df[['id', 'design_to_target_iptm', 'min_interaction_pae', 'design_ptm']],
                    left_on='design_id',
                    right_on='id',
                    how='left'
                )

                comparison = comparison[[
                    'design_id',
                    'boltz_ipSAE',
                    'design_to_target_iptm',
                    'boltz_pDockQ',
                    'design_ptm'
                ]]

                comparison.columns = [
                    'ID',
                    'Boltz ipSAE',
                    'Boltzgen ipTM',
                    'Boltz pDockQ',
                    'Boltzgen PTM'
                ]

                print(comparison.to_string(index=False))

                # Save comparison
                comparison.to_csv("boltz_boltzgen_comparison.csv", index=False)
                print(f"\n📄 Comparison saved to: boltz_boltzgen_comparison.csv")

        else:
            print("\n❌ No candidates were successfully scored")

    finally:
        # Cleanup (optional - keep for debugging)
        try:
            keep_temp = input("\n🗑️  Keep temporary files for inspection? (y/N): ").lower() == 'y'
        except EOFError:
            # Running in background or non-interactive
            keep_temp = True
            print(f"\n📁 Temporary files kept in: {temp_dir}")

        if not keep_temp:
            print(f"Cleaning up {temp_dir}...")
            shutil.rmtree(temp_dir, ignore_errors=True)
        else:
            print(f"Temporary files kept in: {temp_dir}")


if __name__ == "__main__":
    main()
