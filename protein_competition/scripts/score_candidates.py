import os
import subprocess
import pandas as pd
import glob
import re
import sys
import tempfile
import shutil
from pathlib import Path

# --- CONFIGURATION ---
IPSAE_SCRIPT = "/lambda/nfs/Nipah-hackathon/andrei/IPSAE/ipsae.py"
BOLTZ_PREDICT_SCRIPT = "/lambda/nfs/Nipah-hackathon/andrei/boltz/boltz"  # Adjust if needed

# List of your experiment directories to score
# Example: ["results_A_medium", "results_B_small", "20designs_5budget"]
RESULT_DIRS = [
    # "results_A_medium",
    # "results_B_small",
    # "results_C_large",
    "/lambda/nfs/Nipah-hackathon/andrei/boltzgen/protein_competition/results/20designs_5budget"  # Added your example dir
]

OUTPUT_CSV = "final_rankings.csv"
PAE_CUTOFF = 10
DIST_CUTOFF = 10
USE_BOLTZ_PREDICTION = True  # Set to True to run Boltz-2 predictions for full PAE matrices


def parse_ipsae_output(output_text):
    """Extracts scores from IPSAE stdout."""
    scores = {'ipSAE': 0.0, 'pDockQ': 0.0}
    # Regex for "ipSAE: 0.85" or "ipSAE:0.85"
    ipsae_match = re.search(r"ipSAE[:\s]+([\d\.]+)", output_text)
    pdockq_match = re.search(r"pDockQ[:\s]+([\d\.]+)", output_text)

    if ipsae_match: scores['ipSAE'] = float(ipsae_match.group(1))
    if pdockq_match: scores['pDockQ'] = float(pdockq_match.group(1))
    return scores


def run_boltz_prediction(cif_path, output_dir):
    """
    Run Boltz-2 prediction on a CIF file to get full PAE matrix.
    Returns the path to the PAE file (npz format) and CIF file.
    """
    try:
        cif_name = Path(cif_path).stem
        # Run boltz predict
        cmd = [
            "boltz", "predict",
            str(cif_path),
            "--out_dir", str(output_dir),
            "--num_workers", "1",
            "--accelerator", "gpu",
            "--devices", "1"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            print(f"  Warning: Boltz prediction failed for {cif_name}: {result.stderr[:200]}")
            return None, None

        # Find the generated pae and confidence files
        # Boltz-2 output: predictions/<job_name>/pae_<job_name>_model_0.npz
        pred_dir = Path(output_dir) / "predictions" / cif_name
        if not pred_dir.exists():
            # Try alternative naming
            pred_dirs = list(Path(output_dir).glob("predictions/*"))
            if pred_dirs:
                pred_dir = pred_dirs[0]
            else:
                print(f"  Warning: Could not find predictions directory for {cif_name}")
                return None, None

        pae_files = list(pred_dir.glob("pae_*_model_0.npz"))
        cif_files = list(pred_dir.glob("*_model_0.cif"))

        if pae_files and cif_files:
            return str(pae_files[0]), str(cif_files[0])
        else:
            print(f"  Warning: Could not find PAE/CIF files in {pred_dir}")
            return None, None

    except subprocess.TimeoutExpired:
        print(f"  Warning: Boltz prediction timeout for {cif_path}")
        return None, None
    except Exception as e:
        print(f"  Warning: Boltz prediction error for {cif_path}: {e}")
        return None, None


def find_cif_npz_pairs(root_dir):
    """
    Smartly finds matched CIF and NPZ files based on the directory structure.
    Prioritizes Refolded/Inverse Folded results.
    """
    pairs = []  # List of (design_id, cif_path, pae_path)

    # 1. Strategy A: Inverse Folded (The Gold Standard)
    # Structure: root/intermediate_designs_inverse_folded/refold_cif/*.cif
    #            root/intermediate_designs_inverse_folded/fold_out_npz/*.npz
    inv_dir = os.path.join(root_dir, "intermediate_designs_inverse_folded")
    if os.path.exists(inv_dir):
        print(f"  > Detected Inverse Folded results in {os.path.basename(root_dir)}")
        cif_dir = os.path.join(inv_dir, "refold_cif")
        npz_dir = os.path.join(inv_dir, "fold_out_npz")  # Usually fold_out_npz or fold_out_design_npz

        # Fallback if names differ slightly
        if not os.path.exists(npz_dir):
            npz_dir = os.path.join(inv_dir, "fold_out_design_npz")

        if os.path.exists(cif_dir) and os.path.exists(npz_dir):
            cif_files = glob.glob(os.path.join(cif_dir, "*.cif"))
            for cif_path in cif_files:
                filename = os.path.basename(cif_path)
                design_id = os.path.splitext(filename)[0]  # remove .cif

                # Try to find matching NPZ
                npz_path = os.path.join(npz_dir, f"{design_id}.npz")

                if os.path.exists(npz_path):
                    pairs.append((design_id, cif_path, npz_path))

            if pairs: return pairs

    # 2. Strategy B: Intermediate Designs (Raw Generation)
    # Structure: root/intermediate_designs/X.cif and X.npz (Flat)
    inter_dir = os.path.join(root_dir, "intermediate_designs")
    if os.path.exists(inter_dir):
        print(f"  > Detected Flat Intermediate Designs in {os.path.basename(root_dir)}")
        cif_files = glob.glob(os.path.join(inter_dir, "*.cif"))
        for cif_path in cif_files:
            filename = os.path.basename(cif_path)
            design_id = os.path.splitext(filename)[0]
            npz_path = os.path.join(inter_dir, f"{design_id}.npz")

            if os.path.exists(npz_path):
                pairs.append((design_id, cif_path, npz_path))

        if pairs: return pairs

    # 3. Strategy C: Standard Boltz/Predictions Folder (Nested)
    # Structure: root/predictions/design_id/X.cif
    pred_dir = os.path.join(root_dir, "predictions")
    if os.path.exists(pred_dir):
        print(f"  > Detected Nested Predictions in {os.path.basename(root_dir)}")
        for design_folder in os.listdir(pred_dir):
            folder_path = os.path.join(pred_dir, design_folder)
            if not os.path.isdir(folder_path): continue

            cif_files = glob.glob(os.path.join(folder_path, "*_model_0.cif"))
            pae_files = glob.glob(os.path.join(folder_path, "pae_*_model_0.npz"))

            if cif_files and pae_files:
                pairs.append((design_folder, cif_files[0], pae_files[0]))

        if pairs: return pairs

    return []


def score_directory(dir_path):
    results_list = []
    if not os.path.exists(dir_path):
        print(f"Warning: Directory not found: {dir_path}")
        return []

    print(f"--- Scanning {dir_path} ---")

    # Use smart finder
    pairs = find_cif_npz_pairs(dir_path)
    print(f"  > Found {len(pairs)} candidate pairs to score.")

    # Create temporary directory for Boltz predictions if needed
    temp_dir = None
    if USE_BOLTZ_PREDICTION:
        temp_dir = tempfile.mkdtemp(prefix="boltz_predictions_")
        print(f"  > Using Boltz-2 predictions. Temp dir: {temp_dir}")

    count = 0
    try:
        for design_id, cif_path, pae_path in pairs:
            actual_pae_path = pae_path
            actual_cif_path = cif_path

            # If USE_BOLTZ_PREDICTION, run Boltz-2 on the CIF to get full PAE
            if USE_BOLTZ_PREDICTION:
                print(f"  Running Boltz-2 prediction for {design_id}...")
                pae_file, pred_cif = run_boltz_prediction(cif_path, temp_dir)
                if pae_file is None or pred_cif is None:
                    print(f"  Skipping {design_id} - Boltz prediction failed")
                    continue
                actual_pae_path = pae_file
                actual_cif_path = pred_cif

            # Command: python ipsae.py <pae> <cif> <pae_cut> <dist_cut>
            cmd = [
                "python", IPSAE_SCRIPT,
                actual_pae_path,
                actual_cif_path,
                str(PAE_CUTOFF),
                str(DIST_CUTOFF)
            ]

            try:
                process = subprocess.run(cmd, capture_output=True, text=True)
                if process.returncode != 0:
                    print(f"Error on {design_id}: {process.stderr.strip()}")
                    continue

                parsed = parse_ipsae_output(process.stdout.strip())

                results_list.append({
                    "design_id": design_id,
                    "source_experiment": os.path.basename(dir_path),
                    "ipSAE": parsed['ipSAE'],
                    "pDockQ": parsed['pDockQ'],
                    "cif_path": cif_path,
                    "type": "refolded" if "refold" in cif_path else "raw"
                })

                count += 1
                if count % 10 == 0: print(f"  Scored {count}/{len(pairs)}...")

            except Exception as e:
                print(f"Exception {design_id}: {e}")

    finally:
        # Clean up temporary directory
        if temp_dir and os.path.exists(temp_dir):
            print(f"  > Cleaning up temporary directory: {temp_dir}")
            shutil.rmtree(temp_dir, ignore_errors=True)

    return results_list


def main():
    if not os.path.exists(IPSAE_SCRIPT):
        print(f"CRITICAL ERROR: IPSAE script not found at {IPSAE_SCRIPT}")
        sys.exit(1)

    all_scores = []

    print("Starting Scoring Pipeline...")
    for d in RESULT_DIRS:
        scores = score_directory(d)
        all_scores.extend(scores)

    if all_scores:
        df = pd.DataFrame(all_scores)
        # Sort by ipSAE Descending
        df_sorted = df.sort_values(by="ipSAE", ascending=False)
        df_sorted.to_csv(OUTPUT_CSV, index=False)

        print(f"\nSUCCESS! Scored {len(df)} candidates.")
        print(f"Rankings saved to: {os.path.abspath(OUTPUT_CSV)}")
        print("\n--- TOP 10 CANDIDATES ---")
        # Adjust column width for nice printing
        pd.set_option('display.max_colwidth', 30)
        print(df_sorted[['design_id', 'ipSAE', 'pDockQ', 'type']].head(10).to_string(index=False))
    else:
        print("\nNo valid structure/PAE pairs found.")


if __name__ == "__main__":
    main()