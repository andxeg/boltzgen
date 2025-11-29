#!/usr/bin/env python3
"""Test script to score a single candidate with Boltz-2 + IPSAE."""

import subprocess
import tempfile
import shutil
from pathlib import Path
import sys

IPSAE_SCRIPT = "/lambda/nfs/Nipah-hackathon/andrei/IPSAE/ipsae.py"
TEST_CIF = "/lambda/nfs/Nipah-hackathon/andrei/boltzgen/protein_competition/results/20designs_5budget/intermediate_designs_inverse_folded/refold_cif/job_A_medium_00.cif"

print("Testing Boltz-2 + IPSAE scoring pipeline on one candidate...")
print(f"Test file: {Path(TEST_CIF).name}")

# Create temp directory
temp_dir = tempfile.mkdtemp(prefix="boltz_test_")
print(f"Temp directory: {temp_dir}")

try:
    # Step 1: Run Boltz-2 prediction
    print("\n[1/2] Running Boltz-2 prediction...")
    cmd = [
        "boltz", "predict",
        str(TEST_CIF),
        "--out_dir", str(temp_dir),
        "--num_workers", "1",
        "--accelerator", "gpu",
        "--devices", "1",
        "--diffusion_samples", "1"
    ]

    print(f"Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        print(f"❌ Boltz prediction failed!")
        print(f"STDERR: {result.stderr}")
        sys.exit(1)

    print("✓ Boltz prediction completed")

    # Step 2: Find the PAE and CIF files
    print("\n[2/2] Looking for output files...")
    temp_path = Path(temp_dir)
    pred_dirs = list(temp_path.glob("predictions/*"))

    if not pred_dirs:
        print("❌ No prediction directory found!")
        sys.exit(1)

    pred_dir = pred_dirs[0]
    print(f"Prediction directory: {pred_dir}")

    pae_files = list(pred_dir.glob("pae_*_model_0.npz"))
    cif_files = list(pred_dir.glob("*_model_0.cif"))

    if not pae_files or not cif_files:
        print(f"❌ Missing output files!")
        print(f"PAE files found: {pae_files}")
        print(f"CIF files found: {cif_files}")
        sys.exit(1)

    pae_file = pae_files[0]
    cif_file = cif_files[0]

    print(f"✓ Found PAE file: {pae_file.name}")
    print(f"✓ Found CIF file: {cif_file.name}")

    # Step 3: Run IPSAE
    print("\n[3/3] Running IPSAE...")
    cmd = [
        "python", IPSAE_SCRIPT,
        str(pae_file),
        str(cif_file),
        "10",  # PAE cutoff
        "10"   # Distance cutoff
    ]

    print(f"Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"❌ IPSAE failed!")
        print(f"STDERR: {result.stderr}")
        sys.exit(1)

    print("✓ IPSAE completed successfully!")
    print("\n" + "="*80)
    print("IPSAE Output:")
    print("="*80)
    print(result.stdout)
    print("="*80)

finally:
    # Cleanup
    print(f"\nCleaning up temporary directory: {temp_dir}")
    shutil.rmtree(temp_dir, ignore_errors=True)
    print("✓ Done!")
