#!/usr/bin/env python3
import subprocess
import time
import os

def run_pipeline():
    print("Starting KLT Search Pipeline")
    print("============================")
    
    # Ensure directories
    if not os.path.exists("klt_search"):
        os.makedirs("klt_search")
        
    # Step 1: Setup Basis (Already ran manually, but safe to run again)
    print("\n[1/3] Setting up Basis...")
    subprocess.run(["./sage.ps1", "klt_search/01_setup_basis.sage"], check=True)
    
    # Step 2: Generate Data
    # We need ~8000 unknowns. Let's target 8500 points.
    # We can split this into batches.
    total_points = 8500
    batch_size = 500
    num_batches = total_points // batch_size
    
    print(f"\n[2/3] Generating {total_points} data points in {num_batches} batches...")
    
    # Run batches sequentially (could be parallelized, but let's keep it simple)
    for i in range(num_batches):
        print(f"  Running Batch {i+1}/{num_batches}...")
        subprocess.run(["./sage.ps1", "klt_search/02_generate_data.sage", str(i), str(batch_size)], check=True)
        
    # Step 3: Solve
    print("\n[3/3] Solving System...")
    subprocess.run(["./sage.ps1", "klt_search/03_solve_system.sage"], check=True)
    
    print("\nPipeline Complete.")

if __name__ == "__main__":
    run_pipeline()

