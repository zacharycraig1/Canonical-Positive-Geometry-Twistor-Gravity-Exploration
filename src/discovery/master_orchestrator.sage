#!/usr/bin/env sage
"""
Master Orchestrator for Gravity Positive Geometry Discovery
============================================================

Runs all three phases in parallel, collects evidence, and produces verdict.

This is the main entry point for the discovery pipeline.
"""

from sage.all import *
import sys
import os
import time
from multiprocessing import Process, Queue
import json
from datetime import datetime

# Add project root
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, project_root)
sys.path.insert(0, os.path.join(project_root, 'src'))

# Import with explicit path handling
import importlib.util

# Load utils
utils_path = os.path.join(project_root, 'src', 'discovery', 'utils.sage')
if os.path.exists(utils_path):
    spec = importlib.util.spec_from_file_location("discovery_utils", utils_path)
    discovery_utils = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(discovery_utils)
    log_progress = discovery_utils.log_progress
    save_checkpoint = discovery_utils.save_checkpoint
    load_checkpoint = discovery_utils.load_checkpoint
else:
    # Fallback: try direct import
    from src.discovery.utils import log_progress, save_checkpoint, load_checkpoint

# Import phases
from src.discovery.phase1_pushforward import run_phase1_parallel
from src.discovery.phase2_bcfw import run_phase2_parallel
from src.discovery.phase3_positivity import run_phase3_parallel


def run_phase1_wrapper(queue):
    """Wrapper to run Phase 1 in separate process."""
    try:
        result = run_phase1_parallel(n_trials=1000, n_cores=4)
        queue.put(('phase1', result))
    except Exception as e:
        queue.put(('phase1', {'success': False, 'error': str(e)}))


def run_phase2_wrapper(queue):
    """Wrapper to run Phase 2 in separate process."""
    try:
        result = run_phase2_parallel(n_trials=10000, n_cores=8)
        queue.put(('phase2', result))
    except Exception as e:
        queue.put(('phase2', {'success': False, 'error': str(e)}))


def run_phase3_wrapper(queue):
    """Wrapper to run Phase 3 in separate process."""
    try:
        result = run_phase3_parallel(n_brute=10000, n_optimization_restarts=100, n_cores=4)
        queue.put(('phase3', result))
    except Exception as e:
        queue.put(('phase3', {'success': False, 'error': str(e)}))


def run_all_phases_parallel():
    """
    Run all three phases in parallel.
    
    Returns dict with results from all phases.
    """
    log_progress("MASTER", "="*70)
    log_progress("MASTER", "GRAVITY POSITIVE GEOMETRY DISCOVERY")
    log_progress("MASTER", "="*70)
    log_progress("MASTER", "Starting all phases in parallel...")
    
    t_start = time.time()
    
    # Create queue for results
    result_queue = Queue()
    
    # Start all phases
    processes = []
    
    p1 = Process(target=run_phase1_wrapper, args=(result_queue,))
    p1.start()
    processes.append(('phase1', p1))
    log_progress("MASTER", "Started Phase 1: Pushforward Diagnostic")
    
    p2 = Process(target=run_phase2_wrapper, args=(result_queue,))
    p2.start()
    processes.append(('phase2', p2))
    log_progress("MASTER", "Started Phase 2: BCFW Amplituhedron")
    
    p3 = Process(target=run_phase3_wrapper, args=(result_queue,))
    p3.start()
    processes.append(('phase3', p3))
    log_progress("MASTER", "Started Phase 3: Positivity Search")
    
    # Collect results
    results = {}
    completed = set()
    
    while len(completed) < 3:
        try:
            phase_name, result = result_queue.get(timeout=3600)  # 1 hour timeout per phase
            results[phase_name] = result
            completed.add(phase_name)
            log_progress("MASTER", f"Phase {phase_name} completed")
        except:
            # Timeout - check if processes are still alive
            for phase_name, proc in processes:
                if not proc.is_alive() and phase_name not in completed:
                    log_progress("MASTER", f"Phase {phase_name} process died", "WARNING")
                    completed.add(phase_name)
                    results[phase_name] = {'success': False, 'error': 'process_died'}
    
    # Wait for all processes
    for phase_name, proc in processes:
        proc.join(timeout=60)
        if proc.is_alive():
            log_progress("MASTER", f"Force terminating {phase_name}", "WARNING")
            proc.terminate()
    
    elapsed = time.time() - t_start
    
    log_progress("MASTER", f"All phases completed in {elapsed:.1f}s ({elapsed/3600:.2f} hours)")
    
    return results


def run_all_phases_sequential():
    """
    Run phases sequentially (fallback if parallel fails).
    """
    log_progress("MASTER", "Running phases sequentially...")
    
    results = {}
    
    # Phase 1
    log_progress("MASTER", "Starting Phase 1...")
    try:
        results['phase1'] = run_phase1_parallel(n_trials=1000, n_cores=4)
    except Exception as e:
        log_progress("MASTER", f"Phase 1 failed: {e}", "ERROR")
        results['phase1'] = {'success': False, 'error': str(e)}
    
    # Phase 2
    log_progress("MASTER", "Starting Phase 2...")
    try:
        results['phase2'] = run_phase2_parallel(n_trials=10000, n_cores=8)
    except Exception as e:
        log_progress("MASTER", f"Phase 2 failed: {e}", "ERROR")
        results['phase2'] = {'success': False, 'error': str(e)}
    
    # Phase 3
    log_progress("MASTER", "Starting Phase 3...")
    try:
        results['phase3'] = run_phase3_parallel(n_brute=10000, n_optimization_restarts=100, n_cores=4)
    except Exception as e:
        log_progress("MASTER", f"Phase 3 failed: {e}", "ERROR")
        results['phase3'] = {'success': False, 'error': str(e)}
    
    return results


def main():
    """Main entry point."""
    log_progress("MASTER", "Master Orchestrator starting...")
    
    # Try parallel first, fall back to sequential
    try:
        results = run_all_phases_parallel()
    except Exception as e:
        log_progress("MASTER", f"Parallel execution failed: {e}, falling back to sequential", "WARNING")
        results = run_all_phases_sequential()
    
    # Save all results
    save_checkpoint(results, 'all_phases_results.json')
    
    log_progress("MASTER", "All phases complete. Results saved to discovery_results/all_phases_results.json")
    
    # Print summary
    print("\n" + "="*70)
    print("DISCOVERY PIPELINE SUMMARY")
    print("="*70)
    
    for phase_name, result in results.items():
        success = result.get('success', False)
        status = "SUCCESS" if success else "FAILED"
        print(f"{phase_name.upper()}: {status}")
        if not success:
            print(f"  Error: {result.get('error', 'unknown')}")
    
    print("="*70)
    print("\nRun src/discovery/verdict.sage to generate final verdict.")
    
    return results


if __name__ == '__main__':
    main()

