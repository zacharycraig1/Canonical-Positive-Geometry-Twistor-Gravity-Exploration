#!/usr/bin/env python3
"""
Analyze results and suggest optimizations.
Can run without SageMath to analyze existing results.
"""
import json
import os
import glob
from pathlib import Path

def find_latest_run():
    """Find the most recent run directory."""
    runs_dir = Path("runs")
    if not runs_dir.exists():
        return None
    
    run_dirs = [d for d in runs_dir.iterdir() if d.is_dir()]
    if not run_dirs:
        return None
    
    return max(run_dirs, key=lambda x: x.stat().st_mtime)

def analyze_run(run_dir):
    """Analyze a run and return suggestions."""
    report_path = run_dir / "hit_report_intersection.json"
    
    if not report_path.exists():
        return {
            'status': 'no_report',
            'suggestions': ['Run may have failed - check logs']
        }
    
    with open(report_path, 'r') as f:
        report = json.load(f)
    
    final_dim = report.get('final_dim')
    boundaries = report.get('boundaries', [])
    intersection_log = report.get('intersection_log', [])
    invariant_mode = report.get('invariant_mode', 'unknown')
    
    analysis = {
        'final_dim': final_dim,
        'invariant_mode': invariant_mode,
        'num_boundaries': len(boundaries),
        'boundaries_processed': len(intersection_log),
        'success': final_dim == 1,
        'suggestions': []
    }
    
    # Analyze dimension progression
    dims = []
    for entry in intersection_log:
        if 'null_dim' in entry:
            dims.append(entry['null_dim'])
    
    analysis['dimension_progression'] = dims
    
    # Generate suggestions
    if final_dim is None:
        analysis['suggestions'].append("Run incomplete or failed")
    elif final_dim == 1:
        analysis['suggestions'].append("✓ SUCCESS! Found dim=1 candidate")
        analysis['suggestions'].append("Next: Verify Hodge structure")
    elif final_dim > 1:
        analysis['suggestions'].append(f"Dimension {final_dim} > 1 - need more constraints")
        if len(dims) > 0:
            last_dim = dims[-1]
            analysis['suggestions'].append(f"Last boundary reduced to dim={last_dim}")
            if last_dim == final_dim:
                analysis['suggestions'].append("No further reduction - try different invariant mode")
        
        # Suggest next invariant mode
        if invariant_mode == 'S3xS3':
            analysis['suggestions'].append("Try: S3xS3Z2 or S6")
        elif invariant_mode == 'S3xS3Z2':
            analysis['suggestions'].append("Try: S6")
        else:
            analysis['suggestions'].append("Try: Different boundary set or more boundaries")
    elif final_dim == 0:
        analysis['suggestions'].append("Intersection became empty - too many constraints")
        analysis['suggestions'].append("Try: Fewer boundaries or different boundary set")
    
    return analysis

def suggest_next_config(analysis, iteration):
    """Suggest next configuration."""
    config = {
        'iteration': iteration + 1,
        'changes': []
    }
    
    if analysis.get('final_dim') == 1:
        config['action'] = 'SUCCESS - verify candidate'
        return config
    
    if analysis.get('final_dim', 0) > 1:
        # Try different invariant mode
        current_mode = analysis.get('invariant_mode', 'S3xS3')
        if current_mode == 'S3xS3':
            config['INVARIANT_MODE'] = 'S3xS3Z2'
            config['changes'].append('Switch to S3xS3Z2 invariants')
        elif current_mode == 'S3xS3Z2':
            config['INVARIANT_MODE'] = 'S6'
            config['changes'].append('Switch to S6 invariants')
        else:
            # Try more boundaries or different approach
            config['TOTAL_TRIALS'] = 100000  # Increase search
            config['changes'].append('Increase TOTAL_TRIALS for better charts')
    
    if analysis.get('final_dim') == 0:
        # Too many constraints
        config['INTERSECT_BOUNDARY_MODE'] = 'CUSTOM'
        config['INTERSECT_BOUNDARIES'] = [(1,2,3), (1,2,4), (1,3,4), (2,3,4)]
        config['changes'].append('Reduce to 4 boundaries')
    
    config['BASE_SEED'] = 42 + iteration * 1000
    config['changes'].append(f'New seed: {config["BASE_SEED"]}')
    
    return config

def main():
    """Main analysis function."""
    print("="*70)
    print("ANALYZING RESULTS AND SUGGESTING OPTIMIZATIONS")
    print("="*70)
    
    run_dir = find_latest_run()
    if not run_dir:
        print("No runs found. Starting fresh search...")
        return {
            'action': 'start_fresh',
            'config': {
                'INVARIANT_MODE': 'S3xS3',
                'MULTI_STRATEGY_SEARCH': True,
                'BASE_SEED': 42,
                'TOTAL_TRIALS': 50000
            }
        }
    
    print(f"\nAnalyzing: {run_dir}")
    analysis = analyze_run(run_dir)
    
    print(f"\nResults:")
    print(f"  Final dimension: {analysis.get('final_dim', 'unknown')}")
    print(f"  Invariant mode: {analysis.get('invariant_mode', 'unknown')}")
    print(f"  Boundaries processed: {analysis.get('boundaries_processed', 0)}")
    print(f"  Success: {analysis.get('success', False)}")
    
    if analysis.get('dimension_progression'):
        print(f"  Dimension progression: {analysis['dimension_progression']}")
    
    print(f"\nSuggestions:")
    for sug in analysis.get('suggestions', []):
        print(f"  - {sug}")
    
    # Suggest next config
    iteration = int(run_dir.name.split('_')[0].split('-')[-1]) if '_' in run_dir.name else 1
    next_config = suggest_next_config(analysis, iteration)
    
    print(f"\nNext configuration:")
    for key, value in next_config.items():
        if key != 'changes':
            print(f"  {key} = {value}")
    
    # Save suggestions
    suggestions_file = run_dir / "next_config_suggestions.json"
    with open(suggestions_file, 'w') as f:
        json.dump(next_config, f, indent=2)
    print(f"\nSuggestions saved to: {suggestions_file}")
    
    return analysis, next_config

if __name__ == '__main__':
    main()













