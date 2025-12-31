#!/usr/bin/env sage
# Iterative search runner with automatic analysis and parameter adjustment
# Run with: sage run_iterative_search.sage

load('54.sage')

import json
import os
from datetime import datetime

def analyze_results(run_dir):
    """Analyze results from a run and suggest next steps."""
    report_path = os.path.join(run_dir, "hit_report_intersection.json")
    if not os.path.exists(report_path):
        return None
    
    with open(report_path, 'r') as f:
        report = json.load(f)
    
    final_dim = report.get('final_dim', None)
    boundaries = report.get('boundaries', [])
    intersection_log = report.get('intersection_log', [])
    
    analysis = {
        'final_dim': final_dim,
        'num_boundaries': len(boundaries),
        'boundaries_processed': len(intersection_log),
        'success': final_dim == 1,
        'suggestions': []
    }
    
    if final_dim is None:
        analysis['suggestions'].append("Run failed - check logs")
    elif final_dim > 1:
        analysis['suggestions'].append(f"Dimension {final_dim} > 1 - need more constraints")
        analysis['suggestions'].append("Try: More boundaries, different invariant mode, or different seed")
    elif final_dim == 1:
        analysis['suggestions'].append("SUCCESS! Found dim=1 candidate")
        analysis['suggestions'].append("Next: Verify Hodge structure")
    elif final_dim == 0:
        analysis['suggestions'].append("Intersection became empty - too many constraints")
        analysis['suggestions'].append("Try: Fewer boundaries or different boundary set")
    
    # Check which boundaries succeeded
    successful_boundaries = [b for b in intersection_log if b.get('status') == 'HIT']
    analysis['successful_boundaries'] = len(successful_boundaries)
    
    if len(successful_boundaries) < len(boundaries):
        analysis['suggestions'].append(f"Only {len(successful_boundaries)}/{len(boundaries)} boundaries succeeded")
        analysis['suggestions'].append("Consider removing problematic boundaries")
    
    return analysis

def suggest_next_config(analysis, iteration):
    """Suggest next configuration based on analysis."""
    if analysis is None:
        return {
            'INVARIANT_MODE': 'S3xS3',
            'MULTI_STRATEGY_SEARCH': True,
            'BASE_SEED': 42 + iteration * 1000,
            'TOTAL_TRIALS': 100000,
        }
    
    suggestions = {}
    
    if analysis['final_dim'] > 1:
        # Try different invariant mode
        if iteration % 3 == 1:
            suggestions['INVARIANT_MODE'] = 'S3xS3Z2'
        elif iteration % 3 == 2:
            suggestions['INVARIANT_MODE'] = 'S6'
        else:
            suggestions['INVARIANT_MODE'] = 'S3xS3'
        
        suggestions['MULTI_STRATEGY_SEARCH'] = True
        suggestions['BASE_SEED'] = 42 + iteration * 1000
        suggestions['TOTAL_TRIALS'] = 100000
        
    elif analysis['final_dim'] == 0:
        # Too many constraints - reduce boundaries
        suggestions['INTERSECT_BOUNDARY_MODE'] = 'CUSTOM'
        suggestions['INTERSECT_BOUNDARIES'] = [(1,2,3), (1,2,4), (1,3,4), (2,3,4)]  # Start with 4
        suggestions['BASE_SEED'] = 42 + iteration * 1000
    
    return suggestions

def main_iterative():
    """Run iterative search with automatic analysis."""
    print("="*70)
    print("ITERATIVE SEARCH: 6-Point MHV Gravity Positive Geometry")
    print("="*70)
    
    max_iterations = 5
    iteration = 1
    
    while iteration <= max_iterations:
        print(f"\n{'='*70}")
        print(f"ITERATION {iteration}/{max_iterations}")
        print(f"{'='*70}")
        
        # Run the main search
        try:
            main()
        except Exception as e:
            print(f"Error in iteration {iteration}: {e}")
            iteration += 1
            continue
        
        # Find the most recent run directory
        runs_dir = "runs"
        if os.path.exists(runs_dir):
            run_dirs = [d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))]
            if run_dirs:
                latest_run = max(run_dirs, key=lambda x: os.path.getmtime(os.path.join(runs_dir, x)))
                run_path = os.path.join(runs_dir, latest_run)
                
                # Analyze results
                analysis = analyze_results(run_path)
                
                if analysis:
                    print(f"\nAnalysis of iteration {iteration}:")
                    print(f"  Final dimension: {analysis['final_dim']}")
                    print(f"  Boundaries processed: {analysis['boundaries_processed']}")
                    print(f"  Success: {analysis['success']}")
                    print(f"\nSuggestions:")
                    for sug in analysis['suggestions']:
                        print(f"  - {sug}")
                    
                    # If we found dim=1, we're done!
                    if analysis['success']:
                        print(f"\n{'!'*70}")
                        print("SUCCESS! Found dim=1 candidate!")
                        print(f"Results in: {run_path}")
                        print(f"{'!'*70}")
                        break
                
                # Suggest next configuration
                next_config = suggest_next_config(analysis, iteration)
                if next_config:
                    print(f"\nNext iteration will try:")
                    for key, value in next_config.items():
                        print(f"  {key} = {value}")
                    
                    # Apply suggestions (would need to modify global config)
                    # For now, just document
                    config_file = os.path.join(run_path, "next_config_suggestions.json")
                    with open(config_file, 'w') as f:
                        json.dump(next_config, f, indent=2)
                    print(f"  Suggestions saved to: {config_file}")
        
        iteration += 1
    
    print(f"\n{'='*70}")
    print("ITERATIVE SEARCH COMPLETE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main_iterative()




