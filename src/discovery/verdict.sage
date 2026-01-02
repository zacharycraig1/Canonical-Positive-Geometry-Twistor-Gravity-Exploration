#!/usr/bin/env sage
"""
Verdict Generator
=================

Collects evidence from all phases and produces a definitive verdict on whether
gravity can be represented as a positive geometry.
"""

from sage.all import *
import sys
import os
from datetime import datetime

# Add project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.discovery.utils import load_checkpoint, log_progress, save_checkpoint


def collect_evidence():
    """
    Collect all evidence from phase results.
    
    Returns dict with all evidence.
    """
    log_progress("VERDICT", "Collecting evidence from all phases...")
    
    evidence = {
        'timestamp': datetime.now().isoformat(),
        'pushforward': {},
        'bcfw': {},
        'positivity': {},
        'verdict': None,
        'confidence': None
    }
    
    # Try to load from all_phases_results first
    all_results = load_checkpoint('all_phases_results.json')
    
    # Load Phase 1 results
    phase1 = None
    if all_results and 'phase1' in all_results:
        phase1 = all_results['phase1']
    else:
        phase1 = load_checkpoint('phase1_result.json')
    
    if phase1:
        evidence['pushforward'] = {
            'ratio_constant': phase1.get('is_constant', False),
            'normalization': phase1.get('normalization'),
            'ratio_stats': phase1.get('ratio_stats'),
            'verdict': phase1.get('verdict'),
            'success': phase1.get('success', False)
        }
        log_progress("VERDICT", f"Phase 1: {phase1.get('verdict', 'unknown')}")
    else:
        log_progress("VERDICT", "Phase 1 results not found", "WARNING")
    
    # Load Phase 2 results
    phase2 = None
    if all_results and 'phase2' in all_results:
        phase2 = all_results['phase2']
    else:
        phase2 = load_checkpoint('phase2_result.json')
    
    if phase2:
        best_formula = phase2.get('best_formula')
        formula_results = phase2.get('formula_results', {})
        
        evidence['bcfw'] = {
            'best_formula': best_formula,
            'formula_results': formula_results,
            'success': phase2.get('success', False),
            'successful_count': phase2.get('successful_count', 0)
        }
        
        if best_formula and best_formula in formula_results:
            best_data = formula_results[best_formula]
            evidence['bcfw']['best_ratio_mean'] = best_data.get('stats', {}).get('mean')
            evidence['bcfw']['best_exact_matches'] = best_data.get('exact_matches', 0)
            evidence['bcfw']['best_is_constant'] = best_data.get('is_constant', False)
        
        log_progress("VERDICT", f"Phase 2: Best formula = {best_formula}")
    else:
        log_progress("VERDICT", "Phase 2 results not found", "WARNING")
    
    # Load Phase 3 results
    phase3 = None
    if all_results and 'phase3' in all_results:
        phase3 = all_results['phase3']
    else:
        phase3 = load_checkpoint('phase3_result.json')
    
    if phase3:
        evidence['positivity'] = {
            'positive_region_found': phase3.get('positive_region_found', False),
            'brute_force': phase3.get('brute_force', {}),
            'optimization': phase3.get('optimization', {}),
            'success': phase3.get('success', False)
        }
        log_progress("VERDICT", f"Phase 3: Positive region found = {evidence['positivity']['positive_region_found']}")
    else:
        log_progress("VERDICT", "Phase 3 results not found", "WARNING")
    
    return evidence


def determine_verdict(evidence):
    """
    Determine final verdict based on evidence.
    
    Returns: (verdict_string, confidence_level, reasoning)
    """
    pushforward = evidence.get('pushforward', {})
    bcfw = evidence.get('bcfw', {})
    positivity = evidence.get('positivity', {})
    
    # Check if we have sufficient evidence
    has_pushforward = pushforward.get('success', False)
    has_bcfw = bcfw.get('success', False)
    has_positivity = positivity.get('success', False)
    
    if not (has_pushforward or has_bcfw or has_positivity):
        return (
            "INSUFFICIENT_DATA",
            "LOW",
            "Not enough evidence collected from phases"
        )
    
    # Key indicators
    bcfw_matches = bcfw.get('best_exact_matches', 0) > 0
    bcfw_constant = bcfw.get('best_is_constant', False)
    bcfw_ratio = bcfw.get('best_ratio_mean')
    positive_found = positivity.get('positive_region_found', False)
    pushforward_constant = pushforward.get('ratio_constant', False)
    
    # Decision logic
    if bcfw_matches and positive_found:
        return (
            "POSITIVE_GEOMETRY_EXISTS",
            "HIGH",
            "BCFW formula matches Hodges exactly AND positive region found"
        )
    
    if bcfw_constant and abs(bcfw_ratio - 1.0) < 0.1 and positive_found:
        return (
            "POSITIVE_GEOMETRY_EXISTS_WITH_NORMALIZATION",
            "MEDIUM",
            f"BCFW matches Hodges up to constant factor ({bcfw_ratio:.6e}) AND positive region found"
        )
    
    if bcfw_constant and positive_found:
        return (
            "POSITIVE_GEOMETRY_LIKELY_EXISTS",
            "MEDIUM",
            f"BCFW has constant ratio ({bcfw_ratio:.6e}) AND positive region found, but ratio not close to 1"
        )
    
    if positive_found and not bcfw_constant:
        return (
            "POSITIVE_REGION_EXISTS_BUT_FORMULA_UNKNOWN",
            "LOW",
            "Positive region found but no consistent BCFW formula"
        )
    
    if bcfw_matches and not positive_found:
        return (
            "FORMULA_EXISTS_BUT_POSITIVITY_UNCLEAR",
            "MEDIUM",
            "BCFW matches Hodges but positive region not found in search"
        )
    
    if not positive_found and not bcfw_constant:
        return (
            "QUESTION_REMAINS_OPEN",
            "LOW",
            "No positive region found and no consistent BCFW formula"
        )
    
    return (
        "INCONCLUSIVE",
        "LOW",
        "Evidence is mixed or incomplete"
    )


def generate_verdict_report(evidence, verdict, confidence, reasoning):
    """
    Generate final VERDICT.md report.
    """
    os.makedirs('discovery_results', exist_ok=True)
    report_path = os.path.join('discovery_results', 'VERDICT.md')
    
    with open(report_path, 'w') as f:
        f.write("# Gravity Positive Geometry: Definitive Verdict\n\n")
        f.write(f"**Date:** {evidence.get('timestamp', 'unknown')}\n\n")
        f.write("---\n\n")
        
        f.write("## Executive Summary\n\n")
        f.write(f"**VERDICT:** {verdict}\n\n")
        f.write(f"**Confidence:** {confidence}\n\n")
        f.write(f"**Reasoning:** {reasoning}\n\n")
        f.write("---\n\n")
        
        f.write("## Phase 1: Pushforward Diagnostic\n\n")
        pushforward = evidence.get('pushforward', {})
        if pushforward.get('success'):
            f.write(f"- **Status:** Completed\n")
            f.write(f"- **Ratio Constant:** {pushforward.get('ratio_constant', False)}\n")
            if pushforward.get('normalization'):
                f.write(f"- **Normalization Factor:** {pushforward.get('normalization'):.6e}\n")
            f.write(f"- **Verdict:** {pushforward.get('verdict', 'unknown')}\n")
        else:
            f.write(f"- **Status:** Failed or Not Run\n")
        f.write("\n")
        
        f.write("## Phase 2: BCFW Amplituhedron\n\n")
        bcfw = evidence.get('bcfw', {})
        if bcfw.get('success'):
            f.write(f"- **Status:** Completed\n")
            f.write(f"- **Best Formula:** {bcfw.get('best_formula', 'none')}\n")
            f.write(f"- **Exact Matches:** {bcfw.get('best_exact_matches', 0)}\n")
            f.write(f"- **Constant Ratio:** {bcfw.get('best_is_constant', False)}\n")
            if bcfw.get('best_ratio_mean'):
                f.write(f"- **Mean Ratio:** {bcfw.get('best_ratio_mean'):.6e}\n")
            f.write(f"- **Successful Tests:** {bcfw.get('successful_count', 0)}\n")
        else:
            f.write(f"- **Status:** Failed or Not Run\n")
        f.write("\n")
        
        f.write("## Phase 3: Positivity Search\n\n")
        positivity = evidence.get('positivity', {})
        if positivity.get('success'):
            f.write(f"- **Status:** Completed\n")
            f.write(f"- **Positive Region Found:** {positivity.get('positive_region_found', False)}\n")
            brute = positivity.get('brute_force', {})
            f.write(f"- **Brute Force Positive Points:** {brute.get('both_positive_count', 0)}\n")
            opt = positivity.get('optimization', {})
            f.write(f"- **Optimization Found Positive:** {opt.get('positive_found', False)}\n")
        else:
            f.write(f"- **Status:** Failed or Not Run\n")
        f.write("\n")
        
        f.write("---\n\n")
        f.write("## Conclusion\n\n")
        f.write(f"{reasoning}\n\n")
        
        f.write("## Next Steps\n\n")
        if verdict == "POSITIVE_GEOMETRY_EXISTS":
            f.write("1. Publish result - major discovery!\n")
            f.write("2. Extend to higher n (n=7, 8, ...)\n")
            f.write("3. Prove all-n statement\n")
        elif verdict == "POSITIVE_GEOMETRY_EXISTS_WITH_NORMALIZATION":
            f.write("1. Identify the normalization factor mathematically\n")
            f.write("2. Understand its geometric origin\n")
            f.write("3. Extend to higher n\n")
        elif verdict == "QUESTION_REMAINS_OPEN":
            f.write("1. Expand search space (more points, different methods)\n")
            f.write("2. Try alternative approaches (intersection theory, CHY)\n")
            f.write("3. Consider algebraic obstruction proofs\n")
        else:
            f.write("1. Review evidence and refine approach\n")
            f.write("2. Run additional tests\n")
            f.write("3. Consult literature for new methods\n")
    
    log_progress("VERDICT", f"Report written to {report_path}")


def main():
    """Generate final verdict."""
    log_progress("VERDICT", "="*70)
    log_progress("VERDICT", "GENERATING FINAL VERDICT")
    log_progress("VERDICT", "="*70)
    
    # Collect evidence
    evidence = collect_evidence()
    
    # Determine verdict
    verdict, confidence, reasoning = determine_verdict(evidence)
    
    log_progress("VERDICT", f"VERDICT: {verdict}")
    log_progress("VERDICT", f"Confidence: {confidence}")
    log_progress("VERDICT", f"Reasoning: {reasoning}")
    
    # Save evidence
    save_checkpoint(evidence, 'final_evidence.json')
    
    # Generate report
    generate_verdict_report(evidence, verdict, confidence, reasoning)
    
    print("\n" + "="*70)
    print("FINAL VERDICT")
    print("="*70)
    print(f"VERDICT: {verdict}")
    print(f"Confidence: {confidence}")
    print(f"\nReasoning:\n{reasoning}")
    print("="*70)
    print(f"\nFull report: discovery_results/VERDICT.md")
    
    return {
        'evidence': evidence,
        'verdict': verdict,
        'confidence': confidence,
        'reasoning': reasoning
    }


if __name__ == '__main__':
    main()

