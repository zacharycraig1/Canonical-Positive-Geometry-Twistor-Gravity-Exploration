import sys
import os
# sys.path.append(os.getcwd())

from sage.all import *

# Check if we are being loaded or run directly
# If loaded, we assume functions are available or we load dependencies?
# If we use load(), we should be careful about double loading.
# But for now, let's explicit load to be safe if run standalone.
# If imported, load() might fail if paths are wrong.
# We assume this script is run from root or loaded from root.

try:
    # Try to access symbols to see if already loaded
    sample_spinor_helicity_conserving
except NameError:
    # Load dependencies
    # We use relative paths assuming execution from repo root
    # Use explicit checking to avoid re-loading if possible, 
    # but load() usually handles re-execution (overwrites).
    if os.path.exists("src/spinor_sampling.sage"):
        load("src/spinor_sampling.sage")
        load("src/hodges.sage")
        load("src/klt.sage")
        load("src/kinematics_map.sage")

def oracle_M6(lambdas, tilde_lambdas, convention="hodges_reduced"):
    """
    Oracle for 6-point MHV gravity amplitude.
    
    Args:
        lambdas: List of spinor lambda vectors
        tilde_lambdas: List of spinor tilde_lambda vectors
        convention: "hodges_reduced" or "klt"
        
    Returns:
        Dict with keys:
            - value: The computed amplitude value
            - M_hodges_reduced: Raw Hodges reduced value
            - M_klt: Raw KLT value
            - normalization: C(lambda) linking them
            - checks: Dict with diff and ratio
            - meta: skip reasons or warnings
    """
    
    adapter = SpinorHelicityAdapter(lambdas, tilde_lambdas)
    
    # 1. Compute Hodges Reduced
    m_hodges, h_reason = hodges_6pt_mhv_reduced(adapter)
    
    # 2. Compute KLT
    m_klt, k_reason = gravity_6pt_mhv_klt(adapter, mandelstam_invariant)
    
    if m_hodges is None or m_klt is None:
        return {
            "value": None,
            "error": f"Computation failed. Hodges: {h_reason}, KLT: {k_reason}"
        }
        
    # 3. Normalization Factor C = -<0 1>^8
    ang01 = adapter.get_angle(0, 1)
    if ang01 == 0:
         return {
            "value": None,
            "error": "Normalization factor <0 1> is zero"
        }
    
    c_factor = - (ang01 ** 8)
    
    # 4. Consistency Check
    expected_klt = c_factor * m_hodges
    diff = m_klt - expected_klt
    ratio = m_klt / expected_klt if expected_klt != 0 else 0
    
    is_consistent = (diff == 0)
    
    result = {
        "M_hodges_reduced": m_hodges,
        "M_klt": m_klt,
        "normalization": c_factor,
        "checks": {
            "diff": diff,
            "ratio": ratio,
            "consistent": is_consistent
        },
        "meta": {
            "hodges_reason": h_reason,
            "klt_reason": k_reason
        }
    }
    
    if convention == "hodges_reduced":
        result["value"] = m_hodges
    elif convention == "klt":
        result["value"] = m_klt
    else:
        result["value"] = None
        result["error"] = f"Unknown convention: {convention}"
        
    return result

if __name__ == "__main__":
    # Self-test
    print("Running Oracle Self-Test...")
    res = sample_spinor_helicity_conserving(n=6, seed=42)
    if res:
        l, lt = res
        out = oracle_M6(l, lt)
        print("Result:", out)
        if out["checks"]["consistent"]:
            print("Oracle Consistent.")
        else:
            print("Oracle Inconsistent!")

