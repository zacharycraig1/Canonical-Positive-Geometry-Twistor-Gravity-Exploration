#!/usr/bin/env sage
from sage.all import *

def inspect_data():
    data = load("klt_search/data_batch_0.sobj")
    print(f"Loaded {len(data)} points")
    
    p = data[0]
    print("Target:", p['target'])
    print("Target type:", type(p['target']))
    print("Monomials (first 10):", p['mon_vals'][:10])
    print("PT Vals:", p['pt_vals'])
    
    # Check for zeros
    import numpy as np
    t_val = float(p['target'])
    print("Float Target:", t_val)
    
    m_arr = np.array(p['mon_vals'], dtype=float)
    pt_arr = np.array(p['pt_vals'], dtype=float)
    
    print("Max Mon:", np.max(np.abs(m_arr)))
    print("Max PT:", np.max(np.abs(pt_arr)))
    
    # Check if target is zero
    if t_val == 0:
        print("WARNING: Target is zero")

if __name__ == "__main__":
    inspect_data()








