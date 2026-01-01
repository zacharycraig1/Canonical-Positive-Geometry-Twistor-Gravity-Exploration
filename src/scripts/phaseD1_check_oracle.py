import os
import sys
from sage.all import *

# Add src to path
sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor, ang_bracket, sq_bracket

def get_random_spinor():
    return vector(QQ, [randint(-10,10), randint(-10,10)])

def check_oracle_poles():
    print("Checking pole structure of Oracle (hodges_6pt_mhv_spinor)...")
    
    ts = [1/10, 1/100, 1/1000]
    
    # 1. Non-Adjacent Pole <1,3>
    print("\nProbing non-adjacent pole <1,3> -> 0")
    l_base = [get_random_spinor() for _ in range(6)]
    l_base[3] = vector(QQ, [1, 0])
    l_dir = vector(QQ, [0, 1])
    
    # Ensure validity
    # ...
    
    vals = []
    for t in ts:
        L = list(l_base)
        L[1] = L[3] + QQ(t) * l_dir
        res, _ = hodges_6pt_mhv_spinor(L, [get_random_spinor() for _ in range(6)])
        vals.append(res)
        
    print(f"Values: {vals}")
    if vals[-1] is not None and vals[-2] is not None and vals[-1] != 0:
        ratio = vals[-1] / vals[-2]
        print(f"Ratio (t=1e-3 / t=1e-2) = {float(ratio)}")
        if abs(abs(ratio) - 1) < 0.1:
            print("Order 0 (Finite) - CORRECT")
        else:
            print(f"Order != 0. Ratio {ratio}")
    else:
        print("Values are None or 0")

    # 2. Adjacent Pole <1,2>
    print("\nProbing adjacent pole <1,2> -> 0")
    l_base = [get_random_spinor() for _ in range(6)]
    l_base[2] = vector(QQ, [1, 0])
    l_dir = vector(QQ, [0, 1])
    
    vals = []
    lt = [get_random_spinor() for _ in range(6)]
    for t in ts:
        L = list(l_base)
        L[1] = L[2] + QQ(t) * l_dir
        res, _ = hodges_6pt_mhv_spinor(L, lt)
        vals.append(res)
        
    print(f"Values: {vals}")
    if vals[-1] is not None and vals[-2] is not None:
        ratio = vals[-1] / vals[-2]
        print(f"Ratio (t=1e-3 / t=1e-2) = {float(ratio)}")
        # Expecting 1/t^2 => Ratio 100
        if abs(abs(ratio) - 100) < 10:
            print("Order -2 (Double Pole) - CORRECT")
        elif abs(abs(ratio) - 10) < 1:
            print("Order -1 (Simple Pole)")
        else:
            print(f"Unknown Order. Ratio {ratio}")

if __name__ == "__main__":
    check_oracle_poles()






