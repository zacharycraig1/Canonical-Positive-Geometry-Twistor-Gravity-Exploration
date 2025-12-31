import os
import sys
from sage.all import *

# Add src to path
sys.path.append(os.getcwd())

from src.chy_oracle.amplitude_spinor import hodges_6pt_mhv_spinor, ang_bracket, sq_bracket
from src.chy_oracle.matrix_tree import hodges_minor_matrix_tree

def get_random_spinor():
    return vector(QQ, [randint(-10,10), randint(-10,10)])

def generate_kinematics_near_pole(i, j, t_val):
    """
    Generates spinors such that <i,j> ~ t_val, others generic.
    """
    lambdas = [get_random_spinor() for _ in range(6)]
    tilde_lambdas = [get_random_spinor() for _ in range(6)]
    
    # Make l_i = l_j + t * l_dir
    lambdas[j] = get_random_spinor()
    # ensure l_j is not 0
    if lambdas[j] == 0: lambdas[j] = vector(QQ, [1, 0])
        
    l_dir = get_random_spinor()
    # ensure <dir, j> != 0 so <i,j> ~ t
    while ang_bracket(l_dir, lambdas[j]) == 0:
        l_dir = get_random_spinor()
        
    lambdas[i] = lambdas[j] + t_val * l_dir
    
    return lambdas, tilde_lambdas

def check_tree_clearing():
    print("Checking pole structure of Matrix-Tree Sum...")
    
    ts = [1/10, 1/100, 1/1000, 1/10000]
    
    # 1. Non-Adjacent Pole <1,3>
    print("\nProbing non-adjacent pole <1,3> -> 0")
    
    # Use fixed base to isolate t dependence
    l_base = [get_random_spinor() for _ in range(6)]
    lt_base = [get_random_spinor() for _ in range(6)]
    
    # l1 = l3 + t * dir
    l_base[3] = vector(QQ, [1, 0])
    l_dir = vector(QQ, [0, 1])
    # <dir, l3> = <(0,1), (1,0)> = -1. So <1,3> = -t.
    
    # Ensure no other brackets are 0
    valid_base = False
    while not valid_base:
        l_base = [get_random_spinor() for _ in range(6)]
        l_base[3] = vector(QQ, [1, 0])
        valid_base = True
        for a in range(6):
            for b in range(a+1, 6):
                if a==1 and b==3: continue
                if a==1: # check <3, b> roughly
                     if ang_bracket(l_base[3], l_base[b]) == 0: valid_base=False
                elif b==1:
                     if ang_bracket(l_base[a], l_base[3]) == 0: valid_base=False
                else:
                     if ang_bracket(l_base[a], l_base[b]) == 0: valid_base=False
    
    vals = []
    for t in ts:
        L = list(l_base)
        L[1] = L[3] + QQ(t) * l_dir
        try:
            val = hodges_minor_matrix_tree(L, lt_base)
            vals.append(val)
        except ValueError:
            vals.append(None)
            
    print(f"Values: {vals}")
    if vals[-1] is not None and vals[-2] is not None:
        ratio = vals[-1] / vals[-2]
        print(f"Ratio (t=1e-4 / t=1e-3) = {float(ratio)}")
        # If ratio ~ 1, Order 0.
        # If ratio ~ 10, Order -1.
        if abs(abs(ratio) - 1) < 0.1:
            print("Order 0 (Finite) - CORRECT for non-adjacent")
        elif abs(abs(ratio) - 10) < 2:
            print("Order -1 (Simple Pole) - INCORRECT for non-adjacent")
        else:
            print(f"Unknown order (Ratio {ratio})")
            
    # 2. Adjacent Pole <1,2>
    print("\nProbing adjacent pole <1,2> -> 0")
    
    # Reset base
    l_base = [get_random_spinor() for _ in range(6)]
    l_base[2] = vector(QQ, [1, 0])
    l_dir = vector(QQ, [0, 1])
    
    # Ensure validity
    # ... simplify: just retry if error
    
    vals = []
    for t in ts:
        L = list(l_base)
        L[1] = L[2] + QQ(t) * l_dir
        try:
            val = hodges_minor_matrix_tree(L, lt_base)
            vals.append(val)
        except ValueError:
            vals.append(None)
            
    print(f"Values: {vals}")
    if vals[-1] is not None and vals[-2] is not None:
        ratio = vals[-1] / vals[-2]
        print(f"Ratio (t=1e-4 / t=1e-3) = {float(ratio)}")
        
        # Check D_cyc clearing
        # D_cyc contains <1,2>^2 ~ t^2
        # If TreeSum is Order -2, then Product is Order 0.
        if abs(abs(ratio) - 100) < 20:
             print("Order -2 (Double Pole)")
             print("D_cyc (~ t^2) will clear this!")
        elif abs(abs(ratio) - 10) < 2:
             print("Order -1 (Simple Pole)")
             print("D_cyc (~ t^2) will over-clear this (vanish)")
        elif abs(abs(ratio) - 1) < 0.1:
             print("Order 0 (Finite)")
             
if __name__ == "__main__":
    check_tree_clearing()
