#!/usr/bin/env sage
from sage.all import *
load('src/hodges.sage')
load('src/spinor_helicity.sage')
load('src/kinematics_map.sage')

def check_extraction():
    tw = MomentumTwistor(n=6, check_domain=False)
    l, tl, x = extract_spinors_from_twistor(tw)
    
    print(f"Lambdas length: {len(l)}")
    print(f"TildeLambdas length: {len(tl)}")
    
    if len(tl) < 6:
        print("ERROR: TildeLambdas too short")
        return

    adapter = SpinorHelicityAdapter(l, tl)
    print(f"Adapter n: {adapter.n}")
    
    # Try accessing square bracket
    try:
        sq = adapter.get_square(0, 5)
        print(f"Square [0,5]: {sq}")
    except Exception as e:
        print(f"Error accessing square: {e}")

if __name__ == "__main__":
    check_extraction()








