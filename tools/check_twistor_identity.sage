from sage.all import *
load("src/hodges.sage")
# from src.hodges import MomentumTwistor, mandelstam_invariant

def verify_twistor_identity():
    print("Verifying s_{012} vs <5012> identity...")
    
    # Retry loop
    for attempt in range(10):
        # 1. Sample random twistor
        twistor = MomentumTwistor(n=6, seed=42 + attempt)
        
        # 2. Compute s_012
        # s_012 = s_01 + s_12 + s_02
        s01 = mandelstam_invariant(twistor, 0, 1)
        s12 = mandelstam_invariant(twistor, 1, 2)
        s02 = mandelstam_invariant(twistor, 0, 2)
        
        if s01 is None or s12 is None or s02 is None:
            continue

        s012 = s01 + s12 + s02
        
        # 3. Compute <5012>
        det5012 = twistor.get_four_bracket(5, 0, 1, 2)
        
        print(f"Attempt {attempt}:")
        print(f"s_012: {s012}")
        print(f"<5012>: {det5012}")
        
        if det5012 != 0:
            ratio = s012 / det5012
            print(f"Ratio s/<...>: {ratio}")
        else:
            print("Det is 0")
            
        # Check if s_012 is zero when <5012>=0
        # ... (rest of logic)
        break
    else:
        print("Could not find valid sample in 10 attempts")
        return

    print("\nConstructing Z with <5012> = 0...")
    # Make Z5 linear combo of Z0, Z1, Z2
    Z = [v for v in twistor.Z]
    Z[5] = Z[0] + Z[1] + Z[2] # dependent
    
    twistor_deg = MomentumTwistor(n=6, Z=Z, check_domain=False)
    det_deg = twistor_deg.get_four_bracket(5, 0, 1, 2)
    print(f"New <5012>: {det_deg}")
    
    s01_d = mandelstam_invariant(twistor_deg, 0, 1)
    s12_d = mandelstam_invariant(twistor_deg, 1, 2)
    s02_d = mandelstam_invariant(twistor_deg, 0, 2)
    
    if s01_d is None or s12_d is None or s02_d is None:
        print("Domain violation in s_ij calculation (likely denominators)")
    else:
        s012_d = s01_d + s12_d + s02_d
        print(f"New s_012: {s012_d}")
        
    # Check if s_012 is zero when <5012>=0
    
if __name__ == "__main__":
    verify_twistor_identity()

