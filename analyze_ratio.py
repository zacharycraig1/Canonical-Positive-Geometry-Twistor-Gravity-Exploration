
import json
import numpy as np
from itertools import combinations

def analyze():
    with open("ratio_data.json", "r") as f:
        data = json.load(f)
        
    print(f"Loaded {len(data)} points")
    
    # Target variable
    y = np.array([d['ratio'] for d in data])
    
    # Feature matrix
    cr_keys = sorted(data[0]['cross_ratios'].keys())
    print(f"Available cross ratios: {len(cr_keys)}")
    
    X_raw = []
    for d in data:
        row = [d['cross_ratios'][k] for k in cr_keys]
        X_raw.append(row)
    X_raw = np.array(X_raw)
    
    # Add bias term
    X = np.c_[np.ones(X_raw.shape[0]), X_raw]
    
    # Solve least squares: X * w = y
    # w = (X^T X)^-1 X^T y
    try:
        w, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
    except Exception as e:
        print(f"Error in lstsq: {e}")
        return

    # Calculate R^2
    y_pred = X @ w
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot
    
    print(f"Linear Fit R^2: {r2:.6f}")
    print(f"Intercept: {w[0]:.6f}")
    
    # Coefficients
    significant_coeffs = []
    for i, k in enumerate(cr_keys):
        c = w[i+1] # +1 because of bias
        if abs(c) > 1e-4:
            significant_coeffs.append((k, c))
            
    print("\nSignificant Coefficients:")
    significant_coeffs.sort(key=lambda x: abs(x[1]), reverse=True)
    for k, c in significant_coeffs:
        print(f"  {k}: {c:.4f}")
        
    # Check max residual
    max_res = np.max(np.abs(y - y_pred))
    print(f"\nMax Residual: {max_res:.2e}")
    
    if max_res < 1e-5:
        print("SUCCESS: Found a fit!")
    else:
        print("FAILURE: No simple linear fit found.")
        
        # Try Quadratic?
        print("\nTrying Quadratic terms...")
        n_features = X_raw.shape[1]
        X_poly = [np.ones(X_raw.shape[0])]
        X_poly.append(X_raw.T) # Linear terms
        
        # Quadratic terms (cross products)
        # Only take a subset if too many
        # For 15 features, 15*16/2 = 120 + 15 = 135 features.
        # We have 50 points. Too many features.
        # We need more data or fewer features.
        print("Not enough data for full quadratic fit (need > 135 points).")

if __name__ == "__main__":
    analyze()
