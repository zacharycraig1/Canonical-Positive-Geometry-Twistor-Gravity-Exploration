#!/usr/bin/env sage
"""
Shared Utilities for Discovery Pipeline
========================================

Common functions used across all discovery phases.
"""

from sage.all import *
import json
import os
import sys
from datetime import datetime

# Add project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.posgeom.forest_polytope import get_forest_polynomial, get_forest_exponents
from src.posgeom.physics_map import eval_edge_vars_from_spinors
from src.chy_oracle.laplacian_bridge import reconstruct_mhv_from_laplacian
from src.chy_oracle.kinematics_samples import sample_spinors_from_twistor
from src.chy_oracle.amplitude_spinor import ang_bracket


def get_reference_spinors():
    """Standard reference spinors."""
    return vector(QQ, [1, 0]), vector(QQ, [0, 1])


def compute_hodges_oracle(n, seed, roots=(0, 1, 2)):
    """
    Compute Hodges determinant using the trusted laplacian_bridge.
    This is the ORACLE - the ground truth.
    
    Returns: (hodges_value, lambdas, tildes) or (None, None, None) if failed
    """
    try:
        lambdas, tildes = sample_spinors_from_twistor(n=n, seed=seed)
        x, y = get_reference_spinors()
        
        M, status = reconstruct_mhv_from_laplacian(lambdas, tildes, x, y, roots=roots)
        return M, lambdas, tildes
    except Exception as e:
        return None, None, None


def compute_forest_polynomial_value(n, roots, lambdas, tildes, x, y):
    """
    Evaluate F_{n,R}(z) at the given kinematics.
    """
    try:
        F_poly = get_forest_polynomial(n, roots)
        z_map = eval_edge_vars_from_spinors(lambdas, tildes, x, y)
        
        R_ring = F_poly.parent()
        eval_dict = {}
        for var_name in R_ring.variable_names():
            if var_name in z_map:
                eval_dict[R_ring(var_name)] = z_map[var_name]
        
        F_val = F_poly.subs(eval_dict)
        
        if hasattr(F_val, 'is_constant') and F_val.is_constant():
            return F_val.constant_coefficient()
        return F_val
    except Exception:
        return None


def is_singular_kinematics(lambdas, tildes):
    """
    Check if kinematics are singular (degenerate).
    """
    n = len(lambdas)
    for i in range(n):
        ip1 = (i + 1) % n
        ang = ang_bracket(lambdas[i], lambdas[ip1])
        if ang == 0:
            return True
    return False


def save_checkpoint(data, filename):
    """Save checkpoint data to JSON file."""
    os.makedirs('discovery_results', exist_ok=True)
    filepath = os.path.join('discovery_results', filename)
    
    # Convert Sage objects to JSON-serializable
    json_data = {}
    for key, value in data.items():
        if isinstance(value, (list, tuple)):
            json_data[key] = [float(v) if isinstance(v, (int, float, Integer, Rational)) else str(v) for v in value]
        elif isinstance(value, (int, float, Integer, Rational)):
            json_data[key] = float(value)
        elif isinstance(value, bool):
            json_data[key] = value
        elif isinstance(value, str):
            json_data[key] = value
        elif value is None:
            json_data[key] = None
        else:
            json_data[key] = str(value)
    
    with open(filepath, 'w') as f:
        json.dump(json_data, f, indent=2)


def load_checkpoint(filename):
    """Load checkpoint data from JSON file."""
    filepath = os.path.join('discovery_results', filename)
    if not os.path.exists(filepath):
        return None
    
    with open(filepath, 'r') as f:
        return json.load(f)


def log_progress(phase, message, level="INFO"):
    """Log progress message with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_msg = f"[{timestamp}] [{phase}] [{level}] {message}"
    print(log_msg)
    
    # Also write to log file
    os.makedirs('discovery_results', exist_ok=True)
    with open('discovery_results/discovery.log', 'a') as f:
        f.write(log_msg + "\n")


def statistical_analysis(values, name="values"):
    """
    Perform statistical analysis on a list of values.
    
    Returns dict with mean, std, min, max, cv (coefficient of variation)
    """
    if not values:
        return None
    
    values_float = [float(v) for v in values if v is not None]
    
    if not values_float:
        return None
    
    mean_val = sum(values_float) / len(values_float)
    variance = sum((x - mean_val)**2 for x in values_float) / len(values_float)
    std_val = variance ** 0.5
    
    min_val = min(values_float)
    max_val = max(values_float)
    
    cv = std_val / abs(mean_val) if mean_val != 0 else float('inf')
    
    return {
        'mean': mean_val,
        'std': std_val,
        'min': min_val,
        'max': max_val,
        'cv': cv,
        'count': len(values_float)
    }


def is_constant_ratio(ratios, threshold=0.01):
    """
    Determine if a list of ratios is constant (within threshold CV).
    
    Returns: (is_constant, stats_dict)
    """
    stats = statistical_analysis(ratios)
    if stats is None:
        return False, None
    
    is_const = stats['cv'] < threshold
    return is_const, stats


