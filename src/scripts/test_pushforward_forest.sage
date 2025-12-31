import sys
import os
import numpy as np
from sage.all import *

sys.path.append(os.getcwd())

from src.posgeom.stringy_form import StringyCanonicalForm
from src.posgeom.forest_polytope import get_forest_exponents
from src.posgeom.canonical_polytope import eval_canonical_form_dual
from src.posgeom.toric import compute_lattice_basis, get_toric_exponents

def check_stringy_vs_geometric():
    print('Checking Stringy Canonical Form vs Geometric Canonical Form...')
    print('\n--- Calibration: Triangle (2-Simplex) ---')
    exponents = [[0,0], [1,0], [0,1]]
    Z_hom = [vector(QQ, [1] + list(v)) for v in exponents]
    W = [1, 2, 3]
    val_geo = eval_canonical_form_dual(W, exponents)
    print(f'Geometric Value: {val_geo}')
    coeffs = [np.dot(W, z) for z in Z_hom]
    print(f'Coefficients: {coeffs}')
    scf = StringyCanonicalForm(exponents, coeffs)
    s_val = 1e-4
    s_vec = np.ones(2) * s_val
    try:
        val_stringy = scf.compute_canonical_form(s_vec)
        print(f'Stringy Raw: {val_stringy}')
        ratio = val_stringy / val_geo
        print(f'Ratio (Stringy/Geo): {ratio}')
        s_vec2 = s_vec / 10
        val_stringy2 = scf.compute_canonical_form(s_vec2)
        print(f'Stringy Raw (s/10): {val_stringy2}')
        print(f'Scaling factor: {val_stringy2/val_stringy}')
    except Exception as e:
        print(f'Stringy failed: {e}')
    print('\n--- Test: Forest Polytope N=4 ---')
    roots = [0, 1, 2]
    n = 4
    ambient_vertices, _ = get_forest_exponents(n, roots)
    d, basis, origin = compute_lattice_basis(ambient_vertices, saturate=True)
    print(f'Dimension: {d}')
    vertices_u = get_toric_exponents(ambient_vertices, basis, origin)
    W_forest = [1.0 + i*0.5 for i in range(d+1)]
    try:
        val_geo_forest = eval_canonical_form_dual(W_forest, vertices_u)
        print(f'Geometric Forest: {val_geo_forest}')
    except:
        val_geo_forest = 0
        print('Geometric failed (maybe singular W)')
    if val_geo_forest != 0:
        Z_forest_hom = [np.array([1] + list(u)) for u in vertices_u]
        coeffs_forest = [np.dot(W_forest, z) for z in Z_forest_hom]
        scf_forest = StringyCanonicalForm(vertices_u, coeffs_forest)
        s_vec_f = np.ones(d) * 1e-3
        val_str_forest = scf_forest.compute_canonical_form(s_vec_f)
        print(f'Stringy Forest: {val_str_forest}')
        print(f'Ratio: {val_str_forest / val_geo_forest}')
        val_str_forest2 = scf_forest.compute_canonical_form(s_vec_f / 2)
        scaling = val_str_forest2 / val_str_forest
        print(f'Scaling (s -> s/2): {scaling} (Expected ~ {(1/2)**d}?)')
    print('\nConclusion: If ratio is stable (up to s-scaling), the pushforward works.')
    print('For Route K3-A, we need to fix the normalization prefactor.')

if __name__ == '__main__':
    check_stringy_vs_geometric()

