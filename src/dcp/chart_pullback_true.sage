
# src/dcp/chart_pullback_true.sage
from sage.all import *

def compute_dcp_residue(N, Vinv, triples, M6, C, Smask):
    """
    Compute the true DCP pullback residue for chart N and boundary Smask.
    
    Args:
        N: list of flats (masks) defining the chart
        Vinv: list of invariant vectors (defining the form Omega)
        triples: list of triples defining the OS form basis
        M6: the matroid (Sage Matroid object)
        C: list of channels (tuples/canonical form) corresponding to M6 groundset
        Smask: mask of the boundary divisor S
        
    Returns:
        Dictionary with status, residue form, solution, etc.
    """
    import time
    t0 = time.time()
    
    # 1. Setup Nested Set Structure
    N_sorted = sorted([int(x) for x in N], key=lambda x: (x.bit_count(), x))
    
    # Identify children for each flat
    # A flat G is a child of F if G < F and there is no H such that G < H < F.
    children = {F: [] for F in N_sorted}
    for i, F in enumerate(N_sorted):
        for j in range(i-1, -1, -1):
            G = N_sorted[j]
            if (G & F) == G: # G subset F
                # Check if G is already a child of some child of F?
                # Actually, simpler: just use immediate subset logic
                # But transitivity...
                # Better: direct check.
                is_immediate = True
                for other in children[F]:
                    if (G & other) == G:
                        is_immediate = False
                        break
                if is_immediate:
                    children[F].append(G)
                    
    # 2. Build Adapted Basis B
    # r(F) is rank of flat F.
    # We need rank function from M6.
    # M6 groundset indices are 0..len(C)-1.
    
    def get_indices(mask):
        out, mm = [], int(mask)
        while mm:
            lsb = mm & -mm
            out.append(lsb.bit_length() - 1)
            mm -= lsb
        return out
        
    B_indices = [] # List of channel indices in the basis
    B_map = {} # Map basis index -> F (the flat it was chosen for)
    
    # Iterate F in N
    # Also need to cover the full space? N should be maximal, so top element is full space?
    # User says "N not maximal" loop in sampler...
    # If N is maximal, it should cover full rank.
    # The matroid rank is 10.
    
    current_rank_sum = 0
    
    # Ensure N is sorted by size/inclusion
    for F in N_sorted:
        F_idxs = get_indices(F)
        r_F = M6.rank(F_idxs)
        
        # Calculate sum of ranks of children
        r_children = 0
        child_basis_union = []
        for child in children[F]:
            # This logic is tricky. r(union children) != sum r(children) generally.
            # But for nested set in wonderful model, we adapt to the chain.
            # The "adapted basis" usually means: extend basis of children to basis of F.
            pass
            
    # Better approach for Adapted Basis:
    # Maintain a current basis set `basis_so_far`.
    # For F in N (small to large):
    #   Extend `basis_so_far` using elements from F to span F.
    #   The new elements are B_F.
    
    basis_so_far = []
    B_structure = [] # List of (F, [indices])
    
    for F in N_sorted:
        F_idxs = get_indices(F)
        # We need to pick elements from F_idxs that are independent of basis_so_far
        # and extend the span to F.
        
        # Current closure
        current_closure = M6.closure(basis_so_far)
        
        # We want to pick `e` in F such that `e` not in closure(basis_so_far)
        # Greedily.
        
        added = []
        # Try elements in F
        # Optimization: only check elements not in current closure?
        # M6.is_independent(basis_so_far + [e])
        
        # Matroid rank is efficient.
        needed = M6.rank(F_idxs) - M6.rank(basis_so_far)
        if needed > 0:
            for e in F_idxs:
                if len(added) == needed: break
                if e in basis_so_far: continue
                
                # Check independence
                if M6.rank(basis_so_far + added + [e]) > M6.rank(basis_so_far + added):
                    added.append(e)
                    
            if len(added) < needed:
                print(f"Warning: Could not span flat {bin(F)}? Needed {needed}, got {len(added)}")
                
        basis_so_far.extend(added)
        B_structure.append((F, added))
        
    # Check if we span the full space (rank 10)
    if M6.rank(basis_so_far) < 10:
        # Extend to full basis if N is not maximal (e.g. seed only)
        # Pick from groundset
        needed = 10 - M6.rank(basis_so_far)
        added = []
        all_idxs = list(range(len(C)))
        for e in all_idxs:
            if len(added) == needed: break
            if M6.rank(basis_so_far + added + [e]) > M6.rank(basis_so_far + added):
                added.append(e)
        basis_so_far.extend(added)
        B_structure.append((0, added)) # 0 means "remainder"
        
    # 3. Define Chart Variables
    # u_F for F in N (blowups)
    # v_b for b in basis (slopes)
    
    # Use global SR if available, else try import
    try:
        from sage.all import SR
    except ImportError:
        pass # Hope it's there
        
    # Create variables
    u_vars = {}
    for F in N_sorted:
        u_name = f"u_{int(F)}"
        u_vars[F] = SR.var(u_name)
        
    v_vars = {}
    v_fixed = {}

    
    # Collect all basis elements
    all_basis_elements = []
    for F, elems in B_structure:
        all_basis_elements.extend(elems)
        
    # Assign v vars
    used_pivots = set()
    for F, elems in B_structure:
        if F == 0: # Remainder, no u_F associated (or u_root?)
            # If F=0 (full space), we don't blow it up (unless N contains full space?)
            # If N contains full space, we treat it like any other flat.
            # If F=0 is just "everything else", we keep all v's.
            for e in elems:
                v_name = f"v_{e}"
                v_vars[e] = SR.var(v_name)
        else:
            # F is in N. We need to fix one v.
            if not elems:
                print(f"Warning: Flat {bin(F)} has empty adapted basis B_F.")
                continue
                
            pivot = elems[0]
            v_fixed[pivot] = SR(1)
            used_pivots.add(pivot)
            
            for e in elems[1:]:
                v_name = f"v_{e}"
                v_vars[e] = SR.var(v_name)
                
    # 4. Define t_b map
    # t_b = (prod_{F in Chain(b)} u_F) * v_b
    t_map = {}
    for b in all_basis_elements:
        # Find Chain(b) = {F in N : b in F}
        # Be careful: b is an index of a channel.
        # Is b "in" F?
        # Yes, F is a set of channels (mask).
        # Check if b-th bit is set in F.
        
        chain = []
        for F in N_sorted:
            if (int(F) >> b) & 1:
                chain.append(F)
                
        val = 1
        for F in chain:
            val *= u_vars[F]
            
        if b in v_fixed:
            val *= v_fixed[b]
        else:
            val *= v_vars[b]
            
        t_map[b] = val
        
    # 5. Compute l_e(p) for all channels
    # l_e(p) = sum_{b in B} c_{e,b} t_b
    # We need to express e in terms of Basis.
    # M6 allows solving. M6 basis?
    # M6 is represented by Vmat (10 x |C|).
    # Columns of Vmat are the vectors v_e in dual space?
    # Wait, Vmat cols are coordinates of channels in some basis.
    # Let Vmat = [v_0 ... v_{m-1}].
    # Then l_e(p) = < v_e, p >.
    # We defined coordinates t_b corresponding to basis elements B.
    # This means we chose a basis for V (or V*) such that...
    # If B is subset of channels, they form a basis of V*.
    # Then p is defined by coordinates t_b = l_b(p).
    # For any other e, l_e(p) is linear comb of l_b(p).
    # l_e = sum alpha_b l_b.
    # We find alpha_b by solving v_e = sum alpha_b v_b.
    
    # We need matrix of basis vectors.
    # Basis vectors are columns of Vmat indexed by B.
    # Vmat is 10xM.
    # Matrix of basis B: 10x10.
    
    # Get Vmat from M6 (it's not directly exposed but we passed M6... wait, M6 is Matroid object)
    # Matroid object doesn't easily give the matrix back if constructed from it?
    # Actually it does: M6.representation()
    
    V_matrix = M6.representation()
    
    # Basis submatrix
    B_indices_list = [b for b in all_basis_elements]
    # Ensure order matches t_map iteration or whatever
    # all_basis_elements is ordered.
    
    Bas_mat = V_matrix[:, B_indices_list]
    
    # Precompute coefficients for all e
    # l_e = sum c_{e,b} t_b
    # v_e = sum c_{e,b} v_b
    # Bas_mat * c = v_e  ->  c = Bas_mat^-1 * v_e
    
    try:
        Bas_inv = Bas_mat.inverse()
    except ValueError:
        return {'status': 'ERROR', 'msg': 'Basis not invertible'}
        
    # Compute l_e for all e
    l_exprs = []
    for e in range(len(C)):
        v_e = V_matrix[:, e]
        coeffs = Bas_inv * v_e
        
        val = 0
        for i, b in enumerate(B_indices_list):
            c = coeffs[i]
            if c != 0:
                val += c * t_map[b]
        l_exprs.append(val)
        
    # 6. Compute dlog forms
    # We only need dlog l_e for e involved in the invariant.
    # And we specifically want the residue at u_S = 0.
    # Smask is one of the F in N (assuming N extends S or contains S).
    # If S not in N, we can't take residue easily (it's not a chart boundary).
    
    if Smask not in u_vars:
        return {'status': 'ERROR', 'msg': 'Boundary S not in chart N'}
        
    u_S_var = u_vars[Smask]
    
    # Residue Operator:
    # 1. Check if l_e contains u_S.
    #    l_e = u_min(e) * unit_e
    #    If u_min(e) contains u_S (factor), then dlog l_e = dlog u_S + ...
    #    If not, dlog l_e is regular at u_S=0.
    #
    # We need to factor out u variables.
    # For each e, find min_u common to all terms.
    # t_b = (prod u) * v.
    # Terms in l_e are c * (prod u) * v.
    # Common factor is intersection of sets of u's.
    # Since u's correspond to nested set, the intersection is just the u's for the smallest flat containing ALL basis elements involved?
    # No, simple GCD of monomials.
    
    # Symbolic GCD might be slow.
    # Better: compute exponent vector for each u.
    # For each term in l_e (which is a t_b), we know the u-exponents.
    # Take min over all terms with non-zero coeff.
    
    l_factorized = []
    
    for e in range(len(C)):
        # Identify active basis elements
        v_e = V_matrix[:, e]
        coeffs = Bas_inv * v_e
        active_b_indices = [i for i, c in enumerate(coeffs) if c != 0]
        
        if not active_b_indices:
            l_factorized.append((SR(0), SR(0))) # Should not happen for valid arrangement
            continue
            
        # Find common u factors
        # For each u_F, check if it is in Chain(b) for ALL active b.
        common_u = []
        for F in N_sorted:
            in_all = True
            for i in active_b_indices:
                b = B_indices_list[i]
                # Is F in Chain(b)? i.e. b in F?
                if not ((int(F) >> b) & 1):
                    in_all = False
                    break
            if in_all:
                common_u.append(u_vars[F])
                
        # Construct factored form
        # l_e = (prod common_u) * (residual)
        # residual = sum c * t_b / (prod common_u)
        
        u_factor = prod(common_u)
        
        residual = 0
        for i in active_b_indices:
            b = B_indices_list[i]
            c = coeffs[i]
            # t_b / u_factor
            # We know t_b has all common_u factors.
            # Just divide them out symbolically or reconstructed.
            # term = c * (prod u in Chain(b) \ common) * v_b
            
            term_u = []
            for F in N_sorted:
                if ((int(F) >> b) & 1):
                     if u_vars[F] not in common_u:
                         term_u.append(u_vars[F])
            
            val = c * prod(term_u)
            if b in v_fixed:
                val *= v_fixed[b]
            else:
                val *= v_vars[b]
            residual += val
            
        l_factorized.append((u_factor, residual))
        
    # 7. Compute Residue Form
    # We substitute u_S = 0 into the RESIDUAL parts?
    # dlog l_e = dlog(u_factor) + dlog(residual).
    # dlog(u_factor) = sum dlog u_F.
    # If u_S is in u_factor, we have a dlog u_S term.
    # If not, no dlog u_S term.
    # dlog(residual) involves d(residual)/residual.
    # Evaluate at u_S = 0.
    # The residue of form Omega at u_S=0 is coeff of dlog u_S.
    # Omega = sum alpha_I dlog l_I.
    # We linearize: dlog l_e = C_e * dlog u_S + (terms without dlog u_S).
    # C_e = 1 if u_S in common_u, else 0.
    # AND terms from dlog(residual)?
    # residual depends on u_S.
    # residual = A + u_S * B ...
    # at u_S=0, residual = A.
    # d(residual) = dA + du_S * B + ...
    # dlog(residual) = (dA + B du_S) / A = dA/A + (B/A) du_S.
    # So we get dlog u_S term from u_factor (coeff 1) AND from residual (coeff B/A * u_S... wait).
    # dlog u_S = du_S / u_S.
    # The term (B/A) du_S is regular (no 1/u_S pole).
    # So the ONLY pole at u_S=0 comes from u_factor!
    # PROVIDED residual is not zero at u_S=0.
    # Since we factored out maximal u_S, residual should not be divisible by u_S.
    # So residual(u_S=0) != 0.
    # Thus dlog(residual) is regular.
    
    # So: Res_{u_S=0} dlog l_e = 1 if u_S in common_u else 0.
    # Wait, this means Res dlog l_e is a CONSTANT (1 or 0)?
    # No, that's just the coefficient of the pole.
    # But we have a wedge product.
    # Res(w1 ^ w2 ^ ...) = (Res w1) ^ w2 ... + ...
    # If multiple terms have dlog u_S, we use linearity.
    # Replace dlog l_e with: (has_pole ? 1 : 0) * dlog u_S + (regular part).
    # Regular part is dlog(residual)|_{u_S=0} + sum_{F!=S} (has_F ? dlog u_F : 0).
    
    # Let D_e = 1 if u_S in common_u else 0.
    # Let Reg_e = sum_{F!=S, u_F in common} dlog u_F + d(residual|_{u_S=0}) / residual|_{u_S=0}.
    # Then dlog l_e ~ D_e dlog u_S + Reg_e.
    
    # Compute Reg_e for all e.
    # Evaluated at u_S = 0.
    
    reg_forms = []
    
    for e in range(len(C)):
        u_fac, res = l_factorized[e]
        
        # Check pole
        has_pole = (u_S_var in u_fac.operands()) if u_fac.operator() == mul else (u_fac == u_S_var)
        # Handle 1 (empty product)
        if u_fac == 1: has_pole = False
        
        # Evaluate residual at u_S = 0
        res_0 = res.subs({u_S_var: 0})
        
        # If res_0 is 0, then we didn't factor out enough u_S!
        # This shouldn't happen with correct logic.
        if res_0 == 0:
            print(f"Error: Residual for channel {e} vanishes at u_S=0!")
            return {'status': 'ERROR', 'msg': f"Vanishing residual for channel {e}"}
            
        # Compute d(res_0)
        # res_0 depends on other u_F and v_b.
        # d(res_0) = sum (diff(res_0, var) * d_var)
        # We need to represent differential forms.
        # Sage doesn't have a simple built-in exterior algebra for symbolic variables mixed with functions?
        # We can just store it as a linear combination of basis differentials.
        # Basis diffs: {dlog u_F (F!=S)} U {d v_b / v_b ? or just d v_b}.
        # Better: use dlog variables for u, and d v variables for v?
        # Or dlog everywhere?
        # For v_b, if they are "slope" variables, maybe d v_b / v_b is natural?
        # Let's use basis of differentials: D_vars = [u_F (F!=S)] + [v_b (b not fixed)].
        # Form: sum C_k dlog(var_k).
        # d(res_0)/res_0 = d log(res_0).
        # d log f = sum (df/dx * x/f) dlog x.
        
        # Identify variables in res_0
        vars_in_res = res_0.variables()
        
        # We want to express dlog(res_0) = sum A_k dlog x_k.
        # A_k = x_k * diff(res_0, x_k) / res_0.
        
        coeffs = {}
        for var in vars_in_res:
            coeff = var * diff(res_0, var) / res_0
            coeffs[var] = coeff
            
        # Add the u_factor part (excluding u_S)
        # u_factor = prod u_F. dlog = sum dlog u_F.
        if u_fac != 1:
            if u_fac.operator() == mul:
                ops = u_fac.operands()
            else:
                ops = [u_fac]
            for op in ops:
                if op != u_S_var:
                    coeffs[op] = coeffs.get(op, 0) + 1
                    
        reg_forms.append((has_pole, coeffs))

    # 8. Compute Residue of Omega
    # Omega = sum alpha_I dlog l_I1 ^ dlog l_I2 ^ ...
    # We iterate over terms in Vinv.
    # For each term (coeff, [e1, e2...]):
    #   Calculate Wedge( (D_e1 dlog uS + Reg_e1), ... )
    #   Extract coeff of dlog uS.
    #   This is: sum_{k} (-1)^k D_ek * (Reg_e1 ^ ... ^ Reg_ek_removed ^ ...)
    
    # We need to sum these up to get the Total Residue Form.
    # The result is a form in the remaining variables (u_F!=S, v_b).
    # Degree = 2 (since original was 3).
    # We want to check if this form is "Good" (Factorizes).
    # Factorization check:
    # Does Res = Omega_L ^ Omega_R?
    # Or strict support check?
    # User says: "Re-test the 'HIT' chart family...".
    # And "Strict L x R only residue/support tests".
    # So we check if the support of the residue form is compatible with L x R.
    # i.e. NO mixed terms.
    
    # We need to collect the form as sum of wedges of basis differentials.
    # Basis diffs: all u_F (F!=S), all v_b.
    # Let's verify dimensions.
    # Omega is 3-form. Residue is 2-form.
    
    # Accumulate the global residue form
    # We use a dictionary to store coefficients of basis wedges.
    # Key: tuple of sorted variable names (or indices).
    
    residue_form = defaultdict(lambda: 0)
    
    # Basis ordering
    all_vars = sorted([str(v) for v in u_vars.values() if v != u_S_var] + 
                      [str(v) for v in v_vars.values()])
    var_to_idx = {v: i for i, v in enumerate(all_vars)}
    
    for v_idx in range(len(Vinv)):
        # For each basis element of OS3 (w0, w1...)
        # Vinv[v_idx] is a list of terms.
        
        # We compute the residue form for w_v_idx.
        # It will be a symbolic expression or dictionary.
        
        form_dict = defaultdict(lambda: 0)
        
        vec_items = Vinv[v_idx] # list of (col_idx, coeff)
        # Note: Vinv passed to scan_chart_proxy was processed into list of (k, val).
        # Here Vinv is the raw list of vectors.
        # We need to handle it.
        
        # Let's assume Vinv is list of dicts or vectors as usual.
        items = vec_items.dict().items() if hasattr(vec_items, 'dict') else vec_items
        
        for k, val in items:
            k = int(k)
            ci, cj, ck = triples[k]
            
            # Form: dlog l_ci ^ dlog l_cj ^ dlog l_ck
            forms = [reg_forms[ci], reg_forms[cj], reg_forms[ck]]
            # Each is (has_pole, coeffs_dict)
            
            # We want coeff of dlog uS.
            # Term has pole if exactly one of them has pole?
            # Or 3? (dlog uS ^ dlog uS = 0).
            # So exactly 1.
            
            poles = [f[0] for f in forms]
            num_poles = sum(poles)
            
            if num_poles != 1:
                continue
                
            # Identify which position has pole
            p_idx = poles.index(True)
            
            # The residue is (-1)^p_idx * (Reg_other1 ^ Reg_other2)
            # Wedge product of two 1-forms.
            # Reg_a = sum C_a_x dlog x
            # Reg_b = sum C_b_y dlog y
            # Reg_a ^ Reg_b = sum_{x<y} (C_a_x C_b_y - C_a_y C_b_x) dlog x ^ dlog y
            
            other_indices = [i for i in range(3) if i != p_idx]
            fa = forms[other_indices[0]][1] # coeffs dict
            fb = forms[other_indices[1]][1]
            
            sign_wedge = 1 if p_idx == 0 else (-1 if p_idx == 1 else 1)
            # wait, Res(dlog uS ^ w) = w.
            # Res(w1 ^ dlog uS ^ w2) = -w1 ^ w2.
            # Res(w1 ^ w2 ^ dlog uS) = w1 ^ w2.
            
            total_coeff = val * sign_wedge
            
            # Compute wedge fa ^ fb
            # Variables involved
            vars_a = set(fa.keys())
            vars_b = set(fb.keys())
            all_v = sorted(list(vars_a | vars_b), key=lambda v: var_to_idx[str(v)])
            
            for i in range(len(all_v)):
                u = all_v[i]
                for j in range(i+1, len(all_v)):
                    v = all_v[j]
                    
                    c_au = fa.get(u, 0)
                    c_av = fa.get(v, 0)
                    c_bu = fb.get(u, 0)
                    c_bv = fb.get(v, 0)
                    
                    # Coeff of dlog u ^ dlog v
                    term = c_au * c_bv - c_av * c_bu
                    
                    if term != 0:
                        pair_key = (str(u), str(v))
                        form_dict[pair_key] += total_coeff * term
                        
        # Store for this basis vector
        residue_form[v_idx] = form_dict

    # 9. Factorization Test
    # We have residue_form[v_idx] = sum C_{uv} dlog u ^ dlog v.
    # We want to find alpha such that sum alpha[v] * residue_form[v] has "good" support.
    # Good support means: for every pair (u, v) with non-zero coeff,
    # u and v must assume L/R values correctly.
    
    # Classify variables u_F and v_b as L, R, or Mixed.
    # u_F: use Lonly/Ronly masks.
    # v_b: depends on b?
    # b is a channel index.
    # If channel b is in L, v_b is L?
    # If channel b is in R, v_b is R?
    # If channel b crosses, v_b is mixed?
    # This seems reasonable.
    
    # We need Lonly/Ronly masks for channels.
    Lonly_int = int(M6.closure(get_indices(Lonly_int))) if isinstance(Lonly_int, int) else 0 # Wait, Lonly passed as mask?
    # Lonly passed to scan_chart_exact2d was mask of Flats.
    # But here we need to know if a variable is L or R.
    # u_F: check F against L/R masks.
    # v_b: check channel b against L/R masks.
    
    # Recalculate Lonly/Ronly channel masks?
    # User passed Lonly, Ronly as *flat* masks?
    # scan_chart_exact2d uses them as flat masks.
    # But user defined L_set = {1,2,3} (particles).
    # build_LR_masks returns masks of CHANNELS.
    # So Lonly is a mask of channels supported on L.
    
    # Classify vars
    var_cls = {}
    
    # u vars
    for F, u_var in u_vars.items():
        if F == int(Smask): continue # Should not appear
        # F is mask of channels.
        # Check if F subset Lonly
        if (int(F) & ~int(Smask) & ~int(Lonly)) == 0: # Ignoring S?
             # If F contains channels not in L, it's not L.
             # Wait, strict check:
             if (int(F) & ~int(Lonly)) == 0:
                 var_cls[str(u_var)] = 'L'
             elif (int(F) & ~int(Ronly)) == 0:
                 var_cls[str(u_var)] = 'R'
             else:
                 var_cls[str(u_var)] = 'M'
                 
    # v vars
    for b, v_var in v_vars.items():
        # b is channel index.
        # Check if bit b is in Lonly
        if (int(Lonly) >> b) & 1:
            var_cls[str(v_var)] = 'L'
        elif (int(Ronly) >> b) & 1:
            var_cls[str(v_var)] = 'R'
        else:
            var_cls[str(v_var)] = 'M'
            
    # Collect constraints
    # Bad terms: (L, L), (R, R), (M, X), (X, M).
    # Good terms: (L, R) or (R, L).
    
    bad_constraints = defaultdict(lambda: [0]*len(Vinv))
    
    for v_idx, form in residue_form.items():
        for pair, coeff in form.items():
            # pair is (str(u), str(v))
            u, v = pair
            cu = var_cls.get(u, 'M')
            cv = var_cls.get(v, 'M')
            
            is_good = (cu == 'L' and cv == 'R') or (cu == 'R' and cv == 'L')
            
            if not is_good:
                # Add to constraints
                # We need sum alpha_v * coeff = 0.
                # Store coeff in list
                bad_constraints[pair][v_idx] += coeff
                
    # Solve constraints
    # Similar to proxy
    bad_list = [(pair, val) for pair, val in bad_constraints.items() if any(v != 0 for v in val)]
    
    result = {
        'status': 'UNKNOWN',
        'bad_count': len(bad_list),
        'constraints': len(bad_list)
    }
    
    if len(Vinv) == 2:
        # Solve
        sol = None
        if not bad_list:
            sol = (1, 0)
            result['status'] = 'HIT'
        else:
            # Pick one
             for pair, coeffs in bad_list:
                c0, c1 = coeffs
                if c0 == 0 and c1 == 0: continue
                if c0 == 0: sol = (1, 0)
                elif c1 == 0: sol = (0, 1)
                else: sol = (c1, -c0)
                break
                
        if sol:
            valid = True
            for pair, coeffs in bad_list:
                if coeffs[0]*sol[0] + coeffs[1]*sol[1] != 0:
                    valid = False
                    break
            if valid:
                result['status'] = 'HIT'
                result['alpha_sol'] = sol
            else:
                result['status'] = 'NO_HIT'
    
    return result


