from sageall import *

def to_matrix_qq(rows):
    if isinstance(rows, Matrix):
        return rows.change_ring(QQ)
    return matrix(QQ, rows)

def stack_rows(*mats):
    mats_q = [to_matrix_qq(m) for m in mats if m is not None]
    if not mats_q:
        return matrix(QQ, 0, 0)
    out = mats_q[0]
    for m in mats_q[1:]:
        if out.ncols() == 0:
            out = m
        elif m.ncols() == 0:
            continue
        else:
            if out.ncols() != m.ncols():
                # Pad smaller to match
                maxc = max(out.ncols(), m.ncols())
                out = block_matrix([[out, zero_matrix(QQ, out.nrows(), maxc - out.ncols())]])
                m   = block_matrix([[m,   zero_matrix(QQ, m.nrows(),   maxc - m.ncols())]])
            out = block_matrix([[out],[m]])
    return out

def augment_cols(*mats):
    mats_q = [to_matrix_qq(m) for m in mats if m is not None]
    if not mats_q:
        return matrix(QQ, 0, 0)
    out = mats_q[0]
    for m in mats_q[1:]:
        if out.nrows() == 0:
            out = m
        elif m.nrows() == 0:
            continue
        else:
            if out.nrows() != m.nrows():
                # Pad smaller to match
                maxr = max(out.nrows(), m.nrows())
                out = block_matrix([[out],[zero_matrix(QQ, maxr - out.nrows(), out.ncols())]])
                m   = block_matrix([[m],[  zero_matrix(QQ, maxr - m.nrows(),   m.ncols())]])
            out = block_matrix([[out, m]])
    return out

def rank_qq(M):
    return to_matrix_qq(M).rank()

def right_nullspace_matrix(M):
    Mq = to_matrix_qq(M)
    K = Mq.right_kernel()
    return K.basis_matrix().transpose()  # columns span nullspace

def left_nullspace_matrix(M):
    Mq = to_matrix_qq(M)
    K = Mq.left_kernel()
    return K.basis_matrix().transpose()












