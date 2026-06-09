reset()

# HELPER FUNCTIONS ---------------

def compute_G_H(A):
    n = A.nrows()
    m = A.ncols()
    
    Z_n = matrix(QQ, n, n)
    for i in range(n-1): Z_n[i+1, i] = 1
        
    Z_m = matrix(QQ, m, m)
    for i in range(m-1): Z_m[i+1, i] = 1
        
    D = A - Z_n * A * Z_m.T # phi+(A) = D
    rank = D.rank()
    
    if rank == 0:
        G = matrix(QQ, n, 0)
        H = matrix(QQ, m, 0)
        return G, H, 0 
        
    P, L, U = D.LU()
    G_PERM = P * L  
    
    keep = [i for i in range(n) if not U.row(i).is_zero()][:rank] 
    
    G = matrix(QQ, [G_PERM.column(i) for i in keep]).T # n x rank
    H = matrix(QQ, [U.row(i) for i in keep]).T # m x rank
    
    assert G * H.T == D, "data loss"
    return G, H, rank

def reconstruct_disp_plus(G, H):
    nrows_G = G.nrows()
    nrows_H = H.nrows()
    rank = G.ncols()
    A = matrix(QQ, nrows_G, nrows_H)
    for i in range(nrows_G):
        for j in range(nrows_H):
            s = 0
            for k in range(rank):
                for x in range(min(i, j) + 1):
                    s += G[i - x, k] * H[j - x, k]
            A[i, j] = s
    return A

def add_generators(G1, H1, G2, H2):
    return G1.augment(G2), H1.augment(H2)

def neg_generator(G, H):
    return -G, H

def mul_generators(G1, H1, G2, H2):
    n = G1.nrows()
    A = reconstruct_disp_plus(G1, H1)
    B = reconstruct_disp_plus(G2, H2)
    return compute_G_H(A * B)[:2]

def compress_generators(G, H):
    n = G.nrows()
    A = reconstruct_disp_plus(G, H)
    return compute_G_H(A)[:2]


# QUADRANTS (SPLIT & PACK) ---------------

def split_quadrants(G, H):
    n = G.nrows()
    rank = G.ncols()
    n1 = ceil(n/2)
    n2 = floor(n/2)
    
    G_TOP = G[:n1, :] 
    G_BOTTOM = G[n1:, :] 
    H_TOP = H[:n1, :] 
    H_BOTTOM = H[n1:, :] 

    Z_N1 = matrix(QQ, n1, n1); 
    for i in range(n1 - 1): Z_N1[i + 1, i] = 1
        
    Z_N2 = matrix(QQ, n2, n2); 
    for i in range(n2 - 1): Z_N2[i + 1, i] = 1

    e_0_n2 = matrix(QQ, n2, 1, {(0, 0): 1}) 

    v_a = matrix(QQ, n1, 1)
    for i in range(n1):
        v_a[i, 0] = sum(G[i - x, k] * H[n1 - 1 - x, k] for x in range(i + 1) for k in range(rank))

    r_a = matrix(QQ, n1, 1)
    for j in range(n1):
        r_a[j, 0] = sum(G[n1 - 1 - x, k] * H[j - x, k] for x in range(j + 1) for k in range(rank))

    s_a = sum(G[n1 - 1 - x, k] * H[n1 - 1 - x, k] for x in range(n1) for k in range(rank))

    v_c = matrix(QQ, n2, 1)
    for l in range(n2):
        v_c[l, 0] = sum(G[n1 + l - x, k] * H[n1 - 1 - x, k] for x in range(n1) for k in range(rank))

    r_b = matrix(QQ, n2, 1)
    for m in range(n2):
        r_b[m, 0] = sum(G[n1 - 1 - x, k] * H[n1 + m - x, k] for x in range(n1) for k in range(rank))

    G_A = G_TOP
    H_A = H_TOP

    G_B = G_TOP.augment(Z_N1 * v_a)
    H_B = H_BOTTOM.augment(e_0_n2)

    G_C = G_BOTTOM.augment(e_0_n2)
    H_C = H_TOP.augment(Z_N1 * r_a)

    G_D = G_BOTTOM.augment((s_a * e_0_n2) + Z_N2 * v_c).augment(e_0_n2)
    H_D = H_BOTTOM.augment(e_0_n2).augment(Z_N2 * r_b)

    return G_A, H_A, G_B, H_B, G_C, H_C, G_D, H_D

def pack_quadrants(G_A, H_A, G_B, H_B, G_C, H_C, G_D, H_D):
    n1 = G_A.nrows()
    n2 = G_D.nrows()

    Z_N1 = matrix(QQ, n1, n1)
    for i in range(n1 - 1): Z_N1[i + 1, i] = 1
            
    Z_N2 = matrix(QQ, n2, n2)
    for i in range(n2 - 1): Z_N2[i + 1, i] = 1
    
    e_0_n2 = matrix(QQ, n2, 1, {(0, 0): 1}) 
    
    rA, rB, rC, rD = G_A.ncols(), G_B.ncols(), G_C.ncols(), G_D.ncols()

    G_result = block_matrix([
        [G_A, G_B, matrix(QQ, n1, rC), matrix(QQ, n1, rD)],
        [matrix(QQ, n2, rA), matrix(QQ, n2, rB), G_C, G_D]
    ])
    H_result = block_matrix([
        [H_A, matrix(QQ, n1, rB), H_C, matrix(QQ, n1, rD)],
        [matrix(QQ, n2, rA), H_B, matrix(QQ, n2, rC), H_D]
    ])
    
    # v_a: last column of A (col index n1 - 1)
    v_a = matrix(QQ, n1, 1)
    for i in range(n1):
        s = 0
        for k in range(rA):
            for x in range(min(i, n1 - 1) + 1): 
                s += G_A[i - x, k] * H_A[n1 - 1 - x, k]
        v_a[i, 0] = s
    
    # r_a: last row of A (row index n1 - 1)
    r_a = matrix(QQ, n1, 1)
    for j in range(n1):
        s = 0
        for k in range(rA):
            for x in range(min(n1 - 1, j) + 1):
                s += G_A[n1 - 1 - x, k] * H_A[j - x, k]
        r_a[j, 0] = s
    
    # s_a: bottom-right scalar of A (index n1-1, n1-1)
    s_a = 0
    for k in range(rA):
        for x in range(n1): # min(n1-1, n1-1) + 1 = n1
            s_a += G_A[n1 - 1 - x, k] * H_A[n1 - 1 - x, k]

    # v_c: last column of C. (C is n2 x n1, col index n1 - 1)
    v_c = matrix(QQ, n2, 1)
    for l in range(n2):
        s = 0
        for k in range(rC):
            for x in range(min(l, n1 - 1) + 1):
                s += G_C[l - x, k] * H_C[n1 - 1 - x, k]
        v_c[l, 0] = s
    
    # r_b: last row of B. (B is n1 x n2, row index n1 - 1)
    r_b = matrix(QQ, n2, 1)
    for m in range(n2):
        s = 0
        for k in range(rB):
           for x in range(min(n1 - 1, m) + 1):
                s += G_B[n1 - 1 - x, k] * H_B[m - x, k]
        r_b[m, 0] = s
    
    G_c1 = block_matrix([[ -Z_N1 * v_a ], [ matrix(QQ, n2, 1) ]])
    H_c1 = block_matrix([[ matrix(QQ, n1, 1) ], [ e_0_n2 ]])
    
    G_c2 = block_matrix([[ matrix(QQ, n1, 1) ], [ -e_0_n2 ]])
    H_c2 = block_matrix([[ Z_N1 * r_a ], [ matrix(QQ, n2, 1) ]])
    
    G_c3 = block_matrix([[ matrix(QQ, n1, 1) ], [ -(s_a * e_0_n2 + Z_N2 * v_c) ]])
    H_c3 = block_matrix([[ matrix(QQ, n1, 1) ], [ e_0_n2 ]])
    
    G_c4 = block_matrix([[ matrix(QQ, n1, 1) ], [ -e_0_n2 ]])
    H_c4 = block_matrix([[ matrix(QQ, n1, 1) ], [ Z_N2 * r_b ]])
    
    G_final = G_result.augment(G_c1).augment(G_c2).augment(G_c3).augment(G_c4)
    H_final = H_result.augment(H_c1).augment(H_c2).augment(H_c3).augment(H_c4)
    
    return G_final, H_final


# Main ----------------------

def strassen_inverse_generators(G_A, H_A):
    n = G_A.nrows()
    rank = G_A.ncols()

    # BASE CASE
    if n == 1:
        val = sum(G_A[0, i] * H_A[0, i] for i in range(rank))
        if val == 0:
            raise ZeroDivisionError("Matrix is singular")
        G_D = matrix(QQ, 1, 1, [1 / val])
        H_D = matrix(QQ, 1, 1, [1])
        return G_D, H_D

    # SPLIT
    G_a, H_a, G_b, H_b, G_c, H_c, G_d, H_d = split_quadrants(G_A, H_A)

    # RECURSION 1
    G_e, H_e = strassen_inverse_generators(G_a, H_a)

    # SHUR COMPLEMENT
    G_ce, H_ce = mul_generators(G_c, H_c, G_e, H_e)
    G_ce, H_ce = compress_generators(G_ce, H_ce)

    G_eb, H_eb = mul_generators(G_e, H_e, G_b, H_b)
    G_eb, H_eb = compress_generators(G_eb, H_eb)

    G_ceb, H_ceb = mul_generators(G_c, H_c, G_eb, H_eb)
    G_ceb, H_ceb = neg_generator(G_ceb, H_ceb) 
    G_ceb, H_ceb = compress_generators(G_ceb, H_ceb)

    G_S, H_S = add_generators(G_d, H_d, G_ceb, H_ceb)
    G_S, H_S = compress_generators(G_S, H_S)

    # RECURSION 2
    G_t, H_t = strassen_inverse_generators(G_S, H_S)

    # FINAL ASSEMBLY
    G_ebt, H_ebt = mul_generators(G_eb, H_eb, G_t, H_t)
    G_ebt, H_ebt = compress_generators(G_ebt, H_ebt)

    G_tce, H_tce = mul_generators(G_t, H_t, G_ce, H_ce)
    G_tce, H_tce = compress_generators(G_tce, H_tce)

    G_ebtce, H_ebtce = mul_generators(G_ebt, H_ebt, G_ce, H_ce)
    G_ebtce, H_ebtce = compress_generators(G_ebtce, H_ebtce)

    G_x, H_x = add_generators(G_e, H_e, G_ebtce, H_ebtce)
    G_x, H_x = compress_generators(G_x, H_x)

    G_y, H_y = neg_generator(G_ebt, H_ebt)
    G_z, H_z = neg_generator(G_tce, H_tce)

    # PACK
    G_Inv, H_Inv = pack_quadrants(G_x, H_x, G_y, H_y, G_z, H_z, G_t, H_t)
    
    # Final Compression to keep rank minimized
    return compress_generators(G_Inv, H_Inv)


n = 5
success = False
    
while not success:
    A = random_matrix(QQ, n, n)

    if not A.is_invertible():
        continue
            
    G_A, H_A, rank = compute_G_H(A)

    try:
        G_inv, H_inv = strassen_inverse_generators(G_A, H_A)
        success = True 
    except ZeroDivisionError:
        pass
    
print(f"Original Matrix A ({n}x{n})")
print(A)
    
A_inv_reconstructed = reconstruct_disp_plus(G_inv, H_inv)
    
print(f"\nInverse Matrix Reconstructed ({n}x{n})")
print(A_inv_reconstructed)
    
Identity = A * A_inv_reconstructed
print("\nVerification (A * A_inv)")
print(Identity)
    
assert Identity == identity_matrix(QQ, n), "Strassen failed!"
print("\nSuccess")
