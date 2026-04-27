reset()

# Let's start by representing any matrix with their displacement generators:
# square as inversion is always on squares - must test with nxm if for new logic elsewhere
n = 5

# computes G H logic from project, 
# does LU decomposition on displaced matrix and gets the rank nb of cols and rows accordingly
def compute_G_H(A, n):
    Z = matrix(QQ, n, n)#;print(Z)
    for i in range(n-1):
        Z[i+1, i] = 1
        
    D = A - Z * A * Z.T # phi+(A) = D
    print("Phi+(A):");print(D)
    
    P, L, U = D.LU()
    rank = D.rank()
    
    G_PERM = P * L  # P is a list in project
    
    keep = [i for i in range(n) if not U.row(i).is_zero()][:rank] # we keep rank number of non zero rows (sage special)
    
    G = matrix(QQ, [G_PERM.column(i) for i in keep]).T # n x rank
    H = matrix(QQ, [U.row(i) for i in keep]).T # n x rank
    
    assert G * H.T == D, "data loss"
    return G, H, rank

# get any matrix
A = random_matrix(QQ, n, n, num_bound=10, den_bound=10)
print("A:");print(A)

G, H, rank = compute_G_H(A, n)

# print("G:");print(G)
# print("H:");print(H)

# compute the lengths of the halfs (block sizes)
n1 = ceil(n/2)
n2 = floor(n/2)

print("n1 block size= ");print(n1)
print("n2 block size= ");print(n2)

###-------- tools --------###
G_TOP = G[:n1, :]  # first half of G
G_BOTTOM = G[n1 : n1 + n2, :] # second half of G

H_TOP = H[:n1, :] # First half of H
H_BOTTOM = H[n1 : n1 + n2, :] # Second half of H

## UNPACKING a b c d ##

###-------- generator of a --------###
G_A = G_TOP
H_A = H_TOP
a = G_A * H_A.T
print("block a displacement:");print(a)

###-------- generator of b --------###
# Matches the top-right quadrant of Phi+(A) exactly
G_B = G_TOP
H_B = H_BOTTOM

b = G_B * H_B.T
print("block b displacement:");print(b)

###-------- generator of c --------###
# Matches the bottom-left quadrant of Phi+(A) exactly
G_C = G_BOTTOM
H_C = H_TOP

c = G_C * H_C.T
print("block c displacement:"); print(c)

###-------- generator of d --------###
# Matches the bottom-right quadrant of Phi+(A) exactly
G_D = G_BOTTOM
H_D = H_BOTTOM

d = G_D * H_D.T
print("block d displacement:");print(d)

## PACKING a b c d

# creating the matrices
# [G_A     G_B     0       0    ]
# [0       0       G_C     G_D  ]
# and
# [H_A     0    ]
# [0       H_C  ]
# [H_B     0    ]
# [0       H_D  ]
# so the result will be the multiplication of these both matrices will be the displacement matrix which i can the do the reconstruction
# 
rA, rB, rC, rD = G_A.ncols(), G_B.ncols(), G_C.ncols(), G_D.ncols()

G_result = block_matrix([
    [G_A, G_B, matrix(QQ, n1, rC), matrix(QQ, n1, rD)], # n1 rows
    [matrix(QQ, n2, rA), matrix(QQ, n2, rB), G_C, G_D]  # n2 rows
])
H_result = block_matrix([
    [H_A, matrix(QQ, n1, rB), H_C, matrix(QQ, n1, rD)], # n1 rows
    [matrix(QQ, n2, rA), H_B, matrix(QQ, n2, rC), H_D]  # n2 rows
])

assert G_result * H_result.T == G * H.T, "displacement not correct"

# recostruction_A of project
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

RES = reconstruct_disp_plus(G_result,H_result)
print("repacked matrix A:");print(RES)
assert RES == A, "matrix not correct"