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
A = random_matrix(QQ, n, n)
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

# displace operator Z for n1
Z_N1 = matrix(QQ, n1, n1)
for i in range(n1 - 1):
    Z_N1[i + 1, i] = 1

# displace operator Z for n2
Z_N2 = matrix(QQ, n2, n2)
for i in range(n2 - 1):
    Z_N2[i + 1, i] = 1

e_0_n1 = matrix(QQ, n1, 1, {(0, 0): 1}) # 0 vector with its first element as 1 for n1
e_last_n1 = matrix(QQ, n1, 1, {(n1-1, 0): 1}) # 0 vector with its last element as 1 for n1

e_0_n2 = matrix(QQ, n2, 1, {(0, 0): 1}) 
e_last_n2 = matrix(QQ, n2, 1, {(n2-1, 0): 1})

###-------- partial dense matrices --------###

# last column of a
v_a = matrix(QQ, n1, 1)
for i in range(n1):
    s = 0
    for x in range(i + 1):
        for k in range(rank):
            s += G[i - x, k] * H[n1 - 1 - x, k]
    v_a[i, 0] = s

# last row of a
r_a = matrix(QQ, n1, 1)
for j in range(n1):
    s = 0
    for x in range(j + 1):
        for k in range(rank):
            s += G[n1 - 1 - x, k] * H[j - x, k]
    r_a[j, 0] = s

# bottom-right scalar value of a
s_a = 0
for x in range(n1):
    for k in range(rank):
        s_a += G[n1 - 1 - x, k] * H[n1 - 1 - x, k]

# last column of c
v_c = matrix(QQ, n2, 1)
for l in range(n2):
    s = 0
    for x in range(n1):
        for k in range(rank):
            s += G[n1 + l - x, k] * H[n1 - 1 - x, k]
    v_c[l, 0] = s

# last row of b
r_b = matrix(QQ, n2, 1)
for m in range(n2):
    s = 0
    for x in range(n1):
        for k in range(rank):
            s += G[n1 - 1 - x, k] * H[n1 + m - x, k]
    r_b[m, 0] = s

## UNPACKING a b c d ##

###-------- generator of a --------###
G_A = G_TOP
H_A = H_TOP
a = G_A * H_A.T
print("block a displacement:");print(a)

###-------- generator of b --------###
G_B = G_TOP.augment(Z_N1 * v_a)
H_B = H_BOTTOM.augment(e_0_n2)

b = G_B * H_B.T
print("block b displacement:");print(b)

###-------- generator of c --------###
G_C = G_BOTTOM.augment(e_0_n2)
H_C = H_TOP.augment(Z_N1 * r_a)

c = G_C * H_C.T
print("block c displacement:"); print(c)

###-------- generator of d --------###
G_D = G_BOTTOM.augment((s_a * e_0_n2) + Z_N2 * v_c).augment(e_0_n2)
H_D = H_BOTTOM.augment(e_0_n2).augment(Z_N2 * r_b)
d = G_D * H_D.T
print("block d displacement:");print(d)


## PACKING a b c d ##

rA, rB, rC, rD = G_A.ncols(), G_B.ncols(), G_C.ncols(), G_D.ncols()

G_result = block_matrix([
    [G_A, G_B, matrix(QQ, n1, rC), matrix(QQ, n1, rD)],
    [matrix(QQ, n2, rA), matrix(QQ, n2, rB), G_C, G_D]
])
H_result = block_matrix([
    [H_A, matrix(QQ, n1, rB), H_C, matrix(QQ, n1, rD)],
    [matrix(QQ, n2, rA), H_B, matrix(QQ, n2, rC), H_D]
])


# reverse b bleed (-Z_N1 * v_a)
G_c1 = block_matrix([[ -Z_N1 * v_a ], [ matrix(QQ, n2, 1) ]])
H_c1 = block_matrix([[ matrix(QQ, n1, 1) ], [ e_0_n2 ]])
# print("G_c1");print(G_c1)
# print("H_c1");print(H_c1)

# reverse c bleed (-Z_N1 * r_a)
G_c2 = block_matrix([[ matrix(QQ, n1, 1) ], [ -e_0_n2 ]])
H_c2 = block_matrix([[ Z_N1 * r_a ], [ matrix(QQ, n2, 1) ]])
# print("G_c2");print(G_c2)
# print("H_c2");print(H_c2)

# reverse d bleed part 1 (-Z_N2 * v_c and -s_a)
G_c3 = block_matrix([[ matrix(QQ, n1, 1) ], [ -(s_a * e_0_n2 + Z_N2 * v_c) ]])
H_c3 = block_matrix([[ matrix(QQ, n1, 1) ], [ e_0_n2 ]])
# print("G_c3");print(G_c3)
# print("H_c3");print(H_c3)

# reverse d bleed part 2 (-Z_N2 * r_b)
G_c4 = block_matrix([[ matrix(QQ, n1, 1) ], [ -e_0_n2 ]])
H_c4 = block_matrix([[ matrix(QQ, n1, 1) ], [ Z_N2 * r_b ]])
# print("G_c4");print(G_c4)
# print("H_c4");print(H_c4)

G_result = G_result.augment(G_c1).augment(G_c2).augment(G_c3).augment(G_c4)
H_result = H_result.augment(H_c1).augment(H_c2).augment(H_c3).augment(H_c4)

# print("G_result dimensions:", G_result.dimensions());print(G_result)
# print("H_result dimensions:", H_result.dimensions());print(H_result)

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
print("error in a block:"); print(RES[:n1, :n1] - A[:n1, :n1])
print("error in b block:"); print(RES[:n1, n1:] - A[:n1, n1:])
print("error in c block:"); print(RES[n1:, :n1] - A[n1:, :n1])
print("error in d block:"); print(RES[n1:, n1:] - A[n1:, n1:])