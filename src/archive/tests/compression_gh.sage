reset()

# n,k are the number of rows and columns of G and H
n,k = 10,4
# rg and rh are the ranks of G and H
rg,rh = 3,2

# Here I sample G and H in ZZ^{n \times k}
G,H = random_matrix(ZZ,n,rg), random_matrix(ZZ,n,rh)
G = G * random_matrix(ZZ,rg,k)
H = H * random_matrix(ZZ,rh,k)

print("G:")
print(G)
print("H:")
print(H)

# Here we compute G1,H1 using the echelon form of G
# G1 = G * P^T
# G H^T = G1 * P^{-T} * H^T
# Thus H1 = (P^{-T} * H^T)^T = H * P^{-1}
G1,P = (G.transpose()).echelon_form(transformation=True)
G1 = G1.transpose()
H1 = Matrix(ZZ,H * P.inverse()) # Here we cast in ZZ for next echelon form

print("\nG1 = G*P^T before removing columns ?",G1 == G*P.transpose())

indices = []
for i in range(G1.ncols()):
    if G1.column(i).is_zero():
        indices.append(i)

# Here we can delete the columns of G1 that are zero
# BUT, we can also delete the corresponding rows in H1^T (in other words the columns of H1)
G1 = G1.delete_columns(indices)
H1 = H1.delete_columns(indices)

print("\nWe check that everything is fine after first reduction:", G * H.transpose() == G1 * H1.transpose())

# Now we do the same with H
H2,P = (H1.transpose()).echelon_form(transformation=True)
H2 = H2.transpose()
G2 = G1 * P.inverse()

indices = []
for i in range(H2.ncols()):
    if H2.column(i).is_zero():
        indices.append(i)

# Here we can delete the columns of H2 that are zero and the corresponding rows in G2^T
G2 = G2.delete_columns(indices)
H2 = H2.delete_columns(indices)

print("\nWe check that everything is fine after second reduction:", G * H.transpose() == G2 * H2.transpose())

print("G2:")
print(G2)
print("H2:")
print(H2)

print("\nRank of GH^T:",(G*H.transpose()).rank())
print("Number of columns of G2 and H2:",G2.ncols(),',', H2.ncols())
