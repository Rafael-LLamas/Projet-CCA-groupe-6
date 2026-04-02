# I want to compress any matrix A and B such that their
# product has the same result as the product of their
# compressed form.
# 
# Compressed matrices A' and B' should have the same number
# of columns as the rank of their product.
# 
# (In this case, A represents the constructor G and B represents
# the constructor H)

def compress_matrices(G, H):
    # main logic
    E, T = (H.T).echelon_form(transformation=True)
    G_COMP = G * T.inverse()
    H_COMP = E.T
    
    # cut down the excess (based on the zeros on last column of H)
    for i in range(H_COMP.ncols()):
        if H_COMP.column(i).is_zero():
            G_COMP = G_COMP.delete_columns([i])
            H_COMP = H_COMP.delete_columns([i])
        
    return G_COMP, H_COMP

G_EX = matrix(ZZ, 2, 3, [2,3,1,5,4,6])
H_EX = matrix(ZZ, 2, 3, [2,3,1,5,4,6])

print("Matrix G_EX"); print(G_EX)
print("Matrix H_EX"); print(H_EX)

check = G_EX * H_EX.T

G_EX_COMP, H_EX_COMP = compress_matrices(G_EX, H_EX)

print("Matrix G_EX_COMP"); print(G_EX_COMP)
print("Matrix H_EX_COMP"); print(H_EX_COMP)
print("Matrix check"); print(check)

verify = True
verify &= (check == G_EX_COMP * H_EX_COMP.T)
verify &= (check.rank() == G_EX_COMP.ncols())
verify &= (check.rank() == H_EX_COMP.ncols())
print(verify)

print("-" * 50)

# The algorithm is simple in sage cuz I can instantly get
# transformation matrix from echelon_form
# I do not have this priviledge on FLINT unfortunately
# 
# I found this elegant solution on a reddit post made 10 years ago
# 
# As gaussian elimination is a row operation, I can simply
# log the operations made on rows by concatenating an
# identity matrix on the matrix I am the elimination on.
# This identity matrix will have every linear operation done
# multiplied by 1, thus I can simply retrieve the matrix of
# transformation then invert it using FLINT, multiply to G.
# 
# p.s. also known as augmented matrix method (functions name is augment in sage)

def compress_matrices_project(G, H):
    # main logic but less readible
    HT = H.T
    I = identity_matrix(QQ, HT.nrows())
    TEMP = HT.augment(I).echelon_form()
    
    E = TEMP.matrix_from_columns(range(HT.ncols()))
    T = TEMP.matrix_from_columns(range(HT.ncols(), TEMP.ncols()))
    
    G_COMP = G * T.inverse()
    H_COMP = E.transpose()
    
    # cut down the excess (based on the zeros on last column of H)
    for i in range(H_COMP.ncols()):
        if H_COMP.column(i).is_zero():
            G_COMP = G_COMP.delete_columns([i])
            H_COMP = H_COMP.delete_columns([i])
    
    return G_COMP, H_COMP
    
G_EX = matrix(QQ, 2, 3, [2,2,1,5,4,3])
H_EX = matrix(QQ, 2, 3, [1,8,1,9,4,0])    
    
print("Matrix G_EX"); print(G_EX)
print("Matrix H_EX"); print(H_EX)   
    
check = G_EX * H_EX.T   
    
G_EX_COMP, H_EX_COMP = compress_matrices_project(G_EX, H_EX)   
    
print("Matrix G_EX_COMP");print(G_EX_COMP)
print("Matrix H_EX_COMP");print(H_EX_COMP)
print("Matrix check");print(check)
    
verify = True
verify &= (check == G_EX_COMP * H_EX_COMP.T)
verify &= (check.rank() == G_EX_COMP.ncols())
verify &= (check.rank() == H_EX_COMP.ncols())
print(verify)