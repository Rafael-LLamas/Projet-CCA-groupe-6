# Week of 01/02/2026

- implement a random Matrix generator (matrix_aux.c)
    + just normal, iterate each elements and assign a random element in it.

- LU decomposition (matrix_aux.c)
    + gr_mat_lu is not proper, returns a LU that in a row echelon form and not a strict upper and lower form
    + extract the LU with the logic sent from the prof
    ! U is in Row Echelon Form (staircase shape),not necessariliy an upper.
        !! should ensure that pivots dont need to be in (i,i)
    + plug this in to the generation of a displacement matrix -> take matrix A, generate G and H^T

- check memory leaks
    + oh boy

----------------------- if enough time -----------------------

- toeplitz matrix addition
    + apparently not as easy as it seems

- toeplitz matrix multiplication
    + apparently not as easy as it seems