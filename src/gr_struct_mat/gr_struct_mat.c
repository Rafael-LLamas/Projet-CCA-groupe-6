#include "../gr_struct_mat.h"
#include "flint/flint.h"

// Initialize the structured matrix and allocate memory for G and H
void gr_struct_mat_init( gr_struct_mat_t mat, slong n, slong rank, structure_type_t mat_struct, disp_type_t disp, gr_ctx_t ctx )
{
    mat->struct_t = mat_struct;
    mat->disp_t = disp;
    gr_mat_init( mat->G, n, rank, ctx );
    gr_mat_init( mat->H, n, rank, ctx );
}

// Free the memory of G and H
void gr_struct_mat_clear( gr_struct_mat_t mat, gr_ctx_t ctx )
{
    gr_mat_clear( mat->G, ctx );
    gr_mat_clear( mat->H, ctx );
}

// Return the number of rows
slong gr_struct_mat_nrows( gr_struct_mat_srcptr mat, gr_ctx_t ctx )
{
    return gr_mat_nrows( mat->G, ctx );
}

// Return the rank of the matrix (from generators)
slong gr_struct_mat_rank( gr_struct_mat_srcptr mat, gr_ctx_t ctx )
{
    return gr_mat_ncols( mat->G, ctx );
}

// Print the structured matrix with additional information
void gr_struct_mat_print( gr_struct_mat_srcptr mat, gr_ctx_t ctx )
{
    flint_printf(
        "Structured Matrix: %s, Displacement: %s\n",
        mat->struct_t == TOEPLITZ ? "Toeplitz" : "Hankel",
        mat->disp_t == DISP_PLUS ? "Phi_Plus" : "Phi_Minus" );
    flint_printf( "--- Generator G ---\n" );
    gr_mat_print( mat->G, ctx );
    flint_printf( "\n--- Generator H ---\n" );
    gr_mat_print( mat->H, ctx );
    flint_printf( "\n" );
}