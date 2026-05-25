#include "flint/flint.h"
#include "flint/gr_mat.h"
#include "flint/ulong_extras.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "addition.h"
#include "compression.h"
#include "displacement_matrices.h"
#include "inverse_toeplitz.h"
#include "multiplication.h"
#include "random_toeplitz.h"

int gr_mat_ldlt_toeplitz(gr_mat_t G,gr_mat_t L,gr_mat_t d ,gr_ctx_t ctx){
    /* 
    Si j'ai bien compris, on résout Ax=b de façon plus opti que de faire 
    x = A-1 b
    L et d init dans la fonction
    */
   gr_mat_t G_cur;
   int error = GR_SUCCESS;
   int n = gr_mat_nrows(G,ctx);
   gr_mat_init(L,n,n,ctx);
   gr_mat_init(L,n,1,ctx);
   gr_mat_init_set(G_cur,G,ctx);
   for(int k = 0;k <n;k++){
    slong a,b,c,d,e;
    gr_get_ui(&a,gr_mat_entry_srcptr(G_cur,k,0,ctx),ctx);
    gr_get_ui(&b,gr_mat_entry_srcptr(G_cur,k,1,ctx),ctx);
    c = a-b;
    if(c <= 0){
        flint_printf("ERROR C <= 0");
        return GR_TEST_FAIL;
    }
    gr_set_ui(gr_mat_entry_ptr(d,k,0,ctx),c,ctx);
    gr_set_ui(gr_mat_entry_ptr(L,k,k,ctx),1,ctx);
    for(int i = 1 ; i < n-k-1;i++){
    }


   }


}