// pour les arbres
#include "flint/arb.h"
// pour les polys classique
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
// pour les commandes flint
#include "flint/flint.h"
// ATTENTION la on rentre dans le modulo n
#include "flint/fmpz_mod.h"
#include "flint/fmpz_mod_poly.h"
// Bon je ne vais pas le faire a chaque fois
#include "flint/gr.h"
#include "flint/gr_poly.h"

#include "flint/gr_mat.h"

#include "flint/fmpz_mod_mat.h"

/*
OK donc ici c'est mon terrain de jeu ou je vais ecrire des chose qui peuvent ne
pas fonctionner ou fonctionner partiellement. Je test ici FLINT principalement
pour me faire la main dessus et peux etre trouvé un bug. tout ajout est accepté
tant que rien est supprimer, if est possible de trouvé plusieurs facon de faire
la meme chose, si l'une est plus opti merci de ne pas commenter mais de juste
indiquer la version plus opti Je sais que tu ne vas pas lire ce pavé Arda donc
si tu le lis GG je te dois un chocolat chaud du CROUS.Bien sur il faudra avoir
lu jusqu'a la fin car je risque de te poser des qestions.
*/
void naive_Poly_test() {

  // nouveau type
  fmpz_poly_t poly1, poly2, poly3;
  // ici cela permet de faire l'initialisation normal
  fmpz_poly_init(poly1);
  /*
  ici j'alloue un polynome qui aura max coef 8, il est possible de realloc avec
  fmpz_poly_realloc, il existe une meilleur fonction qui est
  fmpz_poly_fit_length qui va comparer le nombre de coef alloc et si la taille
  est plus grande alors il realloc
  */
  fmpz_poly_init2(poly2, 8);

  // ici un int signé au coeff 2
  fmpz_poly_set_coeff_si(poly1, 2, -8);
  fmpz_poly_set_coeff_si(poly2, 0, -2);
  // ici un int non signé au coeff 0
  fmpz_poly_set_coeff_ui(poly2, 3, 12);
  fmpz_poly_set_coeff_ui(poly1, 5, 2);
  flint_printf("poly1 :\n");
  fmpz_poly_print(poly1);
  flint_printf("\n");
  fmpz_poly_init(poly3);
  // permet de faire une multiplication naive
  fmpz_poly_mul_classical(poly3, poly1, poly2);
  fmpz_poly_mul_classical(poly1, poly1, poly2);
  flint_printf("poly1 fois poly2:\n");
  fmpz_poly_print(poly1);
  flint_printf("\n");
  flint_printf("poly2 :\n");
  fmpz_poly_print(poly2);
  flint_printf("\n");

  flint_printf("poly3 :\n");
  fmpz_poly_print(poly3);
  flint_printf("\n");
  flint_printf("equal = 1, 0 sinon :%d\n", fmpz_poly_equal(poly1, poly3));
  // mets le poly a 0
  fmpz_poly_zero(poly1);
  flint_printf("1 si vraie, 0 sinon :%d\n", fmpz_poly_is_zero(poly1));
  flint_printf("poly2 fois 2^15 :\n");
  // il existe des variants pour tout type de multiplication ou addition avec
  // une constante
  fmpz_poly_scalar_mul_2exp(poly1, poly2, 15);
  fmpz_poly_print(poly1);
  flint_printf("\n");
  fmpz_poly_clear(poly1);
  fmpz_poly_clear(poly2);
  fmpz_poly_clear(poly3);
}

void naive_mod_poly_test() {
  /* C'est concretement les memes fonctions que la version générique de poly
   * sauf que la il faut donc le modulo n et que chaque nom de fonction a mod
   * qui apparait en plus !!! ATTENTION !!! il y a vachement moins de fonctions
   * que pour le poly général donc moins de chose faisable*/
  // tes poly particuliers
  fmpz_mod_poly_t poly1, poly2, poly3;
  // ton n pour le modulo
  fmpz_mod_ctx_t n;
  // toujours init
  fmpz_mod_ctx_init_ui(n, 100);
  fmpz_mod_poly_init(poly1, n);
  fmpz_mod_poly_init(poly2, n);
  fmpz_mod_poly_init(poly3, n);
  // Pas de si donc faut ui
  fmpz_mod_poly_set_coeff_ui(poly1, 2, 8, n);
  fmpz_mod_poly_set_coeff_ui(poly2, 0, 2, n);

  fmpz_mod_poly_set_coeff_ui(poly2, 3, 5, n);
  fmpz_mod_poly_set_coeff_ui(poly1, 5, 2, n);
  flint_printf("poly1 :\n");
  flint_printf("%{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly1, n);
  fmpz_mod_poly_init(poly3, n);
  // ne s'appelle plus classical. Pourquoi ??? JE NE SAIS PAS MOI il faut
  // demander au prof
  fmpz_mod_poly_mul(poly3, poly1, poly2, n);
  fmpz_mod_poly_mul(poly1, poly1, poly2, n);
  flint_printf("poly1 fois poly2:\n");
  flint_printf("%{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly1, n);
  flint_printf("poly2 :\n");
  flint_printf("%{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly2, n);

  flint_printf("poly3 :\n");
  flint_printf("%{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly3, n);
  flint_printf("equal = 1, 0 sinon :%d\n", fmpz_mod_poly_equal(poly1, poly3, n));

  fmpz_mod_poly_zero(poly1, n);
  flint_printf("1 si vraie, 0 sinon :%d\n", fmpz_mod_poly_is_zero(poly1, n));
  // LE SWAP utile ??
  fmpz_mod_poly_swap(poly1, poly2, n);
  flint_printf("poly1 = %{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly1, n);
  flint_printf("poly2 = %{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly2, n);
  fmpz_t x0;
  // Pour obtenir un coeff pas besoin d'init la variable x0
  fmpz_mod_poly_get_coeff_fmpz(x0, poly1, 0, n);
  flint_printf("x0 = %{fmpz}\n", x0);
  fmpz_mod_poly_clear(poly1, n);
  fmpz_mod_poly_clear(poly2, n);
  fmpz_mod_poly_clear(poly3, n);
  fmpz_mod_ctx_clear(n);
  fmpz_clear(x0);
}
void arbarb() {
  /*C'etait le code démo donc je l'ai utiliser direct pour test*/
  // arb == arbre ????
  arb_t x;
  // ok donc il faut init l'arbre
  arb_init(x);
  arb_const_pi(x, 50 * 3.33);
  // utile pour print
  arb_printn(x, 50, 0);
  flint_printf("\n");
  // comme mpfr dans la logique
  arb_clear(x);
}

void naive_gr_poly_test() {
  /*Ce que j'ai compris c'est pour de gros polynome qui ne varie pas et qui sont
   * dance un GROS cercle*/
  gr_poly_t poly1, poly2, poly3;
  gr_ctx_t n;
  // toujours init
  gr_ctx_init_fmpz(n);
  gr_poly_init(poly1, n);
  gr_poly_init(poly2, n);
  gr_poly_init(poly3, n);
  // LETS GO DES SI
  gr_poly_set_coeff_si(poly1, 2, -8, n);
  gr_poly_set_coeff_si(poly2, 0, -2, n);

  gr_poly_set_coeff_ui(poly2, 3, 5, n);
  gr_poly_set_coeff_ui(poly1, 5, 2, n);
  flint_printf("poly1 :\n");
  gr_poly_print(poly1, n);
  flint_printf("\n");
  gr_poly_init(poly3, n);
  // il existe une version classical mais ce n'est pas la bonne snif
  gr_poly_mul(poly3, poly1, poly2, n);
  gr_poly_mul(poly1, poly1, poly2, n);
  flint_printf("poly1 fois poly2:\n");
  gr_poly_print(poly1, n);
  flint_printf("\n");
  flint_printf("poly2 :\n");
  gr_poly_print(poly2, n);
  flint_printf("\n");

  flint_printf("poly3 :\n");
  gr_poly_print(poly3, n);
  flint_printf("\n");
  flint_printf("equal : ");
  truth_println(gr_poly_equal(poly1, poly3, n));

  gr_poly_zero(poly1, n);
  flint_printf("vaut zero : ");
  truth_println(gr_poly_is_zero(poly1, n));
  // LE SWAP utile ??
  gr_poly_swap(poly1, poly2, n);
  gr_poly_print(poly1, n);
  flint_printf("\n");
  gr_poly_print(poly2, n);
  flint_printf("\n");
  gr_ptr x0;
  GR_TMP_INIT(x0, n);
  // ICI IL FAUT INIT ATTENTION
  gr_poly_get_coeff_scalar(x0, poly1, 0, n);
  flint_printf("x0 = ");
  gr_println(x0, n);
  gr_poly_clear(poly1, n);
  gr_poly_clear(poly2, n);
  gr_poly_clear(poly3, n);
  GR_TMP_CLEAR(x0, n);
  gr_ctx_clear(n);
}

void gr_mat_test() {
  gr_mat_t mat1, mat2, mat3;
  gr_ctx_t ctx;
  gr_ptr abubu; // labubu

  gr_ctx_init_fmpz(ctx);
  GR_TMP_INIT(abubu, ctx);
  gr_mat_init(mat1, 10, 10, ctx);
  gr_mat_init_set(mat3, mat1, ctx);
  gr_mat_init(mat2, 10, 10, ctx);

  abubu = gr_mat_entry_ptr(mat1, 3, 8, ctx);
  gr_set_si(abubu, -12, ctx);
  abubu = gr_mat_entry_ptr(mat1, 9, 6, ctx);
  gr_set_ui(abubu, 67, ctx);
  abubu = gr_mat_entry_ptr(mat1, 4, 0, ctx);
  gr_set_si(abubu, -8, ctx);

  abubu = gr_mat_entry_ptr(mat2, 0, 0, ctx);
  gr_set_ui(abubu, 1, ctx);
  abubu = gr_mat_entry_ptr(mat2, 4, 1, ctx);
  gr_set_si(abubu, -299, ctx);
  abubu = gr_mat_entry_ptr(mat2, 2, 1, ctx);
  gr_set_ui(abubu, 893, ctx);

  flint_printf("mat1 = ");
  gr_mat_print(mat1, ctx);
  flint_printf("\n");
  flint_printf("mat2 = ");
  gr_mat_print(mat2, ctx);
  flint_printf("\n");
  gr_mat_mul(mat3, mat1, mat2, ctx);
  gr_mat_mul(mat1, mat1, mat2, ctx);
  flint_printf("mat3 = ");
  gr_mat_print(mat3, ctx);
  flint_printf("\n\n\nmat1 = ");
  gr_mat_print(mat1, ctx);
  flint_printf("\nSont egaux = ");
  // pk c'est faux ????? a voir avec les profs
  truth_println(gr_mat_equal(mat1, mat3, ctx));
  gr_mat_zero(mat1, ctx);
  flint_printf("mat1 vaut zero = ");
  truth_println(gr_mat_is_zero(mat1, ctx));
  gr_mat_swap(mat1, mat2, ctx);
  flint_printf("mat1 = ");
  gr_mat_print(mat1, ctx);
  flint_printf("\n");
  flint_printf("mat2 = ");
  gr_mat_print(mat2, ctx);
  flint_printf("\n");

  gr_clear(abubu, ctx);
  gr_mat_clear(mat1, ctx);
  gr_mat_clear(mat2, ctx);
  gr_mat_clear(mat3, ctx);
  gr_ctx_clear(ctx);
}

void fmpz_mod_mat_test() {

  fmpz_mod_mat_t mat1, mat2, mat3;
  fmpz_t zebi;
  fmpz_mod_ctx_t ctx;
  // je peux aussi utiliser fmpz_init_set_ui mais flemme
  fmpz_init(zebi);
  fmpz_set_ui(zebi, 15);

  fmpz_mod_ctx_init(ctx, zebi);
  fmpz_mod_mat_init(mat1, 10, 10, ctx);
  fmpz_mod_mat_init_set(mat3, mat1, ctx);
  fmpz_mod_mat_init(mat2, 10, 10, ctx);
  fmpz_set_ui(zebi, 12);
  fmpz_mod_mat_set_entry(mat1, 0, 1, zebi, ctx);
  fmpz_set_ui(zebi, 67);
  fmpz_mod_mat_set_entry(mat1, 1, 1, zebi, ctx);
  fmpz_set_ui(zebi, 8);
  fmpz_mod_mat_set_entry(mat1, 1, 0, zebi, ctx);
  fmpz_set_ui(zebi, 1);
  fmpz_mod_mat_set_entry(mat2, 0, 0, zebi, ctx);
  fmpz_set_ui(zebi, 299);
  fmpz_mod_mat_set_entry(mat2, 0, 1, zebi, ctx);
  fmpz_set_ui(zebi, 893);
  fmpz_mod_mat_set_entry(mat2, 1, 1, zebi, ctx);

  flint_printf("mat1 = ");
  fmpz_mod_mat_print_pretty(mat1, ctx);
  flint_printf("\nmat2 = ");
  fmpz_mod_mat_print_pretty(mat2, ctx);

  fmpz_mod_mat_mul(mat3, mat1, mat2, ctx);
  fmpz_mod_mat_mul(mat1, mat1, mat2, ctx);

  flint_printf("\nmat3 = ");
  fmpz_mod_mat_print_pretty(mat3, ctx);

  flint_printf("\nequal 1 sinon 0 = %d", fmpz_mod_mat_equal(mat1, mat3, ctx));
  fmpz_mod_mat_zero(mat1, ctx);
  flint_printf("\nmat1 zero 1 sinon 0 = %d", fmpz_mod_mat_is_zero(mat1, ctx));
  fmpz_mod_mat_swap(mat1, mat2, ctx);
  flint_printf("\nmat1 = ");
  fmpz_mod_mat_print_pretty(mat1, ctx);
  flint_printf("\nmat2 = ");
  fmpz_mod_mat_print_pretty(mat2, ctx);
  flint_printf("\n");

  fmpz_mod_mat_clear(mat1, ctx);
  fmpz_mod_mat_clear(mat2, ctx);
  fmpz_mod_mat_clear(mat3, ctx);
  fmpz_mod_ctx_clear(ctx);
  fmpz_clear(zebi);
}

int main() {

  flint_printf("Computed with FLINT-%s\n", flint_version);
  flint_printf("ZEBI UN ARBRE\n");
  arbarb();
  flint_printf("Partie polynome CLASSSSSSSSSSSSSSIQUE\n");
  naive_Poly_test();
  flint_printf("Partie polynome modulo n DONC LE FUNNNNNNNNNNNNNNN\n");
  naive_mod_poly_test();
  flint_printf("Un que ?? Un polynome GRRRRRRRRR\n");
  naive_gr_poly_test();

  flint_printf("ENFIN LES MATRICES et je commence avec les gr\n");
  gr_mat_test();
  flint_printf("Pour finir avec les fmpz_mod_mat\n");
  fmpz_mod_mat_test();
}