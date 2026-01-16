// pour les arbres
#include "flint/arb.h"
// pour les polys classique
#include "flint/fmpz_poly.h"
// pour les commandes flint
#include "flint/flint.h"
// ATTENTION la on rentre dans le modulo n
#include "flint/fmpz_mod.h"
#include "flint/fmpz_mod_poly.h"

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
  fmpz_mod_ctx_init_ui(n, 12);
  fmpz_mod_poly_init(poly1, n);
  fmpz_mod_poly_init(poly2, n);
  fmpz_mod_poly_init(poly3, n);
  // Pas de si donc faut ui
  fmpz_mod_poly_set_coeff_ui(poly1, 2, 8, n);
  fmpz_mod_poly_set_coeff_ui(poly2, 0, 2, n);

  fmpz_mod_poly_set_coeff_ui(poly2, 3, 12, n);
  fmpz_mod_poly_set_coeff_ui(poly1, 5, 2, n);
  flint_printf("poly1 :\n");
  flint_printf("x = %{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly1, n);
  fmpz_mod_poly_init(poly3, n);
  // ne s'appelle plus classical. Pourquoi ??? JE NE SAIS PAS MOI il faut
  // demander au prof
  fmpz_mod_poly_mul(poly3, poly1, poly2, n);
  fmpz_mod_poly_mul(poly1, poly1, poly2, n);
  flint_printf("poly1 fois poly2:\n");
  flint_printf("x = %{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly1, n);
  flint_printf("poly2 :\n");
  flint_printf("x = %{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly2, n);

  flint_printf("poly3 :\n");
  flint_printf("x = %{fmpz_mod_poly} (%{fmpz_mod_ctx})\n", poly3, n);
  flint_printf("equal = 1, 0 sinon :%d\n",
               fmpz_mod_poly_equal(poly1, poly3, n));

  fmpz_mod_poly_zero(poly1, n);
  flint_printf("1 si vraie, 0 sinon :%d\n", fmpz_mod_poly_is_zero(poly1, n));

  fmpz_mod_poly_clear(poly1, n);
  fmpz_mod_poly_clear(poly2, n);
  fmpz_mod_poly_clear(poly3, n);
  fmpz_mod_ctx_clear(n);
}

int main() {
  // arb == arbre ????
  arb_t x;
  // ok donc il faut init l'arbre
  arb_init(x);
  arb_const_pi(x, 50 * 3.33);
  // utile pour print
  arb_printn(x, 50, 0);
  flint_printf("\n");
  flint_printf("Computed with FLINT-%s\n", flint_version);
  // comme mpfr dans la logique
  arb_clear(x);
  naive_Poly_test();
  naive_mod_poly_test();
}