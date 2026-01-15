// pour les arbres
#include "flint/arb.h"
// pour les polys
#include "flint/fmpz_poly.h"
// pour les commandes flint
#include "flint/flint.h"

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
  fmpz_poly_print(poly1);
  flint_printf("\n");
  fmpz_poly_init(poly3);
  // permet de faire une multiplication naive
  fmpz_poly_mul_classical(poly3, poly1, poly2);
  fmpz_poly_mul_classical(poly1, poly1, poly2);

  fmpz_poly_print(poly1);
  flint_printf("\n");
  fmpz_poly_print(poly2);
  flint_printf("\n");

  // mets le poly a 0

  fmpz_poly_print(poly3);
  flint_printf("\n");
  flint_printf("%d\n", fmpz_poly_equal(poly1, poly3));
  fmpz_poly_zero(poly1);
  fmpz_poly_clear(poly1);
  fmpz_poly_clear(poly2);
  fmpz_poly_clear(poly3);
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
}