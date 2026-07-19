From Stdlib Require Import ZArith.
Require Import StructuredMatrix.
From mathcomp Require Import all_ssreflect all_algebra.

Check is_toeplitz.
Check is_hankel.

(* Compression is algorithmically required to generator
representation of structured matrices. As arithmetic operations artifically
increase the rank of the generators G and H, during first conceptualization
of this project, we had to define an independent algebraic problematic where
the how do we find the optimal number of columns for both G and H such that
their products don't change. *)

(* Our implementation is an LU decomposition done to each where the
decomposition operations are inversed to the other. For example, performing
an LU dec. on G^T taking only the rank amount of columns of U as G then we
do H = [H] * [P^T] * [rank number of rows of L]. We do the LU dec. on the
transpose of G as we want to perform the operations on columns as the number
of columns of G actually represent the rank of the generators. Detailed
explanation can be found in docs/archive/report and also project doc. too *)

Section GeneratorCompression.
  Variable n alpha r : nat.

  Hypothesis H_rank : r <= alpha.

  (* We will assume that LU decomposition is correct *)
  Variable LU_decomp : Matrix -> (Matrix * Matrix * Matrix).

  Hypothesis LU_correct : forall (G_T P L U : Matrix),
    LU_decomp G_T = (P, L, U) ->
      mat_mul P G_T = mat_mul L U.


