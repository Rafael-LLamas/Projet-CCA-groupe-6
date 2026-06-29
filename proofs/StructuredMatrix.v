Require Import ZArith.

Definition Matrix := nat -> nat -> Z.

(* A matrix M is Toeplitz if every diagonal is constant: M(i+1, j+1) = M(i, j) *)
Definition is_toeplitz (M : Matrix) (rows cols : nat) : Prop :=
  forall i j : nat, 
    i + 1 < rows -> 
    j + 1 < cols ->
    M (i + 1) (j + 1) = M i j.

(* A matrix M is Hankel if every anti-diagonal is constant: M(i+1, j) = M(i, j+1) *)
Definition is_hankel (M : Matrix) (rows cols : nat) : Prop :=
  forall i j : nat, 
    i + 1 < rows -> 
    j + 1 < cols ->
    M (i + 1) j = M i (j + 1).
    
(* Zero Matrix *)
Definition zero_matrix : Matrix := 
  fun i j => 0%Z.

Lemma zero_matrix_is_toeplitz : forall rows cols : nat,
  is_toeplitz zero_matrix rows cols.
Proof.
  intros rows cols.
  unfold is_toeplitz.
  intros i j H_row H_col.
  unfold zero_matrix.
  reflexivity.
Qed.

Lemma zero_matrix_is_hankel : forall rows cols : nat, 
is_hankel zero_matrix rows cols.
Proof.
  intros rows cols.
  unfold is_hankel.
  intros i j H_row H_col.
  unfold zero_matrix.
  reflexivity.
Qed.
