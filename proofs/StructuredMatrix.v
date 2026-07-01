From Stdlib Require Import ZArith.


(* Matrix Definitions *)

Definition Matrix := nat -> nat -> Z.

(* A matrix M is Toeplitz if every diagonal is constant: M(i+1, j+1) = M(i, j) *)
Definition is_toeplitz (M : Matrix) (rows cols : nat) : Prop :=
  forall i j : nat, 
    i + 1 < rows -> 
    j + 1 < cols ->
    M (i + 1) (j + 1) = M i j.

(* Is quasi toeplitz only when there exists a sub matrix that is Toeplitz.
  We will assume that the sub matrix must be at least 2x2 for demonstration purposes
  Normally, any matrix is quasi something until a point.*)
Definition is_quasi_toeplitz (M : Matrix) (rows cols : nat) : Prop :=
  exists sub_rows sub_cols r_offset c_offset : nat,
    2 <= sub_rows /\ 
    2 <= sub_cols /\
    r_offset + sub_rows <= rows /\
    c_offset + sub_cols <= cols /\
    (forall i j : nat,
      i + 1 < sub_rows ->
      j + 1 < sub_cols ->
      M (r_offset + i + 1) (c_offset + j + 1) = M (r_offset + i) (c_offset + j)).

(* A matrix M is Hankel if every anti-diagonal is constant: M(i+1, j) = M(i, j+1) *)
Definition is_hankel (M : Matrix) (rows cols : nat) : Prop :=
  forall i j : nat, 
    i + 1 < rows -> 
    j + 1 < cols ->
    M (i + 1) j = M i (j + 1).

(* Is quasi hankel only when there exists a sub 2x2 matrix that is Hankel *)
Definition is_quasi_hankel (M : Matrix) (rows cols : nat) : Prop :=
  exists sub_rows sub_cols r_offset c_offset : nat,
  2 <= sub_rows /\
  2 <= sub_cols /\
  r_offset + sub_rows <= rows /\
  c_offset + sub_cols <= cols /\
  (forall i j : nat,
    i + 1 < sub_rows ->
    j + 1 < sub_cols ->
    M (r_offset + i + 1) (c_offset + j) = M (r_offset + i) (c_offset + j + 1)).

(* Zero Matrix *)
Definition zero_matrix : Matrix := 
  fun i j => 0%Z.



(* Some Verifications *)

Lemma zero_matrix_is_toeplitz : forall rows cols : nat,
  is_toeplitz zero_matrix rows cols.
Proof.
  intros rows cols.
  unfold is_toeplitz.
  intros i j H_row H_col.
  unfold zero_matrix.
  reflexivity.
Qed.

Lemma zero_matrix_is_quasi_toeplitz : exists rows cols : nat,
  is_quasi_toeplitz zero_matrix rows cols.
Proof.
  exists 2. exists 2.
  unfold is_quasi_toeplitz.
  exists 2. exists 2. exists 0. exists 0.
  repeat split.
  constructor.
  constructor.
  constructor.
  constructor.
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

Lemma zero_matrix_is_quasi_hankel : exists rows cols : nat,
  is_quasi_hankel zero_matrix rows cols.
Proof.
  exists 2. exists 2.
  unfold is_quasi_toeplitz.
  exists 2. exists 2. exists 0. exists 0.
  repeat split.
  constructor.
  constructor.
  constructor.
  constructor.
Qed.
