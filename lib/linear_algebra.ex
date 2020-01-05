defmodule MatrexNumerix.LinearAlgebra do
  @moduledoc """
  Linear algebra functions used for vector operations,
  matrix factorization and matrix transformation.
  """

  import Matrex
  alias MatrexNumerix.Math

  @doc """
  The L1 norm of a vector, also known as Manhattan norm.
  """
  @spec l1_norm(Common.maybe_vector()) :: Common.maybe_float()
  def l1_norm(vector) do
    norm(1, vector)
  end

  @doc """
  The L2 norm of a vector, also known as Euclidean norm.
  """
  @spec l2_norm(Common.maybe_vector()) :: Common.maybe_float()
  def l2_norm(vector) do
    norm(2, vector)
  end

  @doc """
  The p-norm of a vector.
  """
  @spec norm(integer, Common.maybe_vector()) :: Common.maybe_float()
  def norm(_p, nil), do: nil

  def norm(p, x = %Matrex{}) do
    s = sum(pow(Matrex.apply(x, :abs), p))
    Math.nth_root(s, p)
  end

end
