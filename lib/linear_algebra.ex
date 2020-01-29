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
  def l1_norm(vector, weights \\ nil) do
    norm(1, vector, weights)
  end

  @doc """
  The L2 norm of a vector, also known as Euclidean norm.
  """
  @spec l2_norm(Common.maybe_vector()) :: Common.maybe_float()
  def l2_norm(vector, weights \\ nil) do
    norm(2, vector, weights)
  end

  @doc """
  The p-norm of a vector.
  """
  @spec norm(integer, Common.maybe_vector()) :: Common.maybe_float()
  def norm(p, x), do: norm(p, x, nil)

  def norm(_p, nil, _), do: nil

  def norm(p, x = %Matrex{}, nil) do
    s = sum(pow(Matrex.apply(x, :abs), p))
    Math.nth_root(s, p)
  end
  def norm(p, x = %Matrex{}, w = %Matrex{}) do
    s = sum(pow(Matrex.apply(x, :abs), p) |> Matrex.dot(w))
    Math.nth_root(s, p)
  end

end
