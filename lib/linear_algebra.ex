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

  def ones_upper(n), do: ones_upper(n, n)
  def ones_upper(ni,nj) do
    m1r = Matrex.ones(1, nj)
    for i <- 1..nj, reduce: Matrex.zeros(ni,nj) do
      m! ->
        m1sl = m1r |> Matrex.submatrix(1..1, i..nj)
        m! |> Matrex.set_submatrix(i..i, i..nj, m1sl)
    end
  end

  def ones_upper_offdiag(n), do: ones_upper_offdiag(n, n)
  def ones_upper_offdiag(ni,nj) do
    m1r = Matrex.ones(1, nj)
    for i <- 1..(nj-1), reduce: Matrex.zeros(ni,nj) do
      m! ->
        m1sl = m1r |> Matrex.submatrix(1..1, (i+1)..nj)
        m! |> Matrex.set_submatrix(i..i, (i+1)..nj, m1sl)
    end
  end

  def reverse(m = %Matrex{}) do
    m |> Enum.reverse() |> Matrex.from_list()
  end
  def reverse!(m = %Matrex{}) do
    m |> Enum.reverse() |> Matrex.from_list() |> Matrex.transpose()
  end

  @doc """
  Solve Uppder Right triangular matrix (symmetric).
  """
  def backward_substitution(uu, y) do
    use Matrex.Operators
    {_, n} = size(uu)
    x = zeros(size(y))

    for i <- n..1, reduce: x do # loop over the rows from bottom to top
      x ->
        r = for k <- (i+1)..n |> Enum.reject(& &1 > n), reduce: y[i] do
              r -> r - uu[i][k] * x[k]
            end
        set(x, i, 1, r / uu[i][i] )
    end
  end

  @doc """
  Solve Uppder Right triangular matrix (symmetric).
  """
  def forward_substitution(ll, b) do
    use Matrex.Operators
    {_ni, n} = size(ll)
    c = zeros(size(b))

    for i <- 1..n, reduce: c do
      c ->
        r = for j <- 1..(i-1) |> Enum.reject(& &1 < 1), reduce: b[i] do
              r ->
                r - ll[i][j] * c[j]
            end

        set(c, i, 1, r / ll[i][i] )
    end
  end

end
