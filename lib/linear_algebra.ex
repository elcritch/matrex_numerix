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

  def backward_substitution(u, y) do
    use Matrex.Operators

    {_ni, n} = size(u)

    x = Matrex.zeros(1,n) |> Matrex.set(1, n, y[n] / u[n][n])

    for i <- (n-1)..1, reduce: x do
      x ->
        s = y[i]

        s =
          for j <- n..(i+1), reduce: s do
            s ->
              s - u[i][j] * x[j]
          end

        x |> Matrex.set(1, i, s / u[i][i] )
    end
    |> Matrex.transpose()
  end

  def forward_substitution(l, b) do
    use Matrex.Operators

    {_ni, n} = size(l)

    y = Matrex.zeros(1,n) |> Matrex.set(1, 1, b[1])

    # Y(1) = B(1)
    # DO  I = 2, N-1
    #   S = B(I)
    #   DO  J = 1, I-1
    #     S = S - DPROD(L(I,J), Y(J))
    #   END DO
    # END DO

    for i <- 2..n, reduce: y do
      y ->
        s = b[i]

        s =
          for j <- 1..(i-1), reduce: s do
            s ->
              s - l[i][j] * y[j]
          end

        y |> Matrex.set(1, i, s / l[i][i] )
    end
    |> Matrex.transpose()
  end

end
