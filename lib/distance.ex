defmodule MatrexNumerix.Distance do
  @moduledoc """
  Distance functions between two vectors.
  """

  import Matrex.Guards

  import MatrexNumerix.LinearAlgebra

  alias MatrexNumerix.{Common, Correlation, Statistics}

  @doc """
  Mean squared error, the average of the squares of the errors
  betwen two vectors, i.e. the difference between predicted
  and actual values.
  """
  @spec mse(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def mse(x = %Matrex{}, y = %Matrex{}) do
    p = Matrex.pow(x |> Matrex.subtract(y), 2)
    Statistics.mean(p)
  end
  def mse(x = %Matrex{}, y = %Matrex{}, w = %Matrex{}) do
    p = Matrex.pow(x |> Matrex.subtract(y), 2)
    Statistics.mean(p |> Matrex.dot(w))
  end


  @doc """
  Root mean square error of two vectors, or simply the
  square root of mean squared error of the same set of
  values. It is a measure of the differences between
  predicted and actual values.
  """
  @spec rmse(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def rmse(vector1, vector2) do
    :math.sqrt(mse(vector1, vector2))
  end
  def rmse(vector1, vector2, weights) do
    :math.sqrt(mse(vector1, vector2, weights))
  end

  @doc """
  The Pearson's distance between two vectors.
  """
  @spec pearson(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def pearson(vector1, vector2) do
    1.0 - Correlation.pearson(vector1, vector2)
  end

  @doc """
  The Minkowski distance between two vectors.
  """
  @spec minkowski(Matrex.t(), Matrex.t(), integer) :: Common.maybe_float()
  def minkowski(x = %Matrex{}, y = %Matrex{}, p \\ 3) do
    norm(p, x |> Matrex.subtract(y))
  end
  def minkowski(x = %Matrex{}, y = %Matrex{}, w = %Matrex{}, p) do
    norm(p, x |> Matrex.subtract(y), Matrex.dot(w))
  end

  @doc """
  The Euclidean distance between two vectors.
  """
  @spec euclidean(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def euclidean(x = %Matrex{}, y = %Matrex{}) do
    l2_norm(x |> Matrex.subtract(y))
  end
  def euclidean(x = %Matrex{}, y = %Matrex{}, w = %Matrex{}) do
    l2_norm(x |> Matrex.subtract(y), w)
  end

  @doc """
  The Euclidean distance between two vectors.
  """
  @spec sq_euclidean(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def sq_euclidean(x = %Matrex{}, y = %Matrex{}) do
    l2_norm(x |> Matrex.subtract(y)) |> :math.pow(2)
  end
  def sq_euclidean(x = %Matrex{}, y = %Matrex{}, w = %Matrex{}) do
    l2_norm(x |> Matrex.subtract(y), w) |> :math.pow(2)
  end

  @doc """
  The Manhattan distance between two vectors.
  """
  @spec manhattan(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def manhattan(x = %Matrex{}, y = %Matrex{}) do
    l1_norm(x |> Matrex.subtract(y) )
  end
  def manhattan(x = %Matrex{}, y = %Matrex{}, w = %Matrex{}) do
    l1_norm(x |> Matrex.subtract(y), w)
  end

  @doc """
  The Jaccard distance (1 - Jaccard index) between two vectors.
  """
  @spec jaccard(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def jaccard(
        matrex_data(rows1, columns1, _data1, _first),
        matrex_data(rows2, columns2, _data2, _second)
      ) when rows1 != rows2 or columns1 != columns2,
      do: raise %ArgumentError{message: "incorrect sizes"}

  def jaccard(vector1, vector2) do
    vector1
    |> Stream.zip(vector2)
    |> Enum.reduce({0, 0}, fn {x, y}, {intersection, union} ->
      case {x, y} do
        {x, y} when x == 0 or y == 0 ->
          {intersection, union}

        {x, y} when x == y ->
          {intersection + 1, union + 1}

        _ ->
          {intersection, union + 1}
      end
    end)
    |> to_jaccard_distance
  end

  def diff_conv(method, x, y) do
    w = Matrex.ones(x |> Matrex.size() |> elem(0), 1)
    diff_conv(method, x, y, w)
  end

  def diff_conv(method,
                vector_data(columns1, _data1, _first) = xx,
                vector_data(columns2, _data2, _second) = yy,
                weights) do

    {1, nobsx} = Matrex.size(xx)
    {1, nobsy} = Matrex.size(yy)
    weights = weights || Matrex.ones(1, 1)
    {dimw, nobsw} = Matrex.size(weights)

    (nobsx == nobsy == nobsw) || %ArgumentError{message: "nobs(xx) != nobs(yy) != nobs(weights)"}

    dist_func =
      case method do
        :euclidian -> &euclidean/3
        :sq_euclidian -> &sq_euclidean/3
        :manhattan -> &manhattan/3
        :mse -> &mse/3
        :rmse -> &rmse/3
        :minkowski -> &minkowski/3
        custom when is_function(custom) -> custom
      end

    for i <- 1..nobsx, into: [] do
        for j <- 1..nobsy, into: [] do
          dist_func.(xx |> Matrex.column(i), yy |> Matrex.column(j), weights)
        end
    end
    |> Matrex.new()
  end

  def diff_conv(method, %Matrex{} = xx, %Matrex{} = yy, weights) do

    {dimx, nobsx} = Matrex.size(xx)
    {dimy, nobsy} = Matrex.size(yy)
    weights = weights || Matrex.ones(dimx, 1)
    {dimw, nobsw} = Matrex.size(weights)

    (dimx == dimy == dimw) || %ArgumentError{message: "size(xx, 1) != size(yy, 1) != size(weights, 1)"}
    (nobsx == nobsy == nobsw) || %ArgumentError{message: "nobs(xx) != nobs(yy) != nobs(weights)"}

    for idx <- 1..dimx do
      diff_conv(method, xx[idx], yy[idx])
    end
  end

  defp to_jaccard_distance({intersection, union}) do
    1 - intersection / union
  end
end
