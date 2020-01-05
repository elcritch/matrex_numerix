defmodule MatrexNumerix.Distance do
  @moduledoc """
  Distance functions between two vectors.
  """

  import Matrex

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

  @doc """
  The Pearson's distance between two vectors.
  """
  @spec pearson(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def pearson(vector1, vector2) do
    case Correlation.pearson(vector1, vector2) do
      # _ -> nil
      correlation -> 1.0 - correlation
    end
  end

  @doc """
  The Minkowski distance between two vectors.
  """
  @spec minkowski(Matrex.t(), Matrex.t(), integer) :: Common.maybe_float()
  def minkowski(x = %Matrex{}, y = %Matrex{}, p) do
    norm(p, x |> Matrex.subtract(y))
  end

  def minkowski(vector1, vector2, p) do
    x = Matrex.new(vector1)
    y = Matrex.new(vector2)
    minkowski(x, y, p)
  end

  @doc """
  The Euclidean distance between two vectors.
  """
  @spec euclidean(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def euclidean(x = %Matrex{}, y = %Matrex{}) do
    l2_norm(x |> Matrex.subtract(y))
  end

  def euclidean(vector1, vector2) do
    x = Matrex.new(vector1)
    y = Matrex.new(vector2)
    euclidean(x, y)
  end

  @doc """
  The Manhattan distance between two vectors.
  """
  @spec manhattan(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def manhattan(x = %Matrex{}, y = %Matrex{}) do
    l1_norm(x |> Matrex.subtract(y) )
  end

  def manhattan(vector1, vector2) do
    x = Matrex.new(vector1)
    y = Matrex.new(vector2)
    manhattan(x, y)
  end

  @doc """
  The Jaccard distance (1 - Jaccard index) between two vectors.
  """
  @spec jaccard(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  # def jaccard(%Matrex{items: []}, %Matrex{items: []}), do: 0.0
  # def jaccard(%Matrex{items: []}, _), do: nil
  # def jaccard(_, %Matrex{items: []}), do: nil
  # def jaccard([], []), do: 0.0
  # def jaccard([], _), do: nil
  # def jaccard(_, []), do: nil

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

  defp to_jaccard_distance({intersection, union}) do
    1 - intersection / union
  end
end
