defmodule MatrexNumerix.Correlation do
  @moduledoc """
  Statistical correlation functions between two vectors.
  """

  import Matrex
  import Matrex.Guards
  alias MatrexNumerix.Statistics

  @doc """
  Calculates the Pearson correlation coefficient between two vectors.
  """
  @spec pearson(Matrex.t(), Matrex.t()) :: Common.maybe_float()

  def pearson(
        vector_data(columns1, _data1, _first) = x,
        vector_data(columns2, _data2, _second) = y
      ) when columns1 == columns2 do
    sum1 = sum(x) # |> IO.inspect(label: :sum_of_x))
    sum2 = sum(y) # |> IO.inspect(label: :sum_of_y))
    sum_of_squares1 = sum(pow(x, 2))
    sum_of_squares2 = sum(pow(y, 2))

    sum_of_products = inner_dot(x, y) |> Matrex.scalar()

    size = 1.0 * Enum.count(x)
    num = sum_of_products - (sum1 * sum2 / size)

    density =
      ((sum_of_squares1 - :math.pow(sum1, 2) / size) *
        (sum_of_squares2 - :math.pow(sum2, 2) / size))
      |> :math.sqrt()

    case density do
      0.0 -> 0.0
      _ -> num / density
    end
  end

  def pearson( %Matrex{}, %Matrex{}), do: raise %ArgumentError{message: "incorrect sizes"}

  @doc """
  Calculates the weighted Pearson correlation coefficient between two vectors.
  """
  @spec pearson(Matrex.t(), Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def pearson(vector1, vector2, weights) do
    weighted_covariance_xy = Statistics.weighted_covariance(vector1, vector2, weights)
    weighted_covariance_xx = Statistics.weighted_covariance(vector1, vector1, weights)
    weighted_covariance_yy = Statistics.weighted_covariance(vector2, vector2, weights)
    weighted_covariance_xy / :math.sqrt(weighted_covariance_xx * weighted_covariance_yy)
  end
end
