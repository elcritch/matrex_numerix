defmodule MatrexNumerix.LinearRegression do
  @moduledoc """
  Linear regression functions.

  Taken from `MatrexNumerix` project at https://github.com/safwank/MatrexNumerix under the MIT license.
  """

  import Matrex.Guards
  import MatrexNumerix.Statistics

  alias MatrexNumerix.Correlation

  @doc """
  Least squares best fit for points `{x, y}` to a line `y:x↦a+bx`
  where `x` is the predictor and `y` the response.
  Returns a tuple containing the intercept `a` and slope `b`.
  """
  @spec fit( %Matrex{}, %Matrex{} ) :: {float, float} | no_return()


  def fit(
        vector_data(columns1, _data1, _x) = x,
        vector_data(columns2, _data2, _y) = y
      ) when columns1 == columns2 do

    x_mean = mean(x)
    y_mean = mean(y)
    variance = variance(x)
    covariance = covariance(x, y)
    slope = covariance / variance
    intercept = y_mean - slope * x_mean
    {intercept, slope}
  end

  def fit(%Matrex{} = _x, %Matrex{} = _y),
      do: raise %ArgumentError{message: "mismatched sizes"}

  @doc """
  Estimates a response `y` given a predictor `x`
  and a set of predictors and responses, i.e.
  it calculates `y` in `y:x↦a+bx`.
  """
  @spec predict(number, Matrex.t(), Matrex.t()) :: number
  def predict(x, xs, ys) do
    {intercept, slope} = fit(xs, ys)
    intercept + slope * x
  end

  @doc """
  Measures how close the observed data are to
  the fitted regression line, i.e. how accurate
  the prediction is given the actual data.
  Returns a value between 0 and 1 where 0 indicates
  a prediction that is worse than the mean and 1
  indicates a perfect prediction.
  """
  @spec r_squared(Matrex.t(), Matrex.t()) :: float
  def r_squared(predicted, actual) do
    predicted
    |> Correlation.pearson(actual)
    |> squared
  end

  # defp squared(nil), do: nil
  defp squared(correlation), do: :math.pow(correlation, 2)
end
