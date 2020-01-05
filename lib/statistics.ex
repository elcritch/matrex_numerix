defmodule MatrexNumerix.Statistics do
  @moduledoc """
  Common statistical functions.
  """

  import Matrex.Guards
  alias MatrexNumerix.Common

  @doc """
  The average of a list of numbers.
  """
  @spec mean(Matrex.t()) :: float

  def mean(x = %Matrex{}) do
    Matrex.sum(x) / Enum.count(x)
  end

  # def mean(xs) do
  #   x = Matrex.from_list(xs)
  #   mean(x)
  # end

  @doc """
  The middle value in a list of numbers.
  """
  @spec median(Matrex.t()) :: Common.maybe_float()

  def median(x = %Matrex{}) do
    middle_index = round(Enum.count(x) / 2) - 1
    x |> Enum.sort() |> Enum.at(middle_index)
  end

  def median(xs) do
    x = Matrex.new(xs)
    median(x)
  end

  @doc """
  The most frequent value(s) in a list.
  """
  @spec mode(Matrex.t()) :: float

  def mode(x = %Matrex{}) do
    counts =
      Enum.reduce(x, %{}, fn i, acc ->
        acc |> Map.update(i, 1, fn count -> count + 1 end)
      end)

    {_, max_count} = counts |> Enum.max_by(fn {_x, count} -> count end)

    case max_count do
      1 ->
        nil

      _ ->
        counts
        |> Stream.filter(fn {_x, count} -> count == max_count end)
        |> Enum.map(fn {i, _count} -> i end)
    end
  end

  @doc """
  The difference between the largest and smallest values in a list.
  """
  @spec range(Matrex.t()) :: Common.maybe_float()

  def range(x = %Matrex{}) do
    {minimum, maximum} = Enum.min_max(x)
    maximum - minimum
  end

  def range(xs) do
    x = Matrex.new(xs)
    range(x)
  end

  @doc """
  The unbiased population variance from a sample.
  It measures how far the vector is spread out from the mean.
  """
  @spec variance(Matrex.t()) :: Common.maybe_float()

  def variance(
        matrex_data(rows1, columns1, _data1, _first)
      ) when rows1 <= 1 and columns1 <= 1,
      do: raise %ArgumentError{message: "incorrect sizes"}

  def variance(x = %Matrex{}) do
    powered_deviations(x, 2) / (Enum.count(x) - 1)
  end

  @doc """
  The variance for a full population.
  It measures how far the vector is spread out from the mean.
  """
  @spec population_variance(Matrex.t()) :: Common.maybe_float()
  def population_variance(x = %Matrex{}), do: moment(x, 2)

  @doc """
  The unbiased standard deviation from a sample.
  It measures the amount of variation of the vector.
  """
  @spec std_dev(Matrex.t()) :: Common.maybe_float()
  def std_dev(
        matrex_data(rows1, columns1, _data1, _first)
      ) when rows1 <= 1 and columns1 <= 1,
      do: raise %ArgumentError{message: "incorrect sizes"}

  def std_dev(x = %Matrex{}), do: :math.sqrt(variance(x))

  @doc """
  The standard deviation for a full population.
  It measures the amount of variation of the vector.
  """
  @spec population_std_dev(Matrex.t()) :: float
  def population_std_dev(x = %Matrex{}), do: :math.sqrt(population_variance(x))

  @doc """
  The nth moment about the mean for a sample.
  Used to calculate skewness and kurtosis.
  """
  @spec moment(Matrex.t(), pos_integer) :: float
  def moment(_, 1), do: 0.0
  def moment(x = %Matrex{}, n), do: powered_deviations(x, n) / Enum.count(x)


  @doc """
  The sharpness of the peak of a frequency-distribution curve.
  It defines the extent to which a distribution differs from a normal distribution.
  Like skewness, it describes the shape of a probability distribution.
  """
  @spec kurtosis(Matrex.t()) :: float
  def kurtosis(x = %Matrex{}), do: moment(x, 4) / :math.pow(population_variance(x), 2) - 3


  @doc """
  The skewness of a frequency-distribution curve.
  It defines the extent to which a distribution differs from a normal distribution.
  Like kurtosis, it describes the shape of a probability distribution.
  """
  @spec skewness(Matrex.t()) :: float
  def skewness(x = %Matrex{}), do: moment(x, 3) / :math.pow(population_variance(x), 1.5)


  @doc """
  Calculates the unbiased covariance from two sample vectors.
  It is a measure of how much the two vectors change together.
  """
  @spec covariance(Matrex.t(), Matrex.t()) :: Common.maybe_float()
  def covariance(
        matrex_data(rows1, columns1, _data1, _first),
        matrex_data(rows2, columns2, _data2, _second)
      ) when rows1 <= 1 and columns1 <= 1 or rows2 <= 1 and columns2 <= 1,
      do: raise %ArgumentError{message: "incorrect sizes"}
  def covariance(
        matrex_data(rows1, columns1, _data1, _first),
        matrex_data(rows2, columns2, _data2, _second)
      ) when rows1 != rows2 or columns1 != columns2,
      do: raise %ArgumentError{message: "incorrect sizes"}

  def covariance(x = %Matrex{}, y = %Matrex{}) do
    divisor = Enum.count(x) - 1
    do_covariance(x, y, divisor)
  end

  def covariance(xs, ys) do
    x = Matrex.new(xs)
    y = Matrex.new(ys)
    covariance(x, y)
  end

  @doc """
  Calculates the population covariance from two full population vectors.
  It is a measure of how much the two vectors change together.
  """
  @spec population_covariance(Matrex.t(), Matrex.t()) :: Common.maybe_float()

  def population_covariance(
        matrex_data(rows1, columns1, _data1, _first),
        matrex_data(rows2, columns2, _data2, _second)
      ) when rows1 != rows2 or columns1 != columns2,
      do: raise %ArgumentError{message: "incorrect sizes"}

  def population_covariance(x = %Matrex{}, y = %Matrex{}) do
    divisor = Enum.count(x)
    do_covariance(x, y, divisor)
  end

  def population_covariance(xs, ys) do
    x = Matrex.new(xs)
    y = Matrex.new(ys)
    population_covariance(x, y)
  end

  @doc """
  Estimates the tau-th quantile from the vector.
  Approximately median-unbiased irrespective of the sample distribution.
  This implements the R-8 type of https://en.wikipedia.org/wiki/Quantile.
  """
  @spec quantile(Matrex.t(), number) :: Common.maybe_float()
  def quantile(_xs, tau) when tau < 0 or tau > 1, do: nil

  def quantile(x = %Matrex{}, tau) do
    sorted_x = Enum.sort(x)
    h = (length(sorted_x) + 1 / 3) * tau + 1 / 3
    hf = h |> Float.floor() |> round
    do_quantile(sorted_x, h, hf)
  end

  def quantile(xs, tau) do
    x = Matrex.new(xs)
    quantile(x, tau)
  end

  @doc """
  Estimates the p-Percentile value from the vector.
  Approximately median-unbiased irrespective of the sample distribution.
  This implements the R-8 type of https://en.wikipedia.org/wiki/Quantile.
  """
  @spec percentile(Matrex.t(), integer) :: Common.maybe_float()
  def percentile(_xs, p) when p < 0 or p > 100, do: nil
  def percentile(x = %Matrex{}, p), do: quantile(x, p / 100)

  def percentile(xs, p) do
    x = Matrex.new(xs)
    percentile(x, p)
  end

  @doc """
  Calculates the weighted measure of how much two vectors change together.
  """
  @spec weighted_covariance(Matrex.t(), Matrex.t(), Matrex.t()) ::
          Common.maybe_float()

  def weighted_covariance(
        matrex_data(columns1, _data1, _first),
        matrex_data(columns2, _data2, _second),
        matrex_data(columns3, _data3, _third)
      ) when columns1 != columns2 or columns1 != columns3,
      do: raise %ArgumentError{message: "incorrect sizes"}

  def weighted_covariance(x = %Matrex{}, y = %Matrex{}, w = %Matrex{} ) do
    weighted_covariance(x, y, w, :biased)
  end

  def weighted_covariance(x = %Matrex{}, y = %Matrex{}, w = %Matrex{}, :unbiased) do
    wm1 = weighted_mean(x, w)
    wm2 = weighted_mean(y, w)
    weighted_samples = w |> Matrex.multiply( Matrex.subtract(x, wm1)) |> Matrex.multiply( Matrex.subtract(y, wm2) )
    -1.0 * Matrex.sum(weighted_samples) / (1.0 - Matrex.sum(w |> Matrex.pow(2)))
  end
  def weighted_covariance(x = %Matrex{}, y = %Matrex{}, w = %Matrex{}, :biased) do
    wm1 = weighted_mean(x, w)
    wm2 = weighted_mean(y, w)
    weighted_samples = w |> Matrex.multiply( Matrex.subtract(x, wm1)) |> Matrex.multiply( Matrex.subtract(y, wm2) )
    1.0 * Matrex.sum(weighted_samples) / Matrex.sum(w)
  end

  @doc """
  Calculates the weighted average of a list of numbers.
  """
  @spec weighted_mean(Matrex.t(), Matrex.t()) :: Common.maybe_float()

  def weighted_mean(
    vector_data(columns1, _body1) = x,
    vector_data(columns2, _body2) = w
    ) when columns1 == columns2 do
    Matrex.sum(x |> Matrex.inner_dot(w)) / Matrex.sum(w)
  end
  def weighted_mean( %Matrex{}, %Matrex{}), do: raise %ArgumentError{message: "incorrect sizes"}

  @spec powered_deviations(Matrex.t(), number) :: float
  def powered_deviations(x, n) do
    x_mean = mean(x)
    Matrex.sum(Matrex.pow(x |> Matrex.subtract(x_mean), n))
  end

  defp do_covariance(x, y, divisor) do
    mean_x = mean(x)
    mean_y = mean(y)
    Matrex.sum((x |> Matrex.subtract(mean_x)) |> Matrex.dot_nt(y |> Matrex.subtract(mean_y))) / divisor
  end

  defp do_quantile([head | _], _h, hf) when hf < 1, do: head
  defp do_quantile(xs, _h, hf) when hf >= length(xs), do: List.last(xs)

  defp do_quantile(xs, h, hf) do
    Enum.at(xs, hf - 1) + (h - hf) * (Enum.at(xs, hf) - Enum.at(xs, hf - 1))
  end
end
