defmodule MatrexNumerix.StatisticsTest do
  use ExUnit.Case, async: false
  # use ExCheck
  import ListHelper

  alias MatrexNumerix.Statistics

  test "median is the middle value of a sorted list" do
    numbers = [ Matrex.random(4, 1), Matrex.random(10, 1), ]

    for xs <- numbers do
      xs = Enum.uniq(xs) |> Matrex.from_list()
      median = Statistics.median(xs)
      {first, second} = xs |> Enum.sort() |> Enum.split_while(fn x -> x <= median end)
      Enum.count(first) == Enum.count(second) or Enum.count(first) - 1 == Enum.count(second)
    end
  end

  test "median is between mix and max" do
    numbers = [ Matrex.random(4, 1), Matrex.random(10, 1), ]

    for xs <- numbers do
      Statistics.median(xs) |> between?(Matrex.min(xs), Matrex.max(xs))
    end
  end


  test "mode is nil if no value is repeated" do
    numbers = [ Matrex.random(4, 1), Matrex.random(10, 1) ]
    for xs <- numbers do
      xs |> Enum.to_list() |> Enum.uniq() |> Matrex.from_list() |> Statistics.mode() == nil
    end
  end

  test "mode is the most frequent value" do
    numbers = [ {3.4, Matrex.random(4, 1)}, {23.968, Matrex.random(10, 1)} ]

    # for {x, xs} in {number(), non_empty(list(number()))} do
    for {x, xs} <- numbers do
      frequent = [x]
      frequent_list = frequent |> Stream.cycle() |> Enum.take(Enum.count(xs) + 1)
      xs |> Enum.concat(frequent_list) |> Enum.shuffle() |> Matrex.from_list() |> Statistics.mode() == frequent
    end
  end

  test "mode is the most frequent set of values" do
    numbers = [ {0.17112, 0.34, Matrex.random(4, 1)}, {0.074092, 0.23968, Matrex.random(10, 1)} ]

    # for {x, y, xs} in such_that(
    #           {x_, y_, _} in {number(), number(), non_empty(list(number()))}
    #           when x_ < y_
    #         ) do
    for {x, y, xs} <- numbers do
      frequent_set = [x, y]
      frequent_list = frequent_set |> Stream.cycle() |> Enum.take(2 * (Enum.count(xs) + 1))

      xs
      |> Enum.reject(&Enum.member?(frequent_set, &1))
      |> Enum.concat(frequent_list)
      |> Enum.shuffle()
      |> Matrex.from_list()
      |> Statistics.mode()
      |> Enum.sort() == frequent_set
    end
  end


  # property "range is the difference between the largest and smallest values" do
  #   for_all xs in non_empty(list(number())) do
  #     sorted_xs = Enum.sort(xs)
  #     Statistics.range(xs) == List.last(sorted_xs) - List.first(sorted_xs)
  #   end
  # end


  test "variance is nil when list has only one element" do
    assert_raise ArgumentError, fn ->
      Statistics.variance([42] |> Matrex.from_list)
    end
  end

  # property "variance is the square of standard deviation" do
  #   for_all xs in such_that(xxs in non_empty(list(number())) when length(xxs) > 1) do
  #     xs |> Statistics.variance() |> Float.round(4) ==
  #       xs |> Statistics.std_dev() |> :math.pow(2) |> Float.round(4)
  #   end
  # end


  # property "population variance is the square of population standard deviation" do
  #   for_all xs in such_that(xxs in non_empty(list(number())) when length(xxs) > 1) do
  #     xs |> Statistics.population_variance() |> Float.round(4) ==
  #       xs |> Statistics.population_std_dev() |> :math.pow(2) |> Float.round(4)
  #   end
  # end


  test "std dev is nil when list has only one element" do
    assert_raise ArgumentError, fn ->
      Statistics.std_dev([42] |> Matrex.from_list())
    end
  end

  test "std dev is correct for specific datasets" do
    dataset1 = DataHelper.read("Lew")
    dataset2 = DataHelper.read("Lottery")

    xx =
      dataset1[:data]
      # |> IO.inspect(label: :std_dev_data)
      |> Enum.to_list()
      |> Matrex.from_list()
      # |> IO.inspect(label: :std_dev_data_1)
      |> Statistics.std_dev()

    assert_in_delta(
      xx,
      dataset1[:std_dev],
      0.0001
    )

    assert_in_delta(
      dataset2[:data] |> Enum.to_list() |> Statistics.std_dev(),
      dataset2[:std_dev],
      0.0001
    )
  end

  test "moment for a normal list (for coverage)" do
    assert Statistics.moment([1, 2, 3] |> Matrex.from_list(), 10) == 0.6666666666666666
  end


  test "kurtosis is correct for specific datasets" do
    dataset1 = DataHelper.read("Lew")
    dataset2 = DataHelper.read("Lottery")

    assert_in_delta(
      dataset1[:data] |> Enum.to_list() |> Matrex.from_list() |> Statistics.kurtosis(),
      -1.49604979214447,
      0.01
    )

    assert_in_delta(
      dataset2[:data] |> Enum.to_list() |> Matrex.from_list() |> Statistics.kurtosis(),
      -1.19256091074856,
      0.01
    )
  end


  test "skewness is correct for specific datasets" do
    dataset1 = DataHelper.read("Lew")
    dataset2 = DataHelper.read("Lottery")

    assert_in_delta(
      dataset1[:data] |> Enum.to_list() |> Matrex.from_list() |> Statistics.skewness(),
      -0.050606638756334,
      0.001
    )

    assert_in_delta(
      dataset2[:data] |> Enum.to_list() |> Matrex.from_list() |> Statistics.skewness(),
      -0.09333165310779,
      0.001
    )
  end


  test "covariance is nil when the list lengths do not match" do
    assert_raise ArgumentError, fn ->
     Statistics.covariance([1, 2] |> Matrex.from_list(), [3, 4, 5] |> Matrex.from_list())
    end
  end

  # property "covariance is consistent with variance" do
  #   for_all xs in such_that(xxs in non_empty(list(number())) when length(xxs) > 1) do
  #     assert_in_delta(Statistics.covariance(xs, xs), Statistics.variance(xs), 0.0000000001)
  #   end
  # end

  # property "covariance is symmetric" do
  #   for_all {xs, ys} in such_that(
  #             {xxs, yys} in {non_empty(list(number())), non_empty(list(number()))}
  #             when length(xxs) > 1 and length(yys) > 1
  #           ) do
  #     {xs, ys} = equalize_length(xs, ys)

  #     Statistics.covariance(xs, ys) == Statistics.covariance(ys, xs)
  #   end
  # end

  test "population covariance is nil when the list lengths do not match" do
    assert_raise ArgumentError, fn ->
      Statistics.population_covariance([1, 2] |> Matrex.from_list(), [3, 4, 5] |> Matrex.from_list())
    end
  end

  # property "population covariance is consistent with population variance" do
  #   for_all xs in such_that(xxs in non_empty(list(number())) when length(xxs) > 1) do
  #     assert_in_delta(
  #       Statistics.population_covariance(xs, xs),
  #       Statistics.population_variance(xs),
  #       0.0000000001
  #     )
  #   end
  # end

  # property "population covariance is symmetric" do
  #   for_all {xs, ys} in such_that(
  #             {xxs, yys} in {non_empty(list(number())), non_empty(list(number()))}
  #             when length(xxs) > 1 and length(yys) > 1
  #           ) do
  #     {xs, ys} = equalize_length(xs, ys)

  #     Statistics.population_covariance(xs, ys) == Statistics.population_covariance(ys, xs)
  #   end
  # end


  test "quantile is nil when tau is invalid" do
    refute Statistics.quantile([1, 2, 3] |> Matrex.from_list(), -0.1)
    refute Statistics.quantile([1, 2, 3] |> Matrex.from_list(), -1.1)
  end

  # property "quantile is between mix and max values" do
  #   for_all {xs, tau} in {non_empty(list(number())), int(0, 100)} do
  #     tau = tau / 100
  #     {minimum, maximum} = Enum.min_max(xs)

  #     Statistics.quantile(xs, tau) |> between?(minimum, maximum)
  #   end
  # end

  test "quantile is correct for specific examples" do
    xs = [-1, 5, 0, -3, 10, -0.5, 4, 0.2, 1, 6] |> Matrex.from_list()

    [
      {0, -3},
      {1, 10},
      {0.5, 3 / 5},
      {0.2, -4 / 5},
      {0.7, 137 / 30},
      {0.01, -3},
      {0.99, 10},
      {0.52, 287 / 375},
      {0.325, -37 / 240}
    ]
    |> Enum.each(fn {tau, expected} ->
      assert_in_delta(Statistics.quantile(xs, tau), expected, 0.0001)
    end)
  end


  test "percentile is nil when p is invalid" do
    refute Statistics.percentile([1, 2, 3], -1)
    refute Statistics.percentile([1, 2, 3], -101)
  end

  # property "percentile is between mix and max values" do
  #   for_all {xs, p} in {non_empty(list(number())), int(0, 100)} do
  #     {minimum, maximum} = Enum.min_max(xs)

  #     Statistics.percentile(xs, p) |> between?(minimum, maximum)
  #   end
  # end

  # property "percentile is consistent with quantile" do
  #   for_all {xs, p} in {non_empty(list(number())), int(0, 100)} do
  #     tau = p / 100

  #     Statistics.percentile(xs, p) == Statistics.quantile(xs, tau)
  #   end
  # end


  test "weighted mean is nil when the list lengths do not match" do
    assert_raise ArgumentError, fn ->
      Statistics.weighted_mean([1, 2] |> Matrex.from_list(), [3, 4, 5] |> Matrex.from_list())
    end
  end

  # property "weighted mean is consistent with arithmetic mean" do
  #   for_all {xs, w} in {non_empty(list(int())), pos_integer()} do
  #     weights = [w] |> Stream.cycle() |> Enum.take(length(xs))

  #     Statistics.weighted_mean(xs, weights) == Statistics.mean(xs)
  #   end
  # end

  test "weighted mean is correct for a specific dataset" do
    xs = [1, 3, 5, 6, 8, 9] |> Matrex.from_list()
    weights = [1.0, 0.8, 1.0, 0.9, 1.0, 0.66] |> Matrex.from_list()

    assert_in_delta(Statistics.weighted_mean(xs, weights), 5.175, 0.001)
  end
end
