defmodule MatrexNumerix.CorrelationTest do
  use ExUnit.Case, async: true
  import ListHelper

  alias MatrexNumerix.Correlation

  test "pearson is nil when any vector is empty" do
    refute Correlation.pearson([], [1])
    refute Correlation.pearson([2], [])
    refute Correlation.pearson([], [])
  end

  test "pearson is nil when the vectors are not the same size" do
    refute Correlation.pearson([1, 2, 3], [4, 5, 6, 7])
  end

  test "pearson correlation is zero when the vectors are equal but variance is zero" do
    # for_all {x, len} in {int(), pos_integer()} do
    numbers = [ {34, Matrex.random(4, 1)}, {74092, Matrex.random(10, 1)} ]

    for {x, len} <- numbers  do
      xs = [x] |> Stream.cycle() |> Enum.take(len)

      assert Correlation.pearson(xs, xs) == 0.0
    end
  end

  test "pearson correlation is one when the vectors are equal and variance is not zero" do
    # for_all xs in non_empty(list(int())) do
    numbers = [ Matrex.random(4, 1), Matrex.random(10, 1) ]
    for xs <- numbers do
      xs = xs |> Enum.uniq() |> Matrex.from_list()
      assert Correlation.pearson(xs, xs) == 1.0
    end
  end

  test "pearson correlation is between -1 and 1" do
    # for_all {xs, ys} in {non_empty(list(int())), non_empty(list(int()))} do
    numbers = [ Matrex.random(4, 1), Matrex.random(10, 1) ]
    for {xs, ys} <- numbers do
      {xs, ys} = equalize_length(xs, ys)

      assert Correlation.pearson(xs, ys) |> between?(-1, 1)
    end
  end

  test "pearson correlation is correct for a specific dataset" do
    vector1 = DataHelper.read("Lew") |> Map.get(:data) |> Enum.take(200)
    vector2 = DataHelper.read("Lottery") |> Map.get(:data) |> Enum.take(200)

    assert Correlation.pearson(vector1, vector2) == -0.02947086158072648
  end

  # test "weighted pearson is nil when any vector is empty" do
  #   refute Correlation.pearson([], [1], [2])
  #   refute Correlation.pearson([1], [], [2])
  #   refute Correlation.pearson([1], [2], [])
  #   refute Correlation.pearson([], [], [])
  # end

  test "weighted pearson correlation with constant weights is consistent with unweighted correlation" do
    vector1 = DataHelper.read("Lew") |> Map.get(:data) |> Enum.take(200) |> Matrex.from_list()
    vector2 = DataHelper.read("Lottery") |> Map.get(:data) |> Enum.take(200) |> Matrex.from_list()
    weights = [2.0] |> Stream.cycle() |> Enum.take(200) |> Matrex.from_list()

    # for {v1,v2,w} <- Enum.zip([vector1 |> Enum.to_list(), vector2 |> Enum.to_list(), weights |> Enum.to_list()]), into: [] do
      # IO.puts "#{v1}, #{v2}, #{w}"
    # end
    IO.puts("\n")

    weighted_correlation = Correlation.pearson(vector1, vector2, weights) |> IO.inspect(label: :weighted_corr)
    unweighted_correlation = Correlation.pearson(vector1, vector2) |> IO.inspect(label: :unweighted_corr)

    IO.puts("\n")
    assert_in_delta(weighted_correlation, unweighted_correlation, 1.0e-9)
  end
end
