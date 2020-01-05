defmodule MatrexNumerix.DistanceTest do
  use ExUnit.Case, async: true
  import ListHelper

  alias MatrexNumerix.{Distance, Correlation}

  test "MSE is correct for a specific example" do
    vector1 = [12, 15, 20, 22, 24] |> Matrex.from_list()
    vector2 = [13, 17, 18, 20, 24] |> Matrex.from_list()
    assert Distance.mse(vector1, vector2) == 2.6
  end

  test "MSE is 0 when the vectors are equal" do
    # for_all xs in non_empty(list(number())) do
    numbers = [ Matrex.random(4, 1), Matrex.random(10, 1) ]

    for xs <- numbers do
      assert Distance.mse(xs, xs) == 0
    end
  end

  test "MSE is not 0 when the vectors are different" do
    # for_all {xs, ys} in {non_empty(list(int())), non_empty(list(int()))} do

    data = [ {Matrex.random(4, 1), Matrex.random(4, 1)}, {Matrex.random(10, 1), Matrex.random(10, 1)}, ]

    for {xs, ys} <- data do

      if xs != ys do
        assert Distance.mse(xs, ys) != 0
      end
    end
  end

  test "RMSE is correct for a specific example" do
    vector1 = [7, 10, 12, 10, 10, 8, 7, 8, 11, 13, 10, 8] |> Matrex.from_list()
    vector2 = [6, 10, 14, 16, 7, 5, 5, 13, 12, 13, 8, 5] |> Matrex.from_list()
    assert Distance.rmse(vector1, vector2) == 2.9154759474226504
  end

  test "RMSE is 0 when the vectors are equal" do
    # for_all xs in non_empty(list(number())) do
    data = [ Matrex.random(4, 1), Matrex.random(10, 1) ]
    for xs <- data do
      assert Distance.rmse(xs, xs) == 0
    end
  end

  test "RMSE is not 0 when the vectors are different" do
    # for_all {xs, ys} in {non_empty(list(int())), non_empty(list(int()))} do
    data = [ {Matrex.random(4, 1), Matrex.random(4, 1)}, {Matrex.random(10, 1), Matrex.random(10, 1)}, ]

    for {xs, ys} <- data do

      if xs != ys do
        assert Distance.rmse(xs, ys) != 0
      end
    end
  end

  test "pearson distance is the inverse of its correlation" do
    # for_all {xs, ys} in {non_empty(list(int())), non_empty(list(int()))} do
    data = [ {Matrex.random(1, 4), Matrex.random(1, 4)}, {Matrex.random(1, 10), Matrex.random(1, 10)}, ]
    for {xs, ys} <- data do

      assert Distance.pearson(xs, ys) == 1.0 - Correlation.pearson(xs, ys)
    end
  end

  test "pearson distance is between 0 and 2" do
    # for_all {xs, ys} in {non_empty(list(int())), non_empty(list(int()))} do
    data = [ {Matrex.random(1, 4), Matrex.random(1, 4)}, {Matrex.random(1, 10), Matrex.random(1, 10)}, ]
    for {xs, ys} <- data do

      Distance.pearson(xs, ys) |> between?(0, 2)
    end
  end

  test "minkowski distance is 0 when the vectors are equal" do
    # for_all xs in non_empty(list(number())) do
    data = [ Matrex.random(4, 1), Matrex.random(10, 1) ]
    for xs <- data do
      assert Distance.minkowski(xs, xs) == 0
    end
  end

  test "minkowski distance is correct for a specific dataset when using the default lambda" do
    vector1 = [1, 3, 5, 6, 8, 9] |> Matrex.from_list()
    vector2 = [2, 5, 6, 6, 7, 7] |> Matrex.from_list()

    assert Distance.minkowski(vector1, vector2) == 2.6684016487219897
  end

  test "minkowski distance is correct for a specific dataset when using a different lambda" do
    vector1 = [1, 3, 5, 6, 8, 9] |> Matrex.from_list()
    vector2 = [2, 5, 6, 6, 7, 7] |> Matrex.from_list()
    lambda = 5

    assert Distance.minkowski(vector1, vector2, lambda) == 2.3185419629968713
  end

  test "euclidean distance is 0 when the vectors are equal" do
    # for_all xs in non_empty(list(number())) do
    data = [ Matrex.random(4, 1), Matrex.random(10, 1) ]
    for xs <- data do
      assert Distance.euclidean(xs, xs) == 0
    end
  end

  test "euclidean distance is correct for a specific dataset" do
    vector1 = [1, 3, 5, 6, 8, 9, 6, 4, 3, 2] |> Matrex.from_list()
    vector2 = [2, 5, 6, 6, 7, 7, 5, 3, 1, 1] |> Matrex.from_list()

    assert Distance.euclidean(vector1, vector2) == 4.2426406871196605
  end


  test "manhattan distance is 0 when the vectors are equal" do
    # for_all xs in non_empty(list(number())) do
    data = [ Matrex.random(4, 1), Matrex.random(10, 1) ]
    for xs <- data do
      assert Distance.manhattan(xs, xs) == 0
    end
  end

  test "manhattan distance is correct for a specific dataset" do
    vector1 = [1, 3, 5, 6, 8, 9, 6, 4, 3, 2] |> Matrex.from_list()
    vector2 = [2, 5, 6, 6, 7, 7, 5, 3, 1, 1] |> Matrex.from_list()

    assert Distance.manhattan(vector1, vector2) == 12
  end

  test "jaccard is correct for specific examples" do
    [
      {[0, 0.5], [0.5, 1], 1.0},
      {[4.5, 1], [4, 2], 1.0},
      {[1, 1, 1], [1, 1, 1], 0},
      {[2.5, 3.5, 3.0, 3.5, 2.5, 3.0], [3.0, 3.5, 1.5, 5.0, 3.5, 3.0], 0.6666666666666667},
      {[1, 3, 5, 6, 8, 9, 6, 4, 3, 2], [2, 5, 6, 6, 7, 7, 5, 3, 1, 1], 0.9}
    ]
    |> Enum.each(fn {vector1, vector2, distance} ->
      assert Distance.jaccard(vector1 |> Matrex.from_list(), vector2 |> Matrex.from_list()) == distance
    end)
  end

  test "jaccard is between 0 and 1" do
    # for_all {xs, ys} in {non_empty(list(non_neg_integer())), non_empty(list(non_neg_integer()))} do
    data = [ {Matrex.random(4, 1), Matrex.random(4, 1)}, {Matrex.random(10, 1), Matrex.random(10, 1)}, ]
    for {xs, ys} <- data do
      Distance.jaccard(xs, ys) |> between?(0, 1)
    end
  end
end
