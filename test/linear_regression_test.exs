defmodule MatrexNumerix.LinearRegressionTest do
  use ExUnit.Case, async: true

  import ListHelper
  import MatrexNumerix.LinearRegression


  # test "least square fit is nil when the vectors aren't the same size" do
    # refute fit([1, 2, 3], [4, 5])
  # end

  test "least squares fit is correct for a specific example" do
    {intercept, slope} = fit(
      [1.3, 2.1, 3.7, 4.2] |> Matrex.from_list(),
      [2.2, 5.8, 10.2, 11.8] |> Matrex.from_list()
    )

    assert_in_delta(intercept, -1.52256, 0.00001)
    assert_in_delta(slope, 3.19383, 0.00001)
  end

  test "predict is correct for a specific example" do
    predictors = [1, 2.3, 3.1, 4.8, 5.6, 6.3] |> Matrex.from_list()
    responses = [2.6, 2.8, 3.1, 4.7, 5.1, 5.3] |> Matrex.from_list()

    assert predict(3.5, predictors, responses) == 3.7288688355893296
  end

  test "R^2 is between 0 and 1" do
    data = [
      {Matrex.rand(4,1), Matrex.rand(4,1)},
      {Matrex.rand(10,1), Matrex.rand(10,1)},
    ]

    for {xs, ys} <- data do
      {xs, ys} = equalize_length(xs, ys)

      r_squared(xs, ys) |> between?(0, 1)
    end
  end

  # test "R^2 is 1 when predicted perfectly matches actual" do
  #   data = [
  #     Matrex.rand(4,1),
  #     Matrex.rand(10,1),
  #   ]
  #   for xs <- data do
  #     implies length(xs) > 4 do
  #       r_squared(xs, xs) == 1
  #     end
  #   end
  # end

  test "R^2 is symmetric" do
    data = [
      {Matrex.rand(4,1), Matrex.rand(4,1)},
      {Matrex.rand(10,1), Matrex.rand(10,1)},
    ]
    for {xs, ys} <- data do
      {xs, ys} = equalize_length(xs, ys)

      r_squared(xs, ys) == r_squared(ys, xs)
    end
  end

  test "R^2 is correct for a specific example" do
    actual = [1, 2.3, 3.1, 4.8, 5.6, 6.3]
    predicted = [2.6, 2.8, 3.1, 4.7, 5.1, 5.3]

    assert r_squared(predicted, actual) == 0.9487852070867371
  end
end
