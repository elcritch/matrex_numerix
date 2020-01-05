defmodule MatrexNumerix.OptimizationTest do
  use ExUnit.Case, async: true

  import MatrexNumerix.Optimization

  @tag iterations: 10
  test "genetic optimization returns the solution with the lowest cost" do
    cost_fun = fn(xs) -> Enum.sum(xs) end

    # for_all {min, max, count} in such_that({min_, max_, _} in
      # {non_neg_integer(), non_neg_integer(), non_neg_integer()} when min_ < max_) do
    data = [
      { 419, 877, 489},
      { 1281, 11828, 9934 },
    ]

    for {min, max, count} <- data do

      domain = [min..max] |> Stream.cycle |> Enum.take(count)

      genetic(domain, cost_fun) == [min] |> Stream.cycle |> Enum.take(count)
    end
  end
end
