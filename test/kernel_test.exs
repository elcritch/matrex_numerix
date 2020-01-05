defmodule MatrexNumerix.KernelTest do
  use ExUnit.Case, async: true
  import ListHelper

  alias MatrexNumerix.Kernel

  test "rbf is between 0 and 1" do
    numbers = [ Matrex.random(4, 1), Matrex.random(10, 1) ]
    # for {xs, ys} in such_that({xxs, yys} in {non_empty(list(number())), non_empty(list(number()))}
    for {xs, ys} <- numbers do

      {xs, ys} = equalize_length(xs, ys)

      Kernel.rbf(xs, ys) |> between?(0, 1)
    end
  end
end
