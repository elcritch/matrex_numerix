defmodule MatrexNumerix.MathTest do
  use ExUnit.Case, async: true
  alias MatrexNumerix.Math

  test "nth root is the reverse of power" do
    numbers = [ Enum.random(1..1_000), Enum.random(1..1_000), ]
    # for_all {x, n} in such_that({xx, nn} in {int(), int()} when xx > 0 && nn > 0) do
    for {x, n} <- numbers do
      x |> Math.nth_root(n) |> :math.pow(n) |> Float.round() == x
    end
  end
end
