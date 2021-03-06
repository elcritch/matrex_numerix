defmodule MatrexNumerix.WindowTest do
  use ExUnit.Case, async: true

  import ListHelper

  alias MatrexNumerix.Window

  test "gaussian is correct for specific examples" do
    [
      {0.0, 1.0},
      {0.1, 0.99501247919268232},
      {1.0, 0.60653065971263342},
      {3.0, 0.01110899653824231}
    ]
    |> Enum.each(fn {width, expected} ->
      assert_in_delta(Window.gaussian(width), expected, 0.0000000001)
    end)
  end

  test "gaussian is between 0 and 1" do
    numbers = [ {:random.uniform(), Enum.random(1..1_000)},
                {:random.uniform(), Enum.random(1..1_000)} ]
    # for_all {width, sigma} in {number(), non_neg_integer()} do
    for {width, sigma} <- numbers do
      sigma = sigma / 100

      Window.gaussian(width, sigma) |> between?(0, 1)
    end
  end
end
