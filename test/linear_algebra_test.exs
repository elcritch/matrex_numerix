defmodule MatrexNumerix.LinearAlgebraTest do
  use ExUnit.Case, async: true

  import MatrexNumerix.LinearAlgebra

  test "norm of a normal list (for coverage)" do
    assert norm(2, [1, 2, 3] |> Matrex.from_list()) == 3.7416573867739458
  end
end
