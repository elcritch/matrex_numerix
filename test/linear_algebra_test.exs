defmodule MatrexNumerix.LinearAlgebraTest do
  use ExUnit.Case, async: true

  import MatrexNumerix.LinearAlgebra

  test "norm of a normal list (for coverage)" do
    assert norm(2, [1, 2, 3] |> Matrex.from_list()) == 3.7416573867739458
  end

  test "backward subst" do

    aa! = Matrex.new " 7 -2 1; 0 -3 -5; 0 0 4 "
    bb! = Matrex.new " 12.0 ; -7.0 ; -4.0 "

    res = MatrexNumerix.LinearAlgebra.backward_substitution(aa!, bb!)

    assert res == Matrex.from_list([3.0, 4.0, -1.0])
  end
end
