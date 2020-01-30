defmodule MatrexNumerix.LinearAlgebraTest do
  use ExUnit.Case, async: true

  import MatrexNumerix.LinearAlgebra

  test "norm of a normal list (for coverage)" do
    assert norm(2, [1, 2, 3] |> Matrex.from_list()) == 3.7416573867739458
  end

  test "backward subst" do

    aa! = Matrex.new " 7 -2 1; 0 -3 -5; 0 0 4 "
    bb! = Matrex.new " 12; -7; -4 "

    res = MatrexNumerix.LinearAlgebra.backward_substitution(aa!, bb!)

    assert res == Matrex.new("3; 4; -1")
  end

  test "backward subst 4x4" do

    au = Matrex.new " 1 2 1 -1; 0 -4 1 7; 0 0 -2 1; 0 0 0 -1 "
    bu = Matrex.new " 5; 1; 1; 3 "

    res = MatrexNumerix.LinearAlgebra.backward_substitution(au, bu)

    assert res == Matrex.new("16; -6; -2; -3")
  end

  test "forward subst" do

    aa = Matrex.new [[1.0, 0.0, 0.0], [2.0, 1.0, 0.0], [-1.0, -3.0, 1.0]]
    bb = Matrex.new [[12.0], [17.0], [5.0]]

    res = MatrexNumerix.LinearAlgebra.forward_substitution(aa, bb)

    assert res == Matrex.new("12; -7; -4")
  end
end
