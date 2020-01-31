defmodule MatrexNumerix.LinearAlgebraTest do
  use ExUnit.Case, async: true

  import MatrexNumerix.LinearAlgebra
  alias Matrex.Vector

  test "norm of a normal list (for coverage)" do
    assert norm(2, [1, 2, 3] |> Matrex.from_list()) == 3.7416573867739458
  end

  test "backward subst" do

    aa! = Matrex.new(" 7 -2 1; 0 -3 -5; 0 0 4 ")
    bb! = Matrex.new " 12 -7 -4 "

    res = MatrexNumerix.LinearAlgebra.backward_substitution(aa!, bb!)

    assert res == Matrex.new("3 4 -1")
  end

  test "backward subst 4x4" do

    au = Matrex.new " 1 2 1 -1; 0 -4 1 7; 0 0 -2 1; 0 0 0 -1 "
    bu = Matrex.new " 5 1 1 3 "

    res = MatrexNumerix.LinearAlgebra.backward_substitution(au, bu)

    assert res == Matrex.new("16 -6 -2 -3")
  end

  test "forward subst" do

    aa = Matrex.new [[1.0, 0.0, 0.0], [2.0, 1.0, 0.0], [-1.0, -3.0, 1.0]]
    bb = Vector.new [12.0, 17.0, 5.0]

    res = MatrexNumerix.LinearAlgebra.forward_substitution(aa, bb)

    assert res == Matrex.new("12 -7 -4")
  end

  test "forward subst 4x4" do

    at = Matrex.new """
      3  0  0  0
      2  1  0  0
      1  0  1  0
      1  1  1  1
    """

    bt = Matrex.new " 4.0 2.0 4.0 2.0 "

    res = MatrexNumerix.LinearAlgebra.forward_substitution(at, bt)
    expected = Matrex.new " 1.3333333333333333 -0.6666666666666665 2.666666666666667 -1.3333333333333337 "

    # assert res == expected
    assert abs(res |> Matrex.subtract(expected) |> Matrex.sum()) < 1.0e-5
  end

  test "forward subst 5x5" do

    at = Matrex.new """
      1.0 0.0 0.0 0.0 0.0
      2.25 1.0 0.0 0.0 0.0
      -0.5 6.66667 1.0 0.0 0.0
      2.25 -7.66667 -1.24427 1.0 0.0
      -0.25 -3.66667 -0.473282 33.4951 1.0
    """

    bt = Matrex.new " -7 2 4 -4 -7 "

    res = MatrexNumerix.LinearAlgebra.forward_substitution(at, bt)
    expected = Matrex.new " -7.00 17.750 -117.83339 1.21684 -40.19329 "

    assert abs(res |> Matrex.subtract(expected) |> Matrex.sum()) < 1.0e-5
  end

  test "solve Ax=b " do

    # 1. Factor A = LU by Gaussian elimination (not including row swaps, discussed below!), giving Ax =
      # b =⇒ LUx = L(Ux) = b
      # 2. Let c = Ux. Solve Lc = b for c by forward-substitution.
      # 3. Solve Ux = c for x by backsubstitution.
    _aa = Matrex.new """
       4  -2  -7  -4  -8
       9  -6  -6  -1  -5
      -2  -9   3  -5   2
       9   7  -9   5  -8
      -1   6  -3   9   6
    """

    ll = Matrex.new """
       1.0        0.0        0.0       0.0        0.0
       1.0        1.0        0.0       0.0        0.0
       0.444444   0.0512821  1.0       0.0        0.0
      -0.111111   0.410256   0.582822  1.0        0.0
      -0.222222  -0.794872   0.171779  0.0242696  1.0
    """

    uu = Matrex.new """
      9.0  -6.0  -6.0      -1.0      -5.0
      0.0  13.0  -3.0       6.0      -3.0
      0.0   0.0  -4.17949  -3.86325  -5.62393
      0.0   0.0   0.0       8.67894   9.95297
      0.0   0.0   0.0       0.0      -0.771206
    """

    b = Matrex.new " -7 2 4 -4 -7 "

    c = MatrexNumerix.LinearAlgebra.forward_substitution(ll, b)

    expected_c = Matrex.new " -7.0 9.0 6.649572649572649 -12.345603271983638 -2.244344957587181 "

    x = MatrexNumerix.LinearAlgebra.backward_substitution(uu, c)

    # expected_x = Matrex.new " 0.5059578368469305 -0.9285059578368458 2.1640696608615944 1.4616559731133496 -1.2642835319278931 "
    expected_x = Matrex.new " 1.77543 3.305220 -1.10724 -4.75985 2.91017 "

    assert abs(c |> Matrex.subtract(expected_c) |> Matrex.sum()) < 1.0e-5
    # assert Matrex.transpose(x) == expected_x
    assert abs(x |> Matrex.subtract(expected_x) |> Matrex.sum()) < 1.0e-1

  end
end
