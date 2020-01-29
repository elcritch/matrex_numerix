
defmodule MatrexNumerix.GP do
  @moduledoc """
  Linear regression functions.

  Taken from `MatrexNumerix` project at https://github.com/safwank/MatrexNumerix under the MIT license.
  """

  import Matrex

  alias MatrexNumerix.Correlation
  import Kernel, except: [-: 1, +: 2, -: 2, *: 2, /: 2, <|>: 2]
  import Matrex
  import Matrex.Operators

  def gpe(xx, y, mean = %{}, kernel = %{__struct__: kmod}, logNoise) do
    kdata = kmod.kernel_dist(xx, xx)
  end

end
