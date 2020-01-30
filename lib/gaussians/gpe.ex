
defmodule MatrexNumerix.GPE do
  alias __MODULE__
  alias MatrexNumerix.GP.Mean
  use Matrex.Operators

  @moduledoc """
  Linear regression functions.

  Taken from `MatrexNumerix` project at https://github.com/safwank/MatrexNumerix under the MIT license.
  """

  defstruct [:kdata, :kmod, :kern, :xx, :y]

  def calculate(xx = %Matrex{}, y = %Matrex{}, mean = %Mean{}, kernel = %{__struct__: kmod}, logNoise) do
    kdata = kmod.kern_dist(xx, xx)

    %GPE{kdata: kdata, kmod: kmod, kern: kernel, xx: xx, y: y}
  end

end
