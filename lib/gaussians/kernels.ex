
defmodule MatrexNumerix.GP.KernelData do
  alias __MODULE__
  alias MatrexNumerix.Correlation
  import Kernel, except: [-: 1, +: 2, -: 2, *: 2, /: 2, <|>: 2]
  import Matrex
  import Matrex.Operators

  # defstruct [k: %Isotropic, x1: %Matrex{}, x2: %Matrex{}]

end
