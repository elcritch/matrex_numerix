
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

  # def cov(cK = %Matrex{}, k = %GaussianProcesses.Kernel{}, xx = %Matrex{}, data = %GaussianProcesses.KernelData{}) do
  #    # Turn off original operators
  #   {dim, nobsv} = Matrex.size(xx)

  #   for j <- 1..nobsv, reduce: cK do
  #     Matrex.set(cK, j, j, cov_ij(k, xx, data, j, j, dim))
  #     for i <- 1..j-1, reduce: cK do
  #       ck_ij = k.__struct__.cov_ij(k, xx, data, i, j, dim)

  #       ck
  #       |> Matrex.set(i, j, ck_ij)
  #       |> Matrex.set(j, i, ck_ij)
  #     end
  #   end
  # end

end
