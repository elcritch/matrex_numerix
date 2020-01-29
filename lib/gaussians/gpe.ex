
defmodule MatrexNumerix.GuassianProcess.KernelData do
  alias __MODULE__
  alias MatrexNumerix.Correlation
  import Kernel, except: [-: 1, +: 2, -: 2, *: 2, /: 2, <|>: 2]
  import Matrex
  import Matrex.Operators

  # defstruct [k: %Isotropic, x1: %Matrex{}, x2: %Matrex{}]
  defstruct [:distance_data]

  def new(k = %Isotropic{}, x1 = %Matrex{}, x2 = %Matrex{}) do
    %KernelData{dists: distance(k, x1, x2)}
  end

  def metric(k) do
    0.0
  end

  def distance(k = %Stationary{}, xx = %Matrex{}, yy = %Matrex{}) do
    nobsx = Matrex.size(xx, 2)
    nobsy = Matrex.size(yy, 2)
    dist = Matrex.zeros(nobsx, nobsy)
    m = metric(k)

    {dimx, nobsx} = Matrex.size(xx)
    {dimy, nobsy} = Matrex.size(yy)
    dimx == dimy || %ArgumentError{message: "size(xx, 1) != size(yy, 1)"}

    dist =
      for i <- 1..nobsx, reduce: dist do
        dist ->
          for j <- 1..nobsy, reduce: dist do
            dist ->
              dist_ij = distij(m, xx, yy, i, j, dimx)
              dist |> Matrex.set(i, j, dist_ij)
          end
      end

    dist
  end

 def distij(dist = %Euclidean{}, x1 = %Matrex{}, x2 = %Matrex{}, i, j, dim) do
     if x1 == x2 && i == j do 0.0 else :math.sqrt(sq_euclidean_ij(x1, x2, i, j, dim)) end
 end

 def sq_euclidean_ijk(x1 = %Matrex{}, x2 = %Matrex{}, i, j, k), do: pow(x1[k][i] - x2[k][j], 2.0)

end

"""
20] kerneldata:
  GaussianProcesses.IsotropicData{Array{Float64,2}}(
    [0.0 0.32191679155848085 0.31402052542936065 3.799000865603776 1.369740918063314 0.4693852020631244 2.6317323630880556 4.098969039517268 2.860495105445783 1.3978518515014957;
     0.32191679155848085 0.0 0.6359373169878415 4.120917657162257 1.6916577096217948 0.7913019936216052 2.9536491546465364 4.420885831075749 3.182411897004264 1.7197686430599766;
     0.31402052542936065 0.6359373169878415 0.0 3.4849803401744155 1.0557203926339533 0.15536467663376374 2.317711837658695 3.7849485140879073 2.5464745800164224 1.083831326072135;
     3.799000865603776 4.120917657162257 3.4849803401744155 0.0 2.429259947540462 3.3296156635406517 1.1672685025157206 0.2999681739134916 0.9385057601579929 2.4011490141022804;
     1.369740918063314 1.6916577096217948 1.0557203926339533 2.429259947540462 0.0 0.9003557160001896 1.2619914450247416 2.7292281214539535 1.4907541873824692 0.028110933438181718;
     0.4693852020631244 0.7913019936216052 0.15536467663376374 3.3296156635406517 0.9003557160001896 0.0 2.162347161024931 3.6295838374541436 2.3911099033826586 0.9284666494383713;
     2.6317323630880556 2.9536491546465364 2.317711837658695 1.1672685025157206 1.2619914450247416 2.162347161024931 0.0 1.4672366764292122 0.22876274235772764 1.2338805115865599;
     4.098969039517268 4.420885831075749 3.7849485140879073 0.2999681739134916 2.7292281214539535 3.6295838374541436 1.4672366764292122 0.0 1.2384739340714845 2.701117188015772;
     2.860495105445783 3.182411897004264 2.5464745800164224 0.9385057601579929 1.4907541873824692 2.3911099033826586 0.22876274235772764 1.2384739340714845 0.0 1.4626432539442875;
     1.3978518515014957 1.7197686430599766 1.083831326072135 2.4011490141022804 0.028110933438181718 0.9284666494383713 1.2338805115865599 2.701117188015772 1.4626432539442875 0.0])
"""

defmodule MatrexNumerix.GuassianProcess do
  @moduledoc """
  Linear regression functions.

  Taken from `MatrexNumerix` project at https://github.com/safwank/MatrexNumerix under the MIT license.
  """

  import Matrex
  import MatrexNumerix.GaussianProcesses

  alias MatrexNumerix.Correlation
  import Kernel, except: [-: 1, +: 2, -: 2, *: 2, /: 2, <|>: 2]
  import Matrex
  import Matrex.Operators

  def cov(cK = %Matrex{}, k = %GaussianProcesses.Kernel{}, xx = %Matrex{}, data = %GaussianProcesses.KernelData{}) do
     # Turn off original operators
    {dim, nobsv} = Matrex.size(xx)

    for j <- 1..nobsv, reduce: cK do
      Matrex.set(cK, j, j, cov_ij(k, xx, data, j, j, dim))
      for i <- 1..j-1, reduce: cK do
        ck_ij = k.__struct__.cov_ij(k, xx, data, i, j, dim)

        ck
        |> Matrex.set(i, j, ck_ij)
        |> Matrex.set(j, i, ck_ij)
      end
    end
  end

end
