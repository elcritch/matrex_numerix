defmodule MatrexNumerix.GP.Distance do

  @doc "Calculate distance"

  alias __MODULE__
  alias MatrexNumerix.Correlation
  import Kernel, except: [-: 1, +: 2, -: 2, *: 2, /: 2, <|>: 2]
  import Matrex
  import Matrex.Operators

  # defstruct [k: %Isotropic, x1: %Matrex{}, x2: %Matrex{}]
  # defstruct distance_data: %Matrex{}

  def isotropic(k, x1 = %Matrex{}, x2 = %Matrex{}) do
    # %KernelData{dists: distance(k, x1, x2)}
  end

  def metric(k = %{distance: :euclidian}) do
    0.0
  end

  # def distance(k = %{}, xx = %Matrex{}, yy = %Matrex{}) do
  #   nobsx = Matrex.size(xx, 2)
  #   nobsy = Matrex.size(yy, 2)
  #   dist = Matrex.zeros(nobsx, nobsy)
  #   m = metric(k)

  #   {dimx, nobsx} = Matrex.size(xx)
  #   {dimy, nobsy} = Matrex.size(yy)
  #   dimx == dimy || %ArgumentError{message: "size(xx, 1) != size(yy, 1)"}

  #   dist =
  #     for i <- 1..nobsx, reduce: dist do
  #       dist ->
  #         for j <- 1..nobsy, reduce: dist do
  #           dist ->
  #             dist_ij = dist_ij(m, xx, yy, i, j, dimx)
  #             dist |> Matrex.set(i, j, dist_ij)
  #         end
  #     end

  #   dist
  # end

end
