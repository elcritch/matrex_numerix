defmodule MatrexNumerix.GP.Mat32Ard do
  @moduledoc """
    Matern 32 Isotropic
  """

  alias __MODULE__
  import Matrex
  import MatrexNumerix.GP

  alias MatrexNumerix.Correlation
  import Kernel, except: [-: 1, +: 2, -: 2, *: 2, /: 2, <|>: 2]
  import Matrex
  import Matrex.Operators

  defstruct sigma2: 0.0, scale: 0.0

  @type t :: %Mat32Ard{ sigma2: float, scale: float }


  @sqrt3 :math.sqrt(3.0)

  def calc(kern, xdata) do
    0.0
  end

  @spec cov(GP.Mat32.t(), float) :: float
  def cov(mat, r), do: mat.sigma2 * (1.0 + @sqrt3 * r / mat.scale) * :math.exp(-@sqrt3 * r / mat.scale)

  # def dk_dll(mat, r, wdiffp), do: r / mat.scale * cov(mat, r)
  def dk_dll(mat, r, wdiffp), do: 3.0 * mat.σ2 * wdiffp * exp(-@sqrt3 * r / mat.scale)

  # @spec dKij_dθp( Matrex.t() , Matrex.t(), integer, integer, integer, integer)
  # def dKij_dθp(mat, X, i, j, p, dim) do
  #   r = distij(metric(mat),X,i,j,dim)
  #   if p <= dim do
  #       wdiffp=dist2ijk(metric(mat),X,i,j,p)
  #       if wdiffp == 0.0, do: 0.0, else: dk_dll(mat,r,wdiffp)
  #   else if p == dim+1 do
  #       dk_dlσ(mat, r)
  #   else
  #       :na
  #   end
  # end

  # def dKij_dθp(mat::MaternARD, X::MatF64, data::StationaryARDData, i::Int, j::Int, p::Int, dim::Int) do
  # def dKij_dθp(mat::MaternARD, X::MatF64, data::StationaryARDData, i::Int, j::Int, p::Int, dim::Int) do
      # return dKij_dθp(mat,X,i,j,p,dim)
  # end

  # # def dk_dθp(mat::MaternIso, r::Float64, p::Int) do
  # def dk_dθp(mat::MaternIso, r::Float64, p::Int) do
  #     if p==1 do
  #         if r == 0.0, do: 0.0, else: dk_dll(mat, r)
  #     else if p == 2 do
  #         return dk_dlσ(mat, r)
  #     else
  #         return NaN
  #     end
  # end

end