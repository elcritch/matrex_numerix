defmodule SensorComputations.GP.Mat11Ard do
  alias __MODULE__

  defstruct sigma2: 0.0, scale: 0.0
  @type t :: %Mat11Ard{ sigma2: float, scale: float }

end

# defimpl SensorComputations.GP.Calc, for: SensorComputations.GP.Mat11 do

#   def calc(kern, xdata) do
#     0.0
#   end

#   @spec cov(SensorComputations.GP.Mat11.t(), float) :: float
#   def cov(mat, r), do: mat.sigma2 * :math.exp( -r / mat.scale)

#   def dk_dll(mat, r, wdiffp), do: wdiffp / r / mat.scale * cov(mat, r)

#   @spec dKij_dθp( Matrex.t() , Matrex.t(), integer, integer, integer, integer)
#   def dKij_dθp(mat, X, i, j, p, dim) do
#     r = distij(metric(mat),X,i,j,dim)
#     if p <= dim do
#         wdiffp=dist2ijk(metric(mat),X,i,j,p)
#         if wdiffp == 0.0, do: 0.0, else: dk_dll(mat,r,wdiffp)
#     else if p == dim+1 do
#         dk_dlσ(mat, r)
#     else
#         :na
#     end
#   end

#   # def dKij_dθp(mat::MaternARD, X::MatF64, data::StationaryARDData, i::Int, j::Int, p::Int, dim::Int) do
#   def dKij_dθp(mat::MaternARD, X::MatF64, data::StationaryARDData, i::Int, j::Int, p::Int, dim::Int) do
#       return dKij_dθp(mat,X,i,j,p,dim)
#   end

#   # def dk_dθp(mat::MaternIso, r::Float64, p::Int) do
#   def dk_dθp(mat::MaternIso, r::Float64, p::Int) do
#       if p==1 do
#           if r == 0.0, do: 0.0, else: dk_dll(mat, r)
#       else if p == 2 do
#           return dk_dlσ(mat, r)
#       else
#           return NaN
#       end
#   end

# end
