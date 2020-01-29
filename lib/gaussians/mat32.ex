defmodule MatrexNumerix.GP.Mat32Iso do
  @moduledoc """
    Matern 32 Isotropic
  """
  alias __MODULE__
  use Matrex.Operators

  defstruct sigma2: 0.0, scale: 0.0

  @type t :: %Mat32Iso{ sigma2: float, scale: float }

  @sqrt3 :math.sqrt(3.0)

  # In cov(mat, r) at /Users/elcritch/.julia/packages/GaussianProcesses/9jiDV/src/kernels/mat32_iso.jl:41
  def cov(mat = %Mat32Iso{}, r) when is_number(r) do
    s = @sqrt3 * r / mat.scale
    mat.sigma2 * (1 + s) * :math.exp(-s)
  end

end

# defmodule MatrexNumerix.GP.Mat32Ard do
#   @moduledoc """
#     Matern 32 Isotropic
#   """
#   alias __MODULE__
#   use Matrex.Operators

#   defstruct sigma2: 0.0, scale: 0.0

#   @type t :: %Mat32Ard{ sigma2: float, scale: float }

#   @sqrt3 :math.sqrt(3.0)

#   # In cov(mat, r) at /Users/elcritch/.julia/packages/GaussianProcesses/9jiDV/src/kernels/mat32_iso.jl:41
#   def cov(mat = %Mat32Ard{}, r) do
#     s = √3 * r / mat.ℓ
#      mat.σ2 * (1 + s) * exp(-s)
#   end

# end
