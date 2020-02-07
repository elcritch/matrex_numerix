defmodule MatrexNumerix.GP.Mat11Iso do
  @moduledoc """
  Matern 1/2 ARD covariance Function
      Mat12Ard(ll = %Vector{}, lσ = number() )
  # Arguments
    - `scale = %Vector{}`: vector of length scales (given on log scale)
    - `sigma2`: signal standard deviation (given on log scale)
  """

  alias __MODULE__
  use Matrex.Operators
  import Matrex
  alias MatrexNumerix.GP

  def kern_dist(x1, x2), do: GP.KernelData.isotropic(:euclidian, x1, x2)

  defstruct scale: 0.0, sigma2: 0.0

  @type t :: %Mat11Iso{sigma2: float, scale: float}

  # cov(mat::Mat12Iso, r::Number) = mat.σ2 * exp(-r / mat.ℓ)

  def cov(mat = %Mat11Iso{}, r) when is_number(r) do
    exp(-2.0 * mat.scale) * exp(-r / exp( mat.scale ))
  end

  def dk_dll(mat = %Mat11Iso{}, r, wdiffp) when is_number(r) and is_number(wdiffp) do
    wdiffp / r * cov(mat, r)
  end

end

defmodule MatrexNumerix.GP.Mat11Ard do
  @doc """
  Matern 1/2 ARD covariance Function
      Mat12Ard(ll::Vector{T}, lσ::T)
  # Arguments
    - `ll::Vector{Real}`: vector of length scales (given on log scale)
    - `lσ::Real`: signal standard deviation (given on log scale)
  """

  alias __MODULE__
  use Matrex.Operators
  import Matrex
  alias MatrexNumerix.GP

  def kern_dist(x1, x2), do: GP.KernelData.isotropic(:euclidian, x1, x2)

  defstruct scale: Vector.zeros(1), sigma2: 0.0, priors: Vector.zeros(1)

  @type t :: %Mat11Ard{sigma2: float, scale: Matrex.t(), priors: Matrex.t()}

  def cov(mat = %Mat11Ard{}, r) when is_number(r) do
    exp(2.0 * mat.sigma2) * exp(-r / exp(-2.0 * mat.scale))
  end

  def dk_dll(mat = %Mat11Ard{}, r, wdiffp) when is_number(r) and is_number(wdiffp) do
    wdiffp / r * cov(mat, r)
  end

end
