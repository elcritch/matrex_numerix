defmodule MatrexNumerix.GP.LinIso do
  @moduledoc """
  Linear covariance Function
      Lin(ll = %Vector{}, lσ = number() )
  # Arguments
    - `scale = float`: vector of length scales (given on log scale)
  """

  alias __MODULE__
  use Matrex.Operators
  import Matrex
  alias MatrexNumerix.GP

  def kern_dist(x1, x2), do: GP.KernelData.dot(:square, x1, x2)

  defstruct scale: 0.0

  @type t :: %LinIso{scale: float}

  # cov(mat::Mat12Iso, r::Number) = mat.σ2 * exp(-r / mat.ℓ)

  def cov(mat = %LinIso{}, r) when is_number(r) do
    r / :math.exp(2*mat.scale)
  end

end
