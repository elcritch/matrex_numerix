defmodule MatrexNumerix.Kernel do
  @moduledoc """
  Functions used as kernel methods for classification, regression and clustering.
  """

  import Matrex

  alias MatrexNumerix.Common

  @doc """
  Radial basis function used to approximate given functions.
  It is particularly useful for time series prediction and
  control of nonlinear systems. It can also be seen as a
  simple single-layer type of artificial neural network where
  the function acts as the activation function for the network.
  """
  @spec rbf(Common.vector(), Common.vector(), integer) :: Common.maybe_float()

  def rbf(x = %Matrex{}, y = %Matrex{}, gamma \\ 10) do
    p = pow(x |> Matrex.subtract(y), 2)
    len = sum(p)
    :math.exp(-gamma * len)
  end

end
