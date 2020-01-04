defmodule MatrexNumerix.Activations do
  @moduledoc """
  Activation functions for neural networks.
  """

  use MatrexNumerix.Matrex

  @doc """
  Computes the softmax of a Matrex.
  """
  @spec softmax(%Matrex{}) :: %Matrex{}
  def softmax(x = %Matrex{dims: dims}) when dims > 1 do
    m = max(x)
    e = exp(x - m)
    s = sum(e)
    e / s
  end

  @doc """
  Computes the softplus of a Matrex.
  """
  @spec softplus(%Matrex{}) :: %Matrex{}
  def softplus(x) do
    log(1 + exp(x))
  end

  @doc """
  Computes the softsign of a Matrex.
  """
  @spec softsign(%Matrex{}) :: %Matrex{}
  def softsign(x) do
    x / (1 + abs(x))
  end

  @doc """
  Computes the element-wise sigmoid of a Matrex.
  """
  @spec sigmoid(%Matrex{}) :: %Matrex{}
  def sigmoid(x = %Matrex{dims: 0}) do
    1 / (1 + exp(-x))
  end

  def sigmoid(x = %Matrex{dims: dims}) when dims > 0 do
    z = exp(x)
    z / (1 + z)
  end

  @doc """
  Computes the rectified linear unit of a Matrex.
  """
  @spec relu(%Matrex{}) :: %Matrex{}
  def relu(x) do
    max(0, x)
  end

  @doc """
  Computes the leaky rectified linear unit of a Matrex.
  """
  @spec leaky_relu(%Matrex{}, number) :: %Matrex{}
  def leaky_relu(x, alpha) when alpha != 0 do
    max(alpha * x, x)
  end

  @doc """
  Computes the exponential linear unit of a Matrex.
  """
  @spec elu(%Matrex{}, number) :: %Matrex{}
  def elu(x, alpha \\ 1.0) do
    t_apply(
      fn
        i when i >= 0 -> i
        i -> alpha * (:math.exp(i) - 1)
      end,
      x
    )
  end

  @doc """
  Computes the scaled exponential linear unit of a Matrex.
  """
  @spec selu(%Matrex{}) :: %Matrex{}
  def selu(x) do
    alpha = 1.6732632423543772848170429916717
    scale = 1.0507009873554804934193349852946
    scale * elu(x, alpha)
  end

  @doc """
  Computes the element-wise hyperbolic tangent of a Matrex.
  """
  @spec tanh(%Matrex{}) :: %Matrex{}
  def tanh(x) do
    Matrex.tanh(x)
  end
end
