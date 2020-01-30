
defprotocol MatrexNumerix.GP.Mean do
  @doc "Calculate Kernel"
  def mean(gpm, xdata)
end

# MeanZero
defmodule MatrexNumerix.GP.MeanZero do
  defstruct []
end

defimpl MatrexNumerix.GP.Mean, for: MatrexNumerix.GP.MeanZero do
  def mean(_gpm = %MatrexNumerix.GP.MeanZero{}, xdata = %Matrex{}) do
    {_dim, nobs} = Matrex.size(xdata)
    Matrex.zeros(1, nobs)
  end
end

# MeanConst
defmodule MatrexNumerix.GP.MeanConst do
  defstruct [:const]
end

defimpl MatrexNumerix.GP.Mean, for: MatrexNumerix.GP.MeanConst do
  def mean(gpm = %MatrexNumerix.GP.MeanConst{}, xdata = %Matrex{}) do
    {_dim, nobs} = Matrex.size(xdata)
    Matrex.fill(1, nobs, gpm.const)
  end
end
