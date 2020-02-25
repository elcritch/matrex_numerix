defmodule MatrexNumerix.GP.KernelData do
  alias __MODULE__

  defstruct [:rdata, :ktype, :dtype]

  def dot(dtype, x1 = %Matrex{}, x2 = %Matrex{}) do
    %KernelData{
      rdata: Matrex.transpose(x1) |> Matrex.dot(x2),
      ktype: :dot,
      dtype: dtype
    }
  end

  def isotropic(dtype, x1 = %Matrex{}, x2 = %Matrex{}) do
    %KernelData{
      rdata: MatrexNumerix.Distance.diff_conv(dtype, x1, x2),
      ktype: :isotropic,
      dtype: dtype
    }
  end

  def ard(dtype, x1 = %Matrex{}, x2 = %Matrex{}) do
    %KernelData{
      rdata: MatrexNumerix.Distance.diff_conv(dtype, x1, x2),
      ktype: :ard,
      dtype: dtype
    }
  end

end
