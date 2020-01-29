defmodule MatrexNumerix.GP.KernelData do
  alias __MODULE__

  defstruct [:rdata, :ktype, :dtype]

  def isotropic(dtype, x1 = %Matrex{}, x2 = %Matrex{}) do
    %KernelData{
      rdata: MatrexNumerix.Distance.diff_conv(dtype, x1, x2),
      ktype: :isotropic,
      dtype: dtype
    }
  end

end
