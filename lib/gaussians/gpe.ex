
defmodule MatrexNumerix.GPE do
  alias __MODULE__
  use Matrex.Operators
  alias MatrexNumerix.GPE
  alias MatrexNumerix.GP

  @moduledoc """
  Linear regression functions.

  Taken from `MatrexNumerix` project at https://github.com/safwank/MatrexNumerix under the MIT license.
  """

  defstruct [:kdata, :kmod, :mean, :kern, :xx, :y, :logNoise, :alpha, :cK]

  def calculate(xx = %Matrex{}, y = %Matrex{}, mean = %{}, kernel = %{__struct__: kmod}, logNoise) do
    kdata = kmod.kern_dist(xx, xx)

    %GPE{kdata: kdata, kmod: kmod, kern: kernel, xx: xx, y: y, logNoise: logNoise, mean: mean}
    |> GP.Kernel.update_mll()
  end

  def predict_y(gpe = %GPE{}, x = %Matrex{}) do
    {mu, sigma2} = predict_f(gpe, x)
    {mu, sigma2 + :math.exp(2.0*gpe.logNoise)}
  end

  def predict_f(gpe = %GPE{}, x = %Matrex{}) do
    {dims, n} = size(x)
    {gpdims, _} = gpe.xx |> size()
    IO.inspect(dims, label: :dims)
    IO.inspect(gpdims, label: :gpdims)
    dims == gpdims || throw(%ArgumentError{message: "Gaussian Process object and input observations do not have consistent dimensions"})

    ## Calculate prediction for each point independently
    for k <- 1..n, reduce: {zeros(size(x)), zeros(size(x))} do
      {mu, sigma2} ->
        xk = x |> submatrix(1..dims, k..k)
        # xk = x |> submatrix(k..k, 1..dims)
        {mm, sig} = predictMVN(xk, gpe)

        mu = mu |> set_submatrix(1..1, k..k, mm)

        case sig do
          nil ->
            {mu, sigma2}
          sig ->
            {mu, sigma2 |> set_submatrix(1..1, k..k, sig)}
        end
    end
  end

  def predictMVN(xpred, gpe = %GPE{xx: xtrain , y: ytrain, kern: kern = %{__struct__: kmod}, mean: mf, alpha: alpha, cK: cK}) do
    crossdata = kmod.kern_dist(xtrain, xpred)
    crossdata = %{ crossdata | rdata: crossdata.rdata}
    priordata = kmod.kern_dist(xpred, xpred)
    priordata = %{ priordata | rdata: priordata.rdata }
    # crossdata = KernelData(kernel, xtrain, xpred)
    # priordata = KernelData(kernel, xpred, xpred)

    # 4] crossdata: GaussianProcesses.IsotropicData{Array{Float64,2}}([4.098969039517268; 4.420885831075749; 3.7849485140879073; 0.2999681739134916; 2.7292281214539535; 3.6295838374541436; 1.4672366764292122; 0.0; 1.2384739340714845; 2.701117188015772])
    # 5] priordata: GaussianProcesses.IsotropicData{Array{Float64,2}}([0.0])
    IO.inspect(xpred, label: :XPRED)
    IO.inspect(xtrain, label: :XTRAIN)
    IO.inspect(ytrain, label: :YTRAIN)
    IO.inspect(crossdata, label: :CROSSDATA)
    IO.inspect(priordata, label: :PRIORDATA)

    kcross = GP.Kernel.cov(kern, xpred, xtrain, crossdata)
    IO.inspect(kcross, label: :KCROSS)
    kpred = GP.Kernel.cov(kern, xpred, xpred, priordata)
    mx = GP.Mean.mean(mf, xpred)
    predictMVN!(kpred, cK, kcross, mx, alpha)
  end

  def predictMVN!(kxx, kff, kfx, mx, alphaf) do
    IO.inspect(kxx, label: :PREDICTMVN_kxx)
    IO.inspect(kff, label: :PREDICTMVN_kff)
    IO.inspect(kfx, label: :PREDICTMVN_kfx)
    IO.inspect(mx, label: :PREDICTMVN_mx)
    IO.inspect(alphaf, label: :PREDICTMVN_alphaf)

      mu = mx |> add(kfx |> inner_dot(alphaf))
      # lck = whiten!(kff, kfx)
      # kxx = subtract_Lck(kxx, lck)
      # {mu, kxx}
      {mu, nil}
  end

  # def whiten!(a, x), do: whiten!(x, a, x)
  # def whiten!(r, chol_upper, x) do
  #     v = _rcopy!(r, x)
  #     ldiv!(transpose(chol_upper), v)
  # end
end
