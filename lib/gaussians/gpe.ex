
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
    dims == gpdims || raise(%ArgumentError{message: "Gaussian Process object and input observations do not have consistent dimensions"})

    ## Calculate prediction for each point independently
    for k <- 1..n, reduce: {zeros(size(x)), zeros(size(x))} do
      {mu, sigma2} ->
        xk = x |> submatrix(1..dims, k..k)
        # xk = x |> submatrix(k..k, 1..dims)
        {mm, sig} = predictMVN(xk, gpe)

        mu = mu |> Vector.set_slice(k..k, mm)
        sigma2 = sigma2 |> Vector.set_slice(k..k, sig || Vector.new(" NaN "))

        {mu, sigma2}
    end
  end

  def predictMVN(xpred,
                 %GPE{xx: xtrain , y: ytrain,
                      kern: kern = %{__struct__: kmod},
                      mean: mf,
                      alpha: alpha, cK: cK}) do

    crossdata = kmod.kern_dist(xtrain, xpred)
    priordata = kmod.kern_dist(xpred, xpred)

    kcross = GP.Kernel.cov(kern, xpred, xtrain, crossdata)
    kpred = GP.Kernel.cov(kern, xpred, xpred, priordata)
    mx = GP.Mean.mean(mf, xpred)
    predictMVN!(kpred, cK, kcross, mx, alpha)
  end

  def predictMVN!(kxx, kff, kfx, mx, alphaf) do
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
