
defmodule MatrexNumerix.GP.IsotropicKernelData do
  defstruct [:rdata]
end
defmodule MatrexNumerix.GP.ArdKernelData do
  defstruct [:data]
end

defmodule MatrexNumerix.GP.Kernel do
  alias __MODULE__
  alias MatrexNumerix.Statistics
  alias MatrexNumerix.GP.KernelData
  alias MatrexNumerix.GP.IsotropicKernelData
  use Matrex.Operators

  defstruct [:distance_method]

  @eps 1.0e-32

  # def update_mll!(gp = GPE; noise = Bool=true, domean = Bool=true, kern = Bool=true)
  #     update_cK!(gp)

  #     μ = Statistics.mean(gp.mean, gp.x)
  #     y = gp.y - μ
  #     gp.alpha = gp.cK \ y
  #     # Marginal log-likelihood
  #     gp.mll = - (dot(y, gp.alpha) + logdet(gp.cK) + log2π * gp.nobs) / 2
  #     gp
  # end


  # def update_cK(x = %Matrex{}, kernel, logNoise, data = %KernelData{}) do
  #     {dim, nobs} = size(x)

  #     cK = Matrex.zeros(nobsx, nobsx)

  #     sigma_buffer = cov(kernel, x, x, data)

  #     noise = :math.exp(2*logNoise) + @eps

  #     sigma_buffer = sigma_buffer + (Matrex.eye(nobs) * noise)
  #     # for i <- 1..nobs do
  #       # sigma_buffer[i,i] += noise
  #     # end

  #     sigma_buffer, chol = make_posdef!(sigma_buffer, cholfactors(cK))

  #     wrap_cK(cK, sigma_buffer, chol)
  # end

  def cov(cK = %Matrex{}, k = Kernel, xx1 = %Matrex{}, xx2 = %Matrex{}, kdata) do
      # if xx1 == xx2
        # cov(cK, k, xx1, data)
      # end

      {dim1, nobs1} = size(xx1)
      {dim2, nobs2} = size(xx2)

      (dim1==dim2) || throw(%ArgumentError{message: "xx1 and xx2 must have same dimension"})

      {dim, _} = size(xx1)

      {nobs1, nobs2} == size(cK) || throw(%ArgumentError{message: "cK has size $(size(cK)) xx1 $(size(xx1)) and xx2 $(size(xx2))"})

      for i <- 1..nobs1, into: [] do
        for j <- 1..nobs2, into: [] do
          # cK[i,j] = cov_ij(k, xx1, xx2, data, i, j, dim)
          cov_ij(k, xx1, xx2, kdata, i, j, dim)
        end
      end
      |> Matrex.new()
  end

  # def cov_sym(cK = AbstractMatrix, k = Kernel, xx = AbstractMatrix, kdata) do
  #     {dim, nobs} = size(xx)

  #     {nobs, nobs} == size(cK) || throw(%ArgumentError{message: "cK has size $(size(cK)) and xx has size $(size(xx))"})

  #     # for j <- 1..nobs do
  #     #     cK[j,j] = cov_ij(k, xx, xx, kdata, j, j, dim)
  #     #     for i <- 1..j-1 do
  #     #         cK[i,j] = cov_ij(k, xx, xx, kdata, i, j, dim)
  #     #         cK[j,i] = cK[i,j]
  #     #     end
  #     # end
  #     for j <- 1..nobs do
  #         cK[j,j] = cov_ij(k, xx, xx, kdata, j, j, dim)
  #         for i <- 1..j-1 do
  #             cK[i,j] = cov_ij(k, xx, xx, kdata, i, j, dim)
  #             cK[j,i] = cK[i,j]
  #         end
  #     end
  # end

  def cov_ij(kern, x1 = %Matrex{}, x2 = %Matrex{}, data = %IsotropicKernelData{}, i, j, dim) do
      kern.__struct__.cov(kern, data.rdata[i][j])
  end

end
