
defmodule MatrexNumerix.GP.Mean do
  defstruct [:kind]
end

defmodule MatrexNumerix.GP.Kernel do
  alias __MODULE__
  alias MatrexNumerix.Statistics
  alias MatrexNumerix.LinearAlgebra
  alias MatrexNumerix.GP.KernelData
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


  def update_cK(x = %Matrex{}, kernel, logNoise, data = %KernelData{}) do
    {dim, nobs} = size(x)
    sigma_buffer = cov(kernel, x, x, data)
    noise = :math.exp(2*logNoise) + @eps

    x_idt = Matrex.eye(nobs)
    sigma_buffer = sigma_buffer + x_idt * noise

    {sigma_buffer, chol} = make_posdef!(sigma_buffer)
    # make chol symmetric: get upper triag w/ no diag
    lr_offdiag = LinearAlgebra.ones_upper_offdiag(nobs) |> Matrex.transpose()
    chol = chol |> Matrex.add(sigma_buffer |> Matrex.multiply(lr_offdiag))

    {sigma_buffer, chol, chol_combined}
  end

  def cov(k = %{}, xx1 = %Matrex{}, xx2 = %Matrex{}, kdata = %KernelData{}) do
      # if xx1 == xx2
        # cov(cK, k, xx1, data)
      # end

      {dim1, nobs1} = size(xx1)
      {dim2, nobs2} = size(xx2)

      (dim1==dim2) || throw(%ArgumentError{message: "xx1 and xx2 must have same dimension"})

      {dim, _} = size(xx1)

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

  def cov_ij(kern = %{__struct__: kmod}, x1 = %Matrex{}, x2 = %Matrex{}, kdata = %KernelData{}, i, j, dim) do
      kmod.cov(kern, kdata.rdata[i][j])
  end

  def make_posdef!(mm = %Matrex{}) do
    {m, n} = size(mm)
    m == n || throw(%ArgumentError{message: "Covariance matrix must be square"})

    {mm!, chol} =
      for _ <- 1..10, reduce: {mm, :none} do # 10 chances
        {mm!, chol} ->

        case chol do
          :none ->
            # try cholesky ... assuming inf's mean it's no PD?
            chol! = Matrex.cholesky(mm!) |> Matrex.transpose()

            unless chol! |> Matrex.sum() == :nan do
              IO.inspect(:success, label: :chol_state)
              {mm!, chol!}
            else
              # that wasn't (numerically) positive definite,
              # so let's add some weight to the diagonal
              IO.inspect(:fail, label: :chol_state)
              diag_weight = 1.0e-6 * trace(m)/n * Matrew.eye(Matrex.size(m))
              {mm! + diag_weight, :none}
            end

          chol ->
            {mm!, chol}
        end
      end

    # copyto!(chol_factors, m)
    # chol = cholesky!(Symmetric(chol_factors, :U))
    {mm!, chol}
  end

end
