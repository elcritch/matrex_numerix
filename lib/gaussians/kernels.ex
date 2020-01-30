
defmodule MatrexNumerix.GP.Kernel do
  alias __MODULE__
  alias MatrexNumerix.Statistics
  alias MatrexNumerix.LinearAlgebra
  alias MatrexNumerix.GP
  alias MatrexNumerix.GPE
  use Matrex.Operators

  defstruct [:distance_method]

  @eps 1.0e-32

  def update_mll(gpe = %GPE{}) do
    # [:kdata, :kmod, :kern, :xx, :y]
    {cov, chol_upper, chol_comb} = update_cK(gpe)
    chol_lower = chol_upper |> transpose()

    mu = GP.Mean.mean(gpe.mean, gpe.xx)
    y = gpe.y - mu

    #  alpha = gpe.cK \ y
    z = LinearAlgebra.forward_substitution(chol_lower, y |> transpose())
    alpha = LinearAlgebra.backward_substitution(chol_upper, z) |> transpose()

    # Marginal log-likelihood
    # gp.mll = - (dot(y, gp.alpha) + logdet(gp.cK) + log2π * gp.nobs) / 2
    %{ gpe | y: y, alpha: alpha, cK: [u: chol_upper, l: chol_lower] }
  end

  # def update_cK(x = %Matrex{}, kernel, logNoise, data = %KernelData{}) do
  def update_cK(%GPE{xx: x, kern: kernel, logNoise: logNoise, kdata: kdata}) do
    {dim, nobs} = size(x)
    sigma_buffer = cov(kernel, x, x, kdata)
    noise = :math.exp(2*logNoise) + @eps

    x_idt = Matrex.eye(nobs)
    sigma_buffer = sigma_buffer + x_idt * noise

    {sigma_buffer, chol} = make_posdef!(sigma_buffer)
    # make chol symmetric: get upper triag w/ no diag
    lu_offdiag = LinearAlgebra.ones_upper_offdiag(nobs)
    lr_offdiag = LinearAlgebra.ones_upper_offdiag(nobs) |> Matrex.transpose()
    # chol = chol |> Matrex.add(chol |> Matrex.multiply(lu_offdiag) |> Matrex.transpose())
    chol_combined = chol |> Matrex.add(sigma_buffer |> Matrex.multiply(lr_offdiag))

    {sigma_buffer, chol, chol_combined}
  end

  def cov(k = %{}, xx1 = %Matrex{}, xx2 = %Matrex{}, kdata = %GP.KernelData{}) do
      # if xx1 == xx2
        # cov(cK, k, xx1, data)
      # end

      {dim1, nobs1} = size(xx1)
      {dim2, nobs2} = size(xx2)

      (dim1==dim2) || throw(%ArgumentError{message: "xx1 and xx2 must have same dimension"})

      IO.inspect([nobs1: nobs1, nobs2: nobs2], label: :COV_NOBS)

      {dim, _} = size(xx1)

      for i <- 1..nobs1, into: [] do
        for j <- 1..nobs2, into: [] do
          # cK[i,j] = cov_ij(k, xx1, xx2, data, i, j, dim)
          IO.inspect({i, j, dim}, label: COV_DIMS)
          IO.inspect({k, xx1, xx2, kdata}, label: COVS)
          cov_ij(k, xx1, xx2, kdata, i, j, dim)
        end
      end
      |> Matrex.new()
  end

  def cov_ij(kern = %{__struct__: kmod}, x1 = %Matrex{}, x2 = %Matrex{}, kdata = %GP.KernelData{}, i, j, dim) do
      IO.inspect(kdata.rdata, label: COV_IJ_RDATA)
      IO.inspect({i,j}, label: COV_IJ_RDATA_IJ_IDX)
      # IO.inspect(kdata.rdata |> get(i,j), label: COV_IJ_RDATA_IJ)
      # kmod.cov(kern, kdata.rdata |> get(i,j))
      case kdata.rdata |> size() do
        {_, 1} ->
          kmod.cov(kern, kdata.rdata[j] |> IO.inspect(label: COV_IJ_RDATA_IJ_KMOD))
        {1, _} ->
          kmod.cov(kern, kdata.rdata[j] |> IO.inspect(label: COV_IJ_RDATA_IJ_KMOD))
        {_, _} ->
          kmod.cov(kern, kdata.rdata[i][j] |> IO.inspect(label: COV_IJ_RDATA_IJ_KMOD))
      end
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
              # IO.inspect(:success, label: :chol_state)
              {mm!, chol!}
            else
              # that wasn't (numerically) positive definite,
              # so let's add some weight to the diagonal
              # IO.inspect(:fail, label: :chol_state)
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
