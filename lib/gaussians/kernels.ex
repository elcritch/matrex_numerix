
defmodule MatrexNumerix.GP.Kernel do
  alias MatrexNumerix.LinearAlgebra
  alias MatrexNumerix.GP
  alias MatrexNumerix.GPE
  use Matrex.Operators

  defstruct [:distance_method]

  @eps 1.0e-32

  def update_mll(gpe = %GPE{}) do
    # [:kdata, :kmod, :kern, :xx, :y]
    {_cov, chol_upper, _chol_comb} = update_cK(gpe)
    chol_lower = chol_upper |> transpose()

    mu = GP.Mean.mean(gpe.mean, gpe.xx)
    y = gpe.y - mu

    #  alpha = gpe.cK \ y
    z = LinearAlgebra.forward_substitution(chol_lower, y)
    alpha = LinearAlgebra.backward_substitution(chol_upper, z)

    # Marginal log-likelihood
    # gp.mll = - (dot(y, gp.alpha) + logdet(gp.cK) + log2π * gp.nobs) / 2
    %{ gpe | y: y, alpha: alpha, cK: %{u: chol_upper, l: chol_lower} }
  end

  # def update_cK(x = %Matrex{}, kernel, logNoise, data = %KernelData{}) do
  def update_cK(%GPE{xx: x, kern: kernel, logNoise: logNoise, kdata: kdata}) do
    {_dim, nobs} = size(x)
    sigma_buffer = cov(kernel, x, x, kdata)
    noise = :math.exp(2*logNoise) + @eps

    x_idt = Matrex.eye(nobs)
    sigma_buffer = sigma_buffer + x_idt * noise

    {sigma_buffer, chol} = make_posdef!(sigma_buffer)
    # make chol symmetric: get upper triag w/ no diag
    lr_offdiag = LinearAlgebra.ones_upper_offdiag(nobs) |> Matrex.transpose()
    chol_combined = chol |> Matrex.add(sigma_buffer |> Matrex.multiply(lr_offdiag))

    {sigma_buffer, chol, chol_combined}
  end

  def cov(k = %{}, xx1 = %Matrex{}, xx2 = %Matrex{}, kdata = %GP.KernelData{}) do
      # if xx1 == xx2
        # cov(cK, k, xx1, data)
      # end

      {dim1, nobs1} = size(xx1)
      {dim2, nobs2} = size(xx2)

      (dim1==dim2) || throw(%ArgumentError{message: "xx1 and xx2 must have same dimension: #{inspect size(xx1)} // #{inspect size(xx2)}"})

      # IO.inspect([nobs1: nobs1, nobs2: nobs2], label: :COV_NOBS)

      {dim, _} = size(xx1)

      for i <- 1..nobs1, into: [] do
        for j <- 1..nobs2, into: [] do
          # cK[i,j] = cov_ij(k, xx1, xx2, data, i, j, dim)
          # IO.inspect({i, j, dim}, label: COV_DIMS)
          # IO.inspect({k, xx1, xx2, kdata}, label: COVS)
          cov_ij(k, xx1, xx2, kdata, i, j, dim)
        end
      end
      |> Matrex.new()
  end

  def cov_ij(kern = %{__struct__: kmod}, _x1 = %Matrex{}, _x2 = %Matrex{}, kdata = %GP.KernelData{}, i, j, _dim) do
      case kdata.rdata |> size() do
        {_, 1} ->
          kmod.cov(kern, kdata.rdata[j] )
        {1, _} ->
          kmod.cov(kern, kdata.rdata[j] )
        {_, _} ->
          kmod.cov(kern, kdata.rdata[i][j] )
      end
  end

  def make_posdef!(mm = %Matrex{}) do
    {m, n} = size(mm)
    m == n || raise(%ArgumentError{message: "Covariance matrix must be square"})

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
              diag_weight = 1.0e-6 * trace(m)/n * Matrex.eye(Matrex.size(m))
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
