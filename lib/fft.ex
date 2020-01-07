defmodule MatrexNumerix.Fft do
  @moduledoc """
  Computes the discrete Fourier transform (DFT) of the given complex vector.
  """

  import Matrex.Guards

  def dft_freq_and_amplitude(
        vector_data(columns1, _body1) = x,
        vector_data(columns2, _body2) = y
      ) when columns1 == columns2 do
    sampling_interval = x[2] - x[1] # calculate sampling interval

    Matrex.concat(y, Matrex.zeros(1, columns1), :rows)
    |> dft_complex()
    |> to_freq_and_amplitude(sampling_interval)
  end

  @doc """
  Computes the discrete Fourier transform (DFT) of the given real vector.

    - `y` real and imaginary input

  Returns:
    - a matrix of real and imaginary fourier transform output
  """
  def dft_real(vector_data(columns1, _body1) = y) do
    Matrex.concat(y, Matrex.zeros(1, columns1), :rows)
    |> dft_complex()
  end

  @doc """
  Computes the discrete Fourier transform (DFT) of the given complex vector.

    - `xx` real and imaginary input

  Returns:
    - a matrix of real and imaginary fourier transform output
  """
  def dft_complex(xx) do
    {2, nn} = xx |> Matrex.size()

    # IO.inspect(xx, label: :xx)
    rx = xx[1]
    ry = xx[2]

    output! =
      for k <- 1..nn, reduce: Matrex.zeros(2, nn) do
        output ->
          {sumreal, sumimag} =
            for t <- 1..nn, reduce: {0.0, 0.0} do  # For each input element
              {sumreal, sumimag} ->
                angle = 2 * :math.pi * t * k / nn
                sumreal! = sumreal +  rx[t] * :math.cos(angle) + ry[t] * :math.sin(angle)
                sumimag! = sumimag + -rx[t] * :math.sin(angle) + ry[t] * :math.cos(angle)

                {sumreal!, sumimag!}
            end

          output
          |> Matrex.set(1, k, sumreal)
          |> Matrex.set(2, k, sumimag)
      end

    output!
  end

  def to_freq_and_amplitude(
        matrex_data(2, nn, _data1, _first) = fft,
        sampling_interval) do

    tt = sampling_interval
    # 1/T = frequency
    ff =
      (Enum.to_list(1..nn)
      |> Matrex.from_list()
      |> Matrex.apply(fn x -> (x-1) * 1/(nn*tt) end))[1..div(nn,2)]
    # ff = np.linspace(0, 1 / tt, nn)

    fft_amp =
      (Enum.map(1..nn, fn i -> fft |> Matrex.column_to_list(i) end)
      |> Enum.map(fn [r,i] -> :math.sqrt( :math.pow(r,2) + :math.pow(i,2)) end)
      |> Matrex.from_list()
      |> Matrex.divide(1.0*nn))[1..div(nn,2)]

    [frequency: ff, amplitude: fft_amp]
  end
end
