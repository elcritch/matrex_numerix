defmodule MatrexNumerix.FftTest do
  use ExUnit.Case, async: true


  test "test fft produces accurate numbers" do
    nn = 500
    tt = nn / 0.5
    t = Enum.to_list(1..nn) |> Matrex.from_list() |> Matrex.divide(tt) # np.linspace(0, 0.5, 500)
    tt = t[2] - t[1]  # sampling interval

    s = t |> Matrex.apply(fn t ->
      :math.sin(40 * 2 * :math.pi * t) + 0.5 * :math.sin(90 * 2 * :math.pi * t)
    end)

    fft = s |> MatrexNumerix.Fft.dft_real()

    # IO.inspect(t, label: :t)
    # IO.inspect(s, label: :s)
    # IO.inspect(fft, label: :fft)

    freq_amp =
      fft
      |> MatrexNumerix.Fft.to_freq_and_amplitude(tt)

    # assert Matrex.max(freq_amp[:amplitude]) == 0.50
    # assert Matrex.max(freq_amp[:amplitude][32..250]) == 0.251
    assert abs(Matrex.max(freq_amp[:amplitude]) - 0.5) < 1.0e-5
    assert abs(Matrex.max(freq_amp[:amplitude][32..250]) - 0.25) < 1.0e-5
  end
end
