defmodule Numerix.Mixfile do
  use Mix.Project

  def project do
    [
      app: :numerix,
      name: "Numerix",
      description:
        "A collection of useful mathematical functions in Elixir with a slant towards statistics, linear algebra and machine learning",
      version: "0.5.1",
      elixir: "~> 1.5",
      source_url: "https://github.com/safwank/Numerix",
      deps: deps(),
      package: package(),
      test_coverage: [tool: ExCoveralls],
      preferred_cli_env: [
        credo: :test,
        "coveralls.html": :test,
        commit: :test
      ],
      aliases: [
        commit: ["dialyzer", "credo --strict", "coveralls.html --trace"]
      ],
      default_task: "commit"
    ]
  end

  def application do
    [applications: [:logger]]
  end

  defp deps do
    [
      {:matrex, "~> 0.6.8", github: "elcritch/matrex", branch: "add-numerix-ports"},
      {:ex_doc, "~> 0.18", only: :dev},
      {:earmark, "~> 1.2", only: :dev}
    ]
  end

  defp package do
    [
      maintainers: ["Safwan Kamarrudin"],
      licenses: ["MIT"],
      links: %{github: "https://github.com/safwank/Numerix"}
    ]
  end
end
