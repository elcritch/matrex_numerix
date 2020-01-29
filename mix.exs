defmodule MatrexNumerix.Mixfile do
  use Mix.Project

  @version "0.6.0"

  def project do
    [
      app: :matrex_numerix,
      version: @version,
      elixir: "~> 1.7",
      start_permanent: Mix.env() == :prod,
      deps: deps(),
      package: package(),
      preferred_cli_env: ["bench.matrex": :bench, docs: :docs],
      description:
        "Port of Numerix to Matrex.",
      name: "MatrexNumerix",
      source_url: "https://github.com/elcritch/matrex_numerix",
      test_coverage: [tool: ExCoveralls],
      docs: docs()
    ]
  end

  def application do
    [applications: [:logger]]
  end

  def deps do
    [
      {:matrex, "~> 0.6.9", github: "elcritch/matrex", branch: "master"},
      {:ex_doc, "~> 0.18", only: :dev},
      {:earmark, "~> 1.2", only: :dev}
    ]
  end

  def package do
    [
      maintainers: ["Jaremy Creechley"],
      licenses: ["MIT"],
      links: %{github: "https://github.com/elcritch/matrex_numerix"}
    ]
  end

  def docs() do
    [
      main: "Matrex",
      logo: "docs/matrex_logo_dark_rounded.png",
      source_ref: "v#{@version}",
      canonical: "http://hexdocs.pm/matrex",
      extras: [
        "README.md"
      ]
    ]
  end
end
