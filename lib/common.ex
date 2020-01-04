defmodule MatrexNumerix.Common do
  @moduledoc """
  Common typespecs and functions.
  """

  @typedoc """
  A type representing an unreal number.
  """
  @type unreal_number :: :negative_infinity | :infinity

  @typedoc """
  A type representing the affinely extended real number system.
  """
  @type extended_number :: number | unreal_number

  @typedoc """
  Something that may be a float.
  """
  @type maybe_float :: float | nil

  @typedoc """
  A type representing a vector (1D Matrex) of numbers.
  """
  @type vector :: [number] | %Matrex{}

  @typedoc """
  Something that may be a vector.
  """
  @type maybe_vector :: vector | nil
end
