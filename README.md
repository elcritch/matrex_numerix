[![Build Status](https://travis-ci.org/safwank/Numerix.svg?branch=master)](https://travis-ci.org/safwank/Numerix)

# MatrexNumerix

This is a port of [Numerix](https://github.com/safwank/Numerix) to utilize the [Matrex](https://github.com/versilov/matrex/) library. This library is a collection of useful mathematical functions in Elixir with a slant towards statistics, linear algebra and machine learning. The unit tests have also been converted to standard ExUnit to remove reliance on a (defunct?) testing library. 

## Installation

Add `numerix` to your list of dependencies in `mix.exs`:

```elixir
  def deps do
    [{:matrex_numerix, "~> 0.6", github: "elcritch/matrex_numerix"}]
  end
```

Start `matrex_numerix` and its dependencies if you plan to use `gen_stage` or `flow`:

```elixir
  def application do
    [applications: [:matrex_numerix, :gen_stage, :flow]]
  end
```

## Examples

Check out the [tests](https://github.com/safwank/Numerix/tree/master/test) for examples.

## Documentation

Check out the [API reference](https://hexdocs.pm/numerix/api-reference.html) for the latest documentation.

## Features

### Tensor API

Numerix now includes a [Tensor API](https://hexdocs.pm/numerix/Numerix.Tensor.html) that lets you implement complex math functions with little code, similar to what you get from `numpy`. And since Numerix is written in Elixir, it uses `Flow` to run independent pieces of computation in parallel to speed things up. Depending on the type of calculations you're doing, the bigger the data and the more cores you have, the faster it gets.

NOTE: Parallelization can only get you so far. In terms of raw speed, a pure Elixir solution will always be much slower compared to one that leverages low-level routines like BLAS or similar.

### Statistics

* Mean
* Weighted mean
* Median
* Mode
* Range
* Variance
* Population variance
* Standard deviation
* Population standard deviation
* Moment
* Kurtosis
* Skewness
* Covariance
* Weighted covariance
* Population covariance
* Quantile
* Percentile

### Correlation functions

* Pearson
* Weighted Pearson

### Distance functions

* Mean squared error (MSE)
* Root mean square error (RMSE)
* Pearson
* Minkowski
* Euclidean
* Manhattan
* Jaccard

### General math functions

* nth root

### Special functions

* Logit
* Logistic

### Window functions

* Gaussian

### Linear algebra

* L1-norm
* L2-norm
* p-norm

### Data Processing

* FFT (DFT) for both real and complex values

### Linear regression

* Least squares best fit
* Prediction
* R-squared

### Kernel functions

* RBF

### Optimization

* ~~Genetic algorithms~~ (not ported yet)

### Neural network activation functions

* softmax
* softplus
* softsign
* sigmoid
* ReLU, leaky ReLU, ELU and SELU
* tanh
