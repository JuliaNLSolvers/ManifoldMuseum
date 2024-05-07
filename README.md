<div align="center">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="https://github.com/JuliaManifolds/Manifolds.jl/raw/master/docs/src/assets/logo-text-readme-dark.png">
      <img alt="Manifolds.jl logo with text on the side" src="https://github.com/JuliaManifolds/Manifolds.jl/raw/master/docs/src/assets/logo-text-readme.png">
    </picture>
</div>

[![docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliamanifolds.github.io/Manifolds.jl/stable/)
[![docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliamanifolds.github.io/Manifolds.jl/latest/)  [![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![CI](https://github.com/JuliaManifolds/Manifolds.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaManifolds/Manifolds.jl/actions?query=workflow%3ACI+branch%3Amaster)
[![codecov.io](http://codecov.io/github/JuliaManifolds/Manifolds.jl/coverage.svg?branch=master)](https://codecov.io/gh/JuliaManifolds/Manifolds.jl/)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

[![ACM TOMS](https://img.shields.io/badge/ACM%20TOMS-10.1145%2F3618296-blue.svg)](http://doi.org/10.1145/3618296)
[![DOI](https://zenodo.org/badge/190447542.svg)](https://zenodo.org/badge/latestdoi/190447542)

Package __Manifolds.jl__ aims to provide both a unified interface to define and
use manifolds as well as a library of manifolds to use for your projects.
This package is mostly stable, see https://github.com/JuliaManifolds/Manifolds.jl/issues/438 for planned upcoming changes.

## Getting started

To install the package just type

```julia
using Pkg; Pkg.add("Manifolds")
```

Then you can directly start, for example to stop half way from the north pole on the [`Sphere`](https://juliamanifolds.github.io/Manifolds.jl/stable/manifolds/sphere.html) to a point on the equator, you can generate the [`shortest_geodesic`](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/functions.html#ManifoldsBase.shortest_geodesic-Tuple{AbstractManifold,%20Any,%20Any}).
It internally employs [`exp`](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/functions.html#Base.exp-Tuple{AbstractManifold,%20Any,%20Any}) and [`log`](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/functions.html#Base.log-Tuple{AbstractManifold,%20Any,%20Any}).

```julia
using Manifolds
M = Sphere(2)
γ = shortest_geodesic(M, [0., 0., 1.], [0., 1., 0.])
γ(0.5)
```

## Ecosystem highlights

* A wide selection of Riemannian manifolds like spheres, hyperbolic spaces, Stiefel and Grassmann manifolds.
* Support for optimization on manifolds using [Manopt.jl](https://github.com/JuliaManifolds/Manopt.jl/).
* Support for many operations used in optimization and manifold-valued statistics with a focus on performance and ease of use.
* Connection manifolds.
* Lie groups.
* Atlases, charts and custom metrics (work in progress).
* A lightweight interface package: [ManifoldsBase.jl](https://github.com/JuliaManifolds/ManifoldsBase.jl).
* Differential equations on manifolds: [ManifoldDiffEq.jl](https://github.com/JuliaManifolds/ManifoldDiffEq.jl).
* Finite differences and automatic differentiation on manifolds using [ManifoldDiff.jl](https://github.com/JuliaManifolds/ManifoldDiff.jl) (work in progress).
* Integration and measures on manifolds: [ManifoldMeasures.jl](https://github.com/JuliaManifolds/ManifoldMeasures.jl) (work in progress).
* Functional manifolds: [FunManifolds.jl](https://github.com/JuliaManifolds/FunManifolds.jl) (work in progress).

## Support

If you have any questions regarding the Manifolds.jl ecosystem feel free to reach us using [Github discussion forums](https://github.com/JuliaManifolds/Manifolds.jl/discussions), [Julia Slack](https://julialang.org/slack/), [Julia Zulip](https://julialang.zulipchat.com/) or [Julia discourse](https://discourse.julialang.org/) forums. We are interested in new applications and methods on manifolds -- sharing your work is welcome!

## Citation

If you use `Manifolds.jl` in your work, please cite the following open access article

```biblatex
@article{AxenBaranBergmannRzecki:2023,
    AUTHOR    = {Axen, Seth D. and Baran, Mateusz and Bergmann, Ronny and Rzecki, Krzysztof},
    ARTICLENO = {33},
    DOI       = {10.1145/3618296},
    JOURNAL   = {ACM Transactions on Mathematical Software},
    MONTH     = {dec},
    NUMBER    = {4},
    TITLE     = {Manifolds.Jl: An Extensible Julia Framework for Data Analysis on Manifolds},
    VOLUME    = {49},
    YEAR      = {2023}
}
```

To refer to a certain version we recommend to also cite for example

```biblatex
@software{manifoldsjl-zenodo-mostrecent,
  Author = {Seth D. Axen and Mateusz Baran and Ronny Bergmann},
  Title = {Manifolds.jl},
  Doi = {10.5281/ZENODO.4292129},
  Url = {https://zenodo.org/record/4292129},
  Publisher = {Zenodo},
  Year = {2021},
  Copyright = {MIT License}
}
```

for the most recent version or a corresponding version specific DOI, see [the list of all versions](https://zenodo.org/search?page=1&size=20&q=conceptrecid:%224292129%22&sort=-version&all_versions=True).
Note that both citations are in [BibLaTeX](https://ctan.org/pkg/biblatex) format.
