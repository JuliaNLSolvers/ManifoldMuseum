# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.18] –

### Changed

* Moved ambiguity testing to Aqua.jl

## [0.9.17] – 2024-04-23

### Added

* `Hyperrectangle` manifold with boundary.

## [0.9.16] – 2024-04-01

### Changed

* `NonlinearSolve.jl` and `PythonCall.jl` are no longer an upper bounded dependency (bugs were fixed).

## [0.9.15] – 2024-03-24

### Added

* using `DocumenterInterLinks` for links to other Julia packages documentation.
* Implementation of `sectional_curvature`, `sectional_curvature_min` and `sectional_curvature_max` for several manifolds.
* `sectional_curvature_matrix` function and a tutorial on coordinate-free curvature.

### Changed

* `default_vector_transport_method` for `GeneralUnitaryMatrices` other than `Rotations` was changed to `ProjectionTransport`.

### Fixed

* typographical errors in tutorials/working-in-charts.jl.
* several typographical errors in the docs
* unifies to use two backticks ``` `` ``` for math instead of ` $ ` further in the docs

## [0.9.14] – 2024-01-31

### Added

* `rand` on `UnitaryMatrices`
* `rand` on arbitrary `GroupManifold`s and manifolds with `IsGroupManifold` trait
  generating points and elements from the Lie algebra, respectively

## [0.9.13] – 2024-01-24

### Added

* added the real symplectic Grassmann manifold `SymplecticGrassmann`
* Introduce the manifold of `HamiltonianMatrices` and a wrapper for `Hamiltonian` matrices
* introduce `rand(:HamiltonianMatrices)`
* extend `rand` to also `rand!` for `HamiltonianMatrices`, `SymplecticMatrices` and `SymplecticStiefel`
* implement `riemannian_gradient` conversion for `SymplecticMatrices` and `SymplecticGrassmann`
* the new manifold of `MultinomialSymmetricPositiveDefinite` matrices
* `rand!` for `MultinomialDoublyStochastic` and `MultinomialSymmetric`

### Deprecated

* Rename `Symplectic` to `SimplecticMatrices` in order to have a `Symplectic` wrapper for such matrices as well in the future for the next breaking change.
* Rename `SymplecticMatrix` to `SymplecticElement` to clarify that it is the special matrix ``J_{2n}`` and not an arbitrary symplectic matrix.

### Fixed

* a bug that cause `project` for tangent vectors to return wrong results on `MultinomialDoublyStochastic`

## [0.9.12] – 2024-01-21

### Fixed

* Fixed `var` on `TranslationGroup`.

## [0.9.11] – 2023-12-27

### Fixed

* Fixed mixed array index number in-place `parallel_transport_to!` on zero-index `Euclidean`, on Julia 1.6.

## [0.9.10] – 2023-12-27

### Added

* Compatibility with `RecursiveArrayTools` v3.

### Fixed

* Fixed mixed array index number in-place `parallel_transport_to!` on real `Circle`, on Julia 1.6.
* Violations of MD004 lint rule in this file.

## [0.9.9] – 2023-12-25

### Fixed

* introduced a nonzero `atol` for all point and vector checks that compre to zero.
  This makes those checks a bit more relaxed by default and resolves [#630](https://github.com/JuliaManifolds/Manifolds.jl/issues/630).
* `default_estimation_method(M, f)` is deprecated, use `default_approximation_method(M, f)` for your specific method `f` on the manifold `M`.
* `AbstractEstimationMethod` is deprecated, use `AbstractApproximationMethod` instead.

## [0.9.8] - 2023-11-17

### Fixed

* Improved distribution of random vector generation for rotation matrices and complex circle.

## [0.9.7] – 2023-11-14

### Fixed

* Fixed `is_flat` for `CholeskySpace` and `SymmetricPositiveDefinite` with `LogCholeskyMetric` [https://github.com/JuliaManifolds/Manifolds.jl/issues/684](https://github.com/JuliaManifolds/Manifolds.jl/issues/684).

## [0.9.6] - 2023-11-09

### Fixed

* Fixed real coefficient basis for complex circle (an issue exposed by [https://github.com/JuliaManifolds/ManifoldsBase.jl/pull/173](https://github.com/JuliaManifolds/ManifoldsBase.jl/pull/173)).
* Fixed `VeeOrthogonalBasis` test for non-real manifolds.

## [0.9.5] - 2023-11-08

### Changed

* `identity_element` now returns a complex matrix for unitary group.
* `number_of_coordinates` is now exported.

## [0.9.4] - 2023-11-06

### Added

* Functions `inv_diff`, `inv_diff!`, `adjoint_inv_diff` and `adjoint_inv_diff!` that correspond to differentials and pullbacks of group inversion.
* Julia 1.10-rc CI workflow.

### Changed

* Documentation project files are marked as compatible with `BoundaryValueDiffEq` v5.

### Fixed

* Fixed issue with incorrect implementation of `apply_diff_group` in `GroupOperationAction` with left backward and right forward action [#669](https://github.com/JuliaManifolds/Manifolds.jl/issues/669).

## [0.9.3] - 2023-10-28

### Added

* Support for `BoundaryValueDiffEq` v5.

## [0.9.2] - 2023-10-27

### Added

* `rand(G; vector_at=Identity(G))` now works for translation, special orthogonal and special Euclidean groups `G` (issue [#665](https://github.com/JuliaManifolds/Manifolds.jl/issues/665)).
* `get_embedding` now works for `GeneralUnitaryMultiplicationGroup`.
* Github action that checks for NEWS.md changes.

## [0.9.1] - 2023-10-25

### Added

* a new retraction and its inverse for the fixed Rank Manifolds, the orthographic retraction.

## [0.9.0] - 2023-10-24

### Added

* Vector bundles are generalized to fiber bundles. Old `BundleFibers` functionality was reworked to better match mathematical abstractions. Fiber bundle functionality is experimental and minor changes may happen without a breaking release, with the exception of `TangentBundle` which is considered to be stable.
* `RotationTranslationAction` is introduced.

### Changed

* Sizes of all manifolds can now be either encoded in type or stored in a field to avoid over-specialization.
  The default is set to store the size in type parameter (except for `PowerManifold` and its variants), replicating the previous behavior.
  For field storage, pass the `parameter=:field` keyword argument to manifold constructor.
  For example statically sized `CenteredMatrices{m,n}` is now `CenteredMatrices{TypeParameter{Tuple{m,n}}}`, whereas the type of special Euclidean group with field-stored size is `CenteredMatrices{Tuple{Int,Int}}`. Similar change applies to:
  * `CenteredMatrices{m,n}`,
  * `CholeskySpace{N}`,
  * `Elliptope{N,K}`,
  * `Euclidean`,
  * `FixedRankMatrices{m,n,k}`,
  * `KendallsPreShapeSpace{n,k}`,
  * `KendallsShapeSpace{n,k}`,
  * `GeneralLinear{n}`,
  * `GeneralUnitaryMultiplicationGroup{n}`,
  * `GeneralizedGrassmann{n,k}`,
  * `GeneralizedStiefel{n,k}`,
  * `Grassmann{n,k}`,
  * `Heisenberg{n}`,
  * `Hyperbolic{n}`,
  * `MultinomialMatrices{N,M}`,
  * `MultinomialDoublyStochastic{n}`,
  * `MultinomialSymmetric{n}`,
  * `Orthogonal{n}`,
  * `PowerManifold`,
  * `PositiveArrays`,
  * `PositiveMatrices`,
  * `PositiveNumbers`,
  * `ProbabilitySimplex{n}`,
  * `SPDFixedDeterminant{n}`,
  * `SpecialLinear{n}`,
  * `SpecialOrthogonal{n}`,
  * `SpecialUnitary{n}`,
  * `SpecialEuclidean{n}`,
  * `SpecialEuclideanManifold{n}`,
  * `Spectrahedron{n,k}`,
  * `SphereSymmetricMatrices{N}`,
  * `Stiefel{n,k}`,
  * `SymmetricMatrices{N}`,
  * `SymmetricPositiveDefinite{n}`,
  * `SymmetricPositiveSemidefiniteFixedRank{n,k}`,
  * `Symplectic{n}`,
  * `SymplecticStiefel{n,k}`,
  * `TranslationGroup`,
  * `Tucker`.

  For example

  ```{julia}
  function Base.show(io::IO, ::CenteredMatrices{m,n}) where {m,n}
      return print(io, "CenteredMatrices($m, $n)")
  end
  ```

  needs to be replaced with

  ```{julia}
  function Base.show(io::IO, ::CenteredMatrices{TypeParameter{Tuple{m,n}}}) where {m,n}
      return print(io, "CenteredMatrices($m, $n)")
  end
  ```

  for statically-sized groups and

  ```{julia}
  function Base.show(io::IO, M::CenteredMatrices{Tuple{Int,Int}})
      m, n = get_parameter(M.size)
      return print(io, "CenteredMatrices($m, $n; parameter=:field)")
  end
  ```

  for groups with size stored in field. Alternatively, you can use a single generic method like this:

  ```{julia}
  function Base.show(io::IO, M::CenteredMatrices{T}) where {T}
      m, n = get_parameter(M)
      if T <: TypeParameter
          return print(io, "CenteredMatrices($m, $n)")
      else
          return print(io, "CenteredMatrices($m, $n; parameter=:field)")
      end
  end
  ```

* Argument order for type aliases `RotationActionOnVector` and `RotationTranslationActionOnVector`: most often dispatched on argument is now first.
* A more consistent handling of action direction was introduced. 4-valued `ActionDirection` was split into 2-valued `ActionDirection` (either left or right action) and `GroupActionSide` (action acting from the left or right side). See [https://github.com/JuliaManifolds/Manifolds.jl/issues/637](https://github.com/JuliaManifolds/Manifolds.jl/issues/637) for a design discussion.

### Removed

* `ProductRepr` is removed; please use `ArrayPartition` instead.
* Default methods throwing "not implemented" `ErrorException` for some group-related operations. Standard `MethodError` is now thrown instead.
* `LinearAffineMetric` was deprecated in a previous release and the symbol is now removed.
  Please use `AffineInvariantMetric` instead.
